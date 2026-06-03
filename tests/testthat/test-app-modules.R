skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

test_that("intro_server writes form inputs into the shared meta store", {
  meta <- shiny::reactiveValues()
  shiny::testServer(intro_server, args = list(meta = meta), {
    session$setInputs(expID = "E1", species = "Daphnia", endpoint = "reproduction",
                      chem1 = "CPF", chem2 = "MPs", unit = "mg/L")
    expect_equal(meta$expID, "E1")
    expect_equal(meta$chem1, "CPF")
    expect_equal(meta$chem2, "MPs")
    expect_equal(meta$unit, "mg/L")
    expect_equal(meta$endpoint, "reproduction")
  })
})

test_that("intro_ui builds a Shiny UI fragment", {
  ui <- intro_ui("intro")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})

test_that("single_server surfaces validation errors and withholds a fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(Conc = c(0, 1, 2), Res = c(100, 50, 10)),
                     path, row.names = FALSE)   # only 3 distinct concentrations
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "short.csv"))
    expect_gt(length(errs()), 0)
    expect_match(paste(errs(), collapse = " "), "distinct")
  })
})

test_that("single_ui builds a Shiny UI fragment", {
  expect_true(inherits(single_ui("single"), c("shiny.tag", "shiny.tag.list", "bslib_fragment")))
})

test_that("binary_server: compare-all loop fills the store live and selects a model", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30, model = "reference",
                      file = list(datapath = csv, name = "binary.csv"))
    expect_length(errs(), 0)

    # supply both single curves deterministically via Simulate, distinct max
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)
    expect_true(isTRUE(frozen()))
    expect_null(current_fit())

    # Kick the staged compare-all loop. The first model fits when the queue is
    # armed; the rest advance on timer ticks.
    session$setInputs(compare_all = 1)
    for (i in 1:8) session$elapse(5)             # drain the stepper queue

    # all four staged fits stored, with staged df counts (curves fixed, only the
    # interaction params are free: reference 0, SA 1 (a), DR/DL 2 (a, b))
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    expect_equal(fits_store()[["reference"]]$df, 0L)
    expect_equal(fits_store()[["DR"]]$df, 2L)
    # adding the interaction term cannot worsen the fit (a=0 reproduces the parent)
    expect_lte(fits_store()[["SA"]]$objective,
               fits_store()[["reference"]]$objective + 1e-6)

    # comparison + chosen produced once the chain completed
    expect_setequal(names(last_compare()$fits), c("reference", "SA", "DR", "DL"))
    expect_equal(nrow(last_compare()$comparison), 3)
    chosen <- last_compare()$chosen
    expect_true(chosen %in% c("reference", "SA", "DR", "DL"))

    # picker follows the chosen model (testServer doesn't echo updateSelectInput)
    session$setInputs(model = chosen)
    expect_equal(current_fit()$deviation, chosen)

    # "explore by hand" still works: Autofit the displayed model (staged a/b)
    session$setInputs(model = "SA", autofit = 1)
    expect_equal(current_fit()$deviation, "SA")
    expect_true("a" %in% names(current_fit()$par))

    # editing a single curve keeps the workspace and re-evaluates stored fits
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_true(isTRUE(frozen()))
    expect_gt(length(fits_store()), 0)
    expect_null(last_compare())                  # comparison invalidated by curve change
  })
})

test_that("a bundled example dataset ships and validates as binary continuous", {
  ex <- system.file("extdata", "binary_ia_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled example not installed")
  df <- read_upload(ex)
  expect_length(validate_upload(df, "binary", "continuous"), 0)
  expect_true(all(c("C1", "C2", "Res") %in% names(df)))
})

test_that("binary_server falls back to the bundled example before any upload", {
  skip_on_cran()
  ex <- system.file("extdata", "binary_ia_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled example not installed")
  meta <- shiny::reactiveValues(chem1 = "CPF", chem2 = "IMI")
  shiny::testServer(binary_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_length(errs(), 0)                 # no upload, yet valid
    expect_gt(nrow(engine_df()), 0)          # data is available
  })
})

test_that("binary_ui: comparison table + tune/refine panel with joint optimize, then diagnostics", {
  html <- as.character(binary_ui("binary"))
  expect_false(grepl("Freeze curves", html, fixed = TRUE))
  # Stage 2: staged verdict hero (compare-all + alpha + table)
  expect_match(html, "binary-compare_all", fixed = TRUE)
  expect_match(html, "Fit &amp; compare all models", fixed = TRUE)
  expect_match(html, "binary-results", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "significance threshold for the model comparison", fixed = TRUE)
  # Tune / refine panel (under the table): picker, help, a/b grid, three actions
  expect_match(html, "Tune / refine selected model", fixed = TRUE)
  expect_match(html, "binary-model", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
  expect_match(html, "binary-val_a", fixed = TRUE)
  expect_match(html, "Autofit (a, b)", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "binary-optimize_all", fixed = TRUE)          # NEW joint button
  expect_match(html, "Optimize all params (joint)", fixed = TRUE)
  expect_match(html, "binary-refine_readout", fixed = TRUE)
  expect_match(html, "binary-refined_badge", fixed = TRUE)
  # Stage 3: diagnostics
  expect_match(html, "binary-surface", fixed = TRUE)
  expect_match(html, "binary-isobole", fixed = TRUE)
  expect_match(html, "binary-op", fixed = TRUE)
  expect_match(html, "binary-cis", fixed = TRUE)
  expect_match(html, "binary-n_starts", fixed = TRUE)
  # the old staged Find-best button is still gone
  expect_false(grepl("binary-find_best", html, fixed = TRUE))
})

test_that("single_ui shows the model equation and Autofit/Simulate buttons", {
  html <- as.character(single_ui("single"))
  expect_match(html, "Y = max / (1 + (C / EC50)", fixed = TRUE)
  expect_match(html, "Autofit parameters", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "the curve is fit to the proportion", fixed = TRUE)
})

test_that("interaction_param_row renders a value input with inert bound cells", {
  ns <- shiny::NS("binary")
  html <- as.character(interaction_param_row(ns, "a", "a", "overall strength"))
  expect_match(html, "binary-val_a", fixed = TRUE)   # editable Value input present
  expect_match(html, "—", fixed = TRUE)          # em-dash placeholder present
  expect_false(grepl("binary-lo_a", html, fixed = TRUE))  # no Lower input
  expect_false(grepl("binary-hi_a", html, fixed = TRUE))  # no Upper input
})

test_that("binary_server: joint refine polishes a model but never moves the staged verdict", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30, model = "reference",
                      file = list(datapath = csv, name = "binary.csv"))
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    session$setInputs(compare_all = 1)
    for (i in 1:8) session$elapse(5)
    chosen     <- last_compare()$chosen
    staged_cmp <- last_compare()$comparison
    staged_obj <- fits_store()[[chosen]]$objective
    staged_df  <- fits_store()[[chosen]]$df

    session$setInputs(model = chosen)
    session$setInputs(optimize_all = 1)

    rf <- refined_fits()[[chosen]]
    expect_true(isTRUE(rf$joint))
    expect_lte(rf$objective, staged_obj + 1e-6)
    expect_gt(rf$df, staged_df)

    expect_identical(last_compare()$chosen, chosen)
    expect_equal(last_compare()$comparison, staged_cmp)
    expect_equal(fits_store()[[chosen]]$objective, staged_obj)
    expect_equal(fits_store()[[chosen]]$df, staged_df)
    expect_true(isTRUE(display_fit()$joint))

    # the before->after readout is action-scoped: switching to another model
    # clears it, so it can never show one model's numbers while another is shown.
    # (The refined fit itself persists -- only the readout is transient.)
    other <- setdiff(c("reference", "SA", "DR", "DL"), chosen)[[1]]
    session$setInputs(model = other)
    expect_null(refine_post())                    # readout hidden for the un-refined model
    session$setInputs(model = chosen)             # back to the refined model
    expect_false(is.null(refined_fits()[[chosen]]))  # its refined fit survived the switch

    session$setInputs(autofit = 1)
    expect_null(refined_fits()[[chosen]])
    expect_false(isTRUE(display_fit()$joint))

    session$setInputs(model = chosen, optimize_all = 2)
    expect_false(is.null(refined_fits()[[chosen]]))
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_equal(length(refined_fits()), 0)
  })
})

