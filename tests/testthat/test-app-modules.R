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

test_that("binary_server: table auto-fills, row-click selects, joint refine is isolated", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      file = list(datapath = csv, name = "binary.csv"))
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    # auto-fill: no button, just drain the timer
    for (i in 1:8) session$elapse(5)
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    chosen     <- last_compare()$chosen
    staged_cmp <- last_compare()$comparison

    # with no row clicked, the displayed model defaults to the winner
    expect_equal(selected_model(), chosen)
    expect_equal(current_fit()$deviation, chosen)

    # clicking a row (index in reference,SA,DR,DL order) selects that model
    ord <- intersect(c("reference", "SA", "DR", "DL"), names(fits_store()))
    session$setInputs(results_rows_selected = match("SA", ord))
    expect_equal(selected_model(), "SA")
    expect_equal(current_fit()$deviation, "SA")

    # optimise the selected (SA) model jointly -> refined_fits gets SA, verdict frozen
    staged_sa_obj <- fits_store()[["SA"]]$objective
    session$setInputs(optimize_all = 1)
    expect_true(isTRUE(refined_fits()[["SA"]]$joint))
    expect_lte(refined_fits()[["SA"]]$objective, staged_sa_obj + 1e-6)
    expect_true(isTRUE(display_fit()$joint))
    expect_identical(last_compare()$chosen, chosen)
    expect_equal(last_compare()$comparison, staged_cmp)
    expect_equal(fits_store()[["SA"]]$objective, staged_sa_obj)

    # readout is action-scoped: selecting a different row clears the before->after
    other <- setdiff(ord, "SA")[[1]]
    session$setInputs(results_rows_selected = match(other, ord))
    expect_null(refine_post())

    # a curve edit re-runs the staged loop and clears refined fits
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    for (i in 1:8) session$elapse(5)
    expect_false(is.null(last_compare()))
    expect_equal(length(refined_fits()), 0)
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
  expect_match(html, "binary-val_b", fixed = TRUE)
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


