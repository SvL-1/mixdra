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

test_that("binary_server: Fit-interactions fills the staged block, Optimize-all appends an isolated joint block", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      file = list(datapath = csv, name = "binary.csv"))
    # the comparison table renders from the very start, before any curve/interaction
    # fit (all four rows present, empty) -- it is no longer gated on `frozen`.
    expect_false(is.null(output$results))
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    # "Fit interaction models" runs the whole staged chain in one (blocking) go.
    session$setInputs(fit_interactions = 1)
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    chosen     <- last_compare()$chosen
    staged_cmp <- last_compare()$comparison

    # with no row clicked, the displayed model defaults to the winner
    expect_equal(selected_model(), chosen)
    expect_equal(current_fit()$deviation, chosen)

    # only the staged block exists so far: 4 rows, all variant "staged"
    expect_equal(length(table_rows()), 4L)
    expect_true(all(vapply(table_rows(), function(e) e$variant == "staged", logical(1))))

    # clicking a staged row (index in reference,SA,DR,DL order) selects that model
    ord <- intersect(c("reference", "SA", "DR", "DL"), names(fits_store()))
    session$setInputs(results_rows_selected = match("SA", ord))
    expect_equal(selected_model(), "SA")
    expect_equal(current_fit()$deviation, "SA")
    expect_false(isTRUE(display_fit()$joint))      # a staged row shows the staged fit

    # "Optimize all params" jointly refines EVERY model -> a separator + joint block
    staged_sa_obj <- fits_store()[["SA"]]$objective
    session$setInputs(optimize_all = 1)
    expect_setequal(names(refined_fits()), c("reference", "SA", "DR", "DL"))
    expect_true(all(vapply(refined_fits(), function(f) isTRUE(f$joint), logical(1))))
    expect_lte(refined_fits()[["SA"]]$objective, staged_sa_obj + 1e-6)

    # the table now has the staged block, a separator, and the joint block
    expect_equal(length(table_rows()), 4L + 1L + 4L)
    sep <- vapply(table_rows(), function(e) isTRUE(e$separator), logical(1))
    expect_equal(sum(sep), 1L)

    # the staged verdict is untouched (still built from fits_store only)
    expect_identical(last_compare()$chosen, chosen)
    expect_equal(last_compare()$comparison, staged_cmp)
    expect_equal(fits_store()[["SA"]]$objective, staged_sa_obj)

    # clicking the JOINT SA row shows the joint-refined fit in the diagnostics
    sa_joint <- which(vapply(table_rows(), function(e)
      identical(e$model, "SA") && identical(e$variant, "joint"), logical(1)))
    session$setInputs(results_rows_selected = sa_joint)
    expect_equal(selected_model(), "SA")
    expect_true(isTRUE(display_fit()$joint))

    # a curve edit clears BOTH the staged and joint fits -- the user re-runs the
    # buttons (curve columns still show the new values, read live in the renderer).
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_null(last_compare())
    expect_equal(length(fits_store()), 0)
    expect_equal(length(refined_fits()), 0)
  })
})

test_that("a bundled example dataset ships and validates as binary continuous", {
  ex <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled example not installed")
  df <- read_upload(ex)
  expect_length(validate_upload(df, "binary", "continuous"), 0)
  expect_true(all(c("C1", "C2", "Res") %in% names(df)))
})

test_that("binary_server falls back to the bundled example before any upload", {
  skip_on_cran()
  ex <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled example not installed")
  meta <- shiny::reactiveValues(chem1 = "CPF", chem2 = "IMI")
  shiny::testServer(binary_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_length(errs(), 0)                 # no upload, yet valid
    expect_gt(nrow(engine_df()), 0)          # data is available
  })
})

test_that("binary_ui: always-visible table + Fit/Optimize buttons, then diagnostics", {
  html <- as.character(binary_ui("binary"))
  expect_false(grepl("Freeze curves", html, fixed = TRUE))
  # Stage 2: the table + alpha + Fit-interactions + the single Optimize button
  expect_match(html, "binary-results", fixed = TRUE)
  expect_match(html, "binary-fit_interactions", fixed = TRUE)
  expect_match(html, "Fit interaction models", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "threshold for the model comparison", fixed = TRUE)
  expect_match(html, "binary-optimize_all", fixed = TRUE)
  expect_match(html, "Optimize all params (joint)", fixed = TRUE)
  expect_match(html, "binary-refine_readout", fixed = TRUE)
  expect_match(html, "binary-refined_badge", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
  expect_match(html, "binary-objective", fixed = TRUE)
  # Stage 3: diagnostics
  expect_match(html, "binary-surface", fixed = TRUE)
  expect_match(html, "binary-isobole", fixed = TRUE)
  expect_match(html, "binary-op", fixed = TRUE)
  expect_match(html, "binary-cis", fixed = TRUE)
  expect_match(html, "binary-n_starts", fixed = TRUE)
  # the manual path + compare-all button + model picker are gone
  expect_false(grepl("binary-compare_all", html, fixed = TRUE))
  expect_false(grepl("binary-model", html, fixed = TRUE))
  expect_false(grepl("binary-autofit", html, fixed = TRUE))
  expect_false(grepl("binary-simulate", html, fixed = TRUE))
  expect_false(grepl("binary-val_a", html, fixed = TRUE))
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

test_that("a bundled example dataset ships and validates as ternary continuous", {
  ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled ternary example not installed")
  df <- read_upload(ex)
  expect_length(validate_upload(df, "ternary", "continuous"), 0)
  expect_true(all(c("C1", "C2", "C3", "Res") %in% names(df)))
})

test_that("intro_server writes chem3 into the shared meta store", {
  meta <- shiny::reactiveValues()
  shiny::testServer(intro_server, args = list(meta = meta), {
    session$setInputs(chem1 = "CPF", chem2 = "FBSA", chem3 = "IMI")
    expect_equal(meta$chem3, "IMI")
  })
})

test_that("ternary_ui builds the sidebar, three curve panels, and the three stages", {
  html <- as.character(ternary_ui("ternary"))
  # sidebar controls + disabled-radio script
  expect_match(html, "ternary-response", fixed = TRUE)
  expect_match(html, "ternary-reference", fixed = TRUE)
  expect_match(html, "coming soon", fixed = TRUE)
  expect_match(html, "prop('disabled', true)", fixed = TRUE)
  # Stage 1: three curve panels
  expect_match(html, "ternary-chem1-autofit", fixed = TRUE)
  expect_match(html, "ternary-chem2-autofit", fixed = TRUE)
  expect_match(html, "ternary-chem3-autofit", fixed = TRUE)
  # Stage 2: fit button + hub table
  expect_match(html, "ternary-fit_asa", fixed = TRUE)
  expect_match(html, "ternary-hub", fixed = TRUE)
  # Stage 3: plots + effect readout
  expect_match(html, "ternary-isoplane", fixed = TRUE)
  expect_match(html, "ternary-sigma_tu", fixed = TRUE)
  expect_match(html, "ternary-effect", fixed = TRUE)
  # no alpha / joint controls (deliberately dropped vs binary)
  expect_false(grepl("ternary-alpha", html, fixed = TRUE))
  expect_false(grepl("ternary-optimize_all", html, fixed = TRUE))
})
