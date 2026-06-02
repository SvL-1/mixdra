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

test_that("binary_server workspace: freeze gate, autofit/simulate/find-best, invalidate", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi", "binary_ca_mps_cpf_imi_continuous.csv")
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
    expect_equal(unname(curve_params()[["max"]]), 650)

    # freeze locks curves but does NOT fit
    session$setInputs(freeze = 1)
    expect_true(isTRUE(frozen()))
    expect_null(current_fit())

    # Autofit the selected (S/A) model -> a fit for SA is stored
    session$setInputs(model = "SA")
    session$setInputs(autofit = 1)
    expect_equal(current_fit()$deviation, "SA")
    expect_true("a" %in% names(current_fit()$par))

    # Simulate S/A with an entered a -> a simulated fit replaces the stored one
    session$setInputs(val_a = 3, simulate = 1)
    expect_true(isTRUE(current_fit()$simulated))
    expect_equal(current_fit()$par[["a"]], 3)

    # Find best model -> comparison of all four + a chosen model
    session$setInputs(find_best = 1)
    expect_setequal(names(last_compare()$fits), c("reference", "SA", "DR", "DL"))
    expect_equal(nrow(last_compare()$comparison), 3)
    expect_true(last_compare()$chosen %in% c("reference", "SA", "DR", "DL"))

    # in a live session Find best switches the picker to the chosen model; testServer
    # does not echo updateSelectInput, so set it explicitly and confirm current_fit
    # follows the picker to that stored (fitted) model -- the per-model store design.
    session$setInputs(model = last_compare()$chosen)
    expect_equal(current_fit()$deviation, last_compare()$chosen)
    expect_false(isTRUE(current_fit()$simulated))

    # Optimize-all: jointly refine the displayed model; objective must not worsen.
    chosen <- last_compare()$chosen
    pre_obj <- current_fit()$objective
    session$setInputs(optimize_all = 1)
    expect_true(isTRUE(current_fit()$joint))
    expect_equal(current_fit()$deviation, chosen)
    expect_lte(current_fit()$objective, pre_obj + 1e-6)

    # Pinning max via equal bounds holds it across a re-optimize.
    held <- unname(current_fit()$par[["max"]])
    session$setInputs(olo_max = held, ohi_max = held, optimize_all = 2)
    expect_equal(unname(current_fit()$par[["max"]]), held, tolerance = 1e-8)
    expect_true("max" %in% current_fit()$fixed)

    # CI table: a pinned parameter shows no interval (it was not estimated),
    # but its held estimate is preserved. Uses the same helper the renderer does.
    cis <- blank_pinned_ci(
      param_ci(current_fit(), engine_df(), current_fit()$reference,
               current_fit()$deviation, current_fit()$response),
      current_fit()$fixed)
    expect_true(is.na(cis$lower[cis$parameter == "max"]))
    expect_true(is.na(cis$upper[cis$parameter == "max"]))
    expect_false(is.na(cis$estimate[cis$parameter == "max"]))  # held value preserved

    # editing a single curve after freezing invalidates everything
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_false(isTRUE(frozen()))
    expect_equal(length(fits_store()), 0)
    expect_null(last_compare())
  })
})

test_that("binary_ui: Freeze-only button + Stage 2 workspace", {
  html <- as.character(binary_ui("binary"))
  # Freeze no longer fits
  expect_match(html, "Freeze curves", fixed = TRUE)
  expect_false(grepl("Freeze curves → fit interactions", html, fixed = TRUE))
  # Stage 2 workspace: model picker, a/b value input, the three actions
  expect_match(html, "binary-model", fixed = TRUE)
  expect_match(html, "binary-val_a", fixed = TRUE)
  expect_match(html, "Autofit (a, b)", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "Find best model", fixed = TRUE)
  # relocated fit options + the per-model explanation slot
  expect_match(html, "binary-n_starts", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
  # alpha now sits next to Find best model and carries an explanation
  expect_match(html, "significance threshold for the model comparison", fixed = TRUE)
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

test_that("binary_ui: Optimize-all card with the 7-param bounds grid", {
  html <- as.character(binary_ui("binary"))
  expect_match(html, "Optimize all params", fixed = TRUE)     # button
  expect_match(html, "binary-olo_max", fixed = TRUE)          # lower input
  expect_match(html, "binary-ohi_ec502", fixed = TRUE)        # upper input
  expect_match(html, "binary-oval_a", fixed = TRUE)           # start value input
  expect_match(html, "binary-optimize_all", fixed = TRUE)     # action id
})
