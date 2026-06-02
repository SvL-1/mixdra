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

test_that("binary_server stages the fit: split, freeze gate, average max, override, invalidate", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      file = list(datapath = csv, name = "binary.csv"))
    expect_length(errs(), 0)

    # marginal split: each panel's series drops the other chemical's column and
    # exposes the kept chemical's concentration as C1
    expect_false("C2" %in% names(m1()))
    expect_false("C2" %in% names(m2()))
    expect_true("C1" %in% names(m1()))
    expect_true("C1" %in% names(m2()))

    # nothing frozen yet
    expect_false(isTRUE(frozen()))

    # supply both single curves deterministically via Simulate, with distinct max
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    # frozen curve parameters: shared max is the average; per-chemical slope/ec50 kept
    expect_equal(unname(curve_params()[["max"]]), 650)
    expect_equal(unname(curve_params()[["slope1"]]), 2)
    expect_equal(unname(curve_params()[["ec502"]]), 5)

    # freeze -> fits the four models and reveals the result
    session$setInputs(freeze = 1)
    expect_true(isTRUE(frozen()))
    expect_setequal(names(res_r()$fits), c("reference", "SA", "DR", "DL"))
    expect_equal(shown_fit()$deviation, res_r()$chosen)

    # picker overrides the displayed model
    session$setInputs(model = "DR")
    expect_equal(shown_fit()$deviation, "DR")

    # the per-model explanation follows the picker. With DR shown, it is the
    # DR block; the DL block label is absent.
    expect_match(output$interaction_help$html, "Dose-ratio dependent", fixed = TRUE)
    expect_false(grepl("Dose-level dependent", output$interaction_help$html, fixed = TRUE))
    # switching the picker to DL switches the block; with reference CA the DL
    # formula is the &Sigma;TU variant (distinct from the intro's TU/&Sigma;TU).
    session$setInputs(model = "DL")
    expect_equal(shown_fit()$deviation, "DL")
    expect_match(output$interaction_help$html, "Dose-level dependent", fixed = TRUE)
    expect_match(output$interaction_help$html, "b&middot;&Sigma;TU)", fixed = TRUE)

    # editing a single curve after freezing invalidates the freeze
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_false(isTRUE(frozen()))
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
