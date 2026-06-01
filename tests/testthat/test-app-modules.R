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

test_that("binary_server fits the four models and exposes the chosen model", {
  skip_on_cran()
  meta <- shiny::reactiveValues()
  shiny::testServer(binary_server, args = list(meta = meta), {
    # Use the validated continuous binary fixture (C1, C2, Res).
    csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA",
                      thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      lo_max = NA, hi_max = NA, lo_slope1 = NA, hi_slope1 = NA,
                      lo_slope2 = NA, hi_slope2 = NA, lo_ec501 = NA, hi_ec501 = NA,
                      lo_ec502 = NA, hi_ec502 = NA,
                      file = list(datapath = csv, name = "binary.csv"))
    expect_length(errs(), 0)
    session$setInputs(fit = 1)
    res <- res_r()
    expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
    expect_true(res$chosen %in% c("reference", "SA", "DR", "DL"))
    # the displayed fit defaults to the chosen model (deviation tag matches)
    expect_equal(shown_fit()$deviation, res$chosen)
  })
})

test_that("binary_ui builds a Shiny UI fragment", {
  expect_true(inherits(binary_ui("binary"), c("shiny.tag", "shiny.tag.list", "bslib_fragment")))
})

test_that("single_ui shows the model equation and Autofit/Simulate buttons", {
  html <- as.character(single_ui("single"))
  expect_match(html, "Y = max / (1 + (C / EC50)", fixed = TRUE)
  expect_match(html, "Autofit parameters", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
})
