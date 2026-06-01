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

test_that("single_server validates, fits, and exposes a single fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    # Valid continuous single-chemical CSV with >= 4 distinct concentrations.
    path <- tempfile(fileext = ".csv")
    utils::write.csv(
      data.frame(Conc = c(0, 0.1, 0.3, 1, 3, 10),
                 Res  = ll3_predict(c(0, 0.1, 0.3, 1, 3, 10), 100, 2, 0.5)),
      path, row.names = FALSE)
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "single.csv"))
    expect_length(errs(), 0)               # passes validation
    session$setInputs(fit = 1)
    expect_equal(fit_r()$kind, "single")
    expect_equal(unname(round(fit_r()$par[["ec50"]], 2)), 0.5)
  })
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
