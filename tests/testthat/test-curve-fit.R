skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

# A clean continuous single-chemical frame with a known curve (ec50 = 0.5).
.curve_df <- function() {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  data.frame(C1 = conc, Res = ll3_predict(conc, 100, 2, 0.5))
}

test_that("curve_fit_server Autofit fits and exposes the fit via the returned reactive", {
  meta <- shiny::reactiveValues()
  df <- .curve_df()
  shiny::testServer(curve_fit_server,
                    args = list(fit_df = shiny::reactive(df), meta = meta), {
    session$setInputs(val_max = NA, val_slope = NA, val_ec50 = NA,
                      lo_max = NA, hi_max = NA, lo_slope = NA, hi_slope = NA,
                      lo_ec50 = NA, hi_ec50 = NA, autofit = 1)
    expect_equal(current_fit()$kind, "single")
    expect_equal(unname(round(current_fit()$par[["ec50"]], 2)), 0.5)
  })
})

test_that("curve_fit_server Autofit blocks when a lower bound exceeds its upper", {
  meta <- shiny::reactiveValues()
  df <- .curve_df()
  shiny::testServer(curve_fit_server,
                    args = list(fit_df = shiny::reactive(df), meta = meta), {
    session$setInputs(val_max = NA, val_slope = NA, val_ec50 = NA,
                      lo_max = NA, hi_max = NA, lo_slope = 5, hi_slope = 1,
                      lo_ec50 = NA, hi_ec50 = NA, autofit = 1)
    expect_null(current_fit())
  })
})

test_that("curve_fit_server Simulate evaluates the typed parameter values", {
  meta <- shiny::reactiveValues()
  df <- .curve_df()
  shiny::testServer(curve_fit_server,
                    args = list(fit_df = shiny::reactive(df), meta = meta), {
    session$setInputs(val_max = 100, val_slope = 2, val_ec50 = 0.5, simulate = 1)
    expect_equal(unname(current_fit()$par[c("max", "slope", "ec50")]), c(100, 2, 0.5))
    expect_equal(current_fit()$ssr, 0)
  })
})

test_that("curve_fit_ui shows the model equation and Autofit/Simulate buttons", {
  html <- as.character(curve_fit_ui("curve"))
  expect_match(html, "Y = max / (1 + (C / EC50)", fixed = TRUE)
  expect_match(html, "Autofit parameters", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "model prediction", fixed = TRUE)
})
