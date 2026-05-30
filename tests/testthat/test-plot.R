test_that("plot_dose_response builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "mixture", conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 2, slope2 = 1, ec501 = 0.5, ec502 = 0.3))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), C2 = 0,
                   Affected = c(10, 6, 3, 1), Exposed = rep(10, 4))
  p <- plot_dose_response(fit, df, chem = 1)
  expect_s3_class(p, "plotly")
})
