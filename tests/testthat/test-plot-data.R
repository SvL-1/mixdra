test_that("dr_curve_data curve matches ll3_predict for the chosen chemical", {
  fit <- list(kind = "mixture", conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 2, slope2 = 1, ec501 = 0.5, ec502 = 0.3))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), C2 = c(0, 0, 0, 0),
                   Affected = c(10, 6, 3, 1), Exposed = rep(10, 4))
  d <- dr_curve_data(fit, df, chem = 1)
  expect_equal(d$curve$response, ll3_predict(d$curve$conc, 1, 2, 0.5))
  expect_equal(d$observed$response, df$Affected / df$Exposed)
  expect_equal(d$chem, "C1")
})

test_that("dr_curve_data handles a single-chemical fit", {
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0.5, 1, 2, 4), Res = c(8, 5, 2, 1))
  d <- dr_curve_data(fit, df)
  expect_equal(d$curve$response, ll3_predict(d$curve$conc, 10, 2, 1))
  expect_equal(d$observed$response, df$Res)
})
