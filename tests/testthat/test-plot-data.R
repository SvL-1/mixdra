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

test_that("obs_pred_data returns fit$pred for a mixture fit", {
  fit <- list(kind = "mixture", pred = c(0.9, 0.5, 0.2))
  df <- data.frame(C1 = c(0, 1, 2), C2 = 0, Affected = c(9, 5, 2), Exposed = rep(10, 3))
  d <- obs_pred_data(fit, df)
  expect_equal(d$predicted, c(0.9, 0.5, 0.2))
  expect_equal(d$observed, c(0.9, 0.5, 0.2))
})

test_that("obs_pred_data computes predictions for a single fit", {
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0.5, 1, 2), Res = c(8, 5, 2))
  d <- obs_pred_data(fit, df)
  expect_equal(d$predicted, ll3_predict(df$C1, 10, 2, 1))
  expect_equal(d$observed, df$Res)
})

test_that("surface_grid_data has correct dims, orientation, and control corner", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 0.9, slope1 = 1.5, slope2 = 1.2, ec501 = 0.1, ec502 = 0.2))
  df <- data.frame(C1 = c(0, 0.1, 0.3), C2 = c(0, 0.2, 0.4),
                   Affected = c(9, 5, 2), Exposed = rep(10, 3))
  g <- surface_grid_data(fit, df, n = 5)
  expect_equal(dim(g$z), c(5, 5))                 # nrow = length(y), ncol = length(x)
  i <- 3; j <- 4                                  # orientation: z[i,j] == f(x[j], y[i])
  expect_equal(g$z[i, j], predict_grid(fit, g$x_vals[j], g$y_vals[i]))
  expect_equal(g$z[1, 1], 0.9)                    # control corner (0,0) == max
})

test_that("surface_grid_data errors on a non-binary fit", {
  fit <- list(n_chem = 1, conc_cols = "C1")
  df <- data.frame(C1 = 1:3, Res = 1:3)
  expect_error(surface_grid_data(fit, df), "binary")
})

test_that("surface_grid_data does not error when some cells are NaN", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 0.5, slope1 = 5, slope2 = 5, ec501 = 0.01, ec502 = 0.01))
  df <- data.frame(C1 = c(0, 1), C2 = c(0, 1), Affected = c(5, 0), Exposed = c(10, 10))
  expect_silent(g <- surface_grid_data(fit, df, n = 10))
  expect_equal(dim(g$z), c(10, 10))
})
