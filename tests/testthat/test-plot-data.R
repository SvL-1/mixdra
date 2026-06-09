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

test_that("dr_curve_data with no fit returns raw points only (no curve)", {
  df <- data.frame(C1 = c(0, 0.5, 1, 2), Res = c(10, 8, 5, 2))
  d <- dr_curve_data(NULL, df)
  expect_null(d$curve)
  expect_equal(d$observed$conc, df$C1)
  expect_equal(d$observed$response, df$Res)
  # also works for quantal data (proportion)
  q <- data.frame(C1 = c(0, 1, 2), Affected = c(0, 5, 9), Exposed = rep(10, 3))
  expect_equal(dr_curve_data(NULL, q)$observed$response, q$Affected / q$Exposed)
})

test_that("dr_curve_data offsets controls to min(positive)/10 for the log axis", {
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), Res = c(10, 8, 5, 2))
  d <- dr_curve_data(fit, df)
  expect_equal(d$control_x, 0.5 / 10)        # min positive conc is 0.5
  # the fitted curve reaches down to the control marker, not just the lowest dose
  expect_equal(min(d$curve$conc), 0.5 / 10)
  # observed concentrations are untouched (0 preserved for the linear axis)
  expect_equal(d$observed$conc, df$C1)
})

test_that("dr_curve_data sets control_x to NA when there is no control row", {
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0.5, 1, 2, 4), Res = c(8, 5, 2, 1))
  d <- dr_curve_data(fit, df)
  expect_true(is.na(d$control_x))
  expect_equal(min(d$curve$conc), 0.5)
})

test_that("dr_curve_data reports control_x with no fit (raw points)", {
  df <- data.frame(C1 = c(0, 0.5, 1, 2), Res = c(10, 8, 5, 2))
  d <- dr_curve_data(NULL, df)
  expect_null(d$curve)
  expect_equal(d$control_x, 0.5 / 10)
  expect_equal(d$observed$conc, df$C1)
})

test_that("dr_curve_data samples the curve densely in the low-dose bend (issue #6)", {
  # Sam's chlorfenapyr data spans ~9 decades (1e-6 .. 900) but the whole sigmoid
  # lives below conc ~5. A linear grid put only ~2 of 200 points there, so the
  # bend was drawn as a single straight segment on the log x-axis. A log-spaced
  # grid keeps the low-dose region well sampled regardless of range.
  csv <- testthat::test_path("fixtures", "single", "single_chlorfenapyr_continuous.csv")
  df  <- to_engine_df(read_upload(csv), "single")
  fit <- analyse_single(df)
  d   <- dr_curve_data(fit, df)

  expect_gt(sum(d$curve$conc < 5), 100)            # was 2 with a linear grid
  # grid is geometric: successive ratios are (near) constant
  ratios <- d$curve$conc[-1] / d$curve$conc[-length(d$curve$conc)]
  expect_lt(diff(range(ratios)), 1e-6)
  # bounds unchanged: still spans lowest dose to highest dose
  expect_equal(min(d$curve$conc), min(df$C1))
  expect_equal(max(d$curve$conc), max(df$C1))
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

test_that("isobole_data returns contour paths at the requested effect levels", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 0.5, 1), C2 = c(0, 0.5, 1),
                   Affected = c(10, 5, 1), Exposed = rep(10, 3))
  d <- isobole_data(fit, df, levels = c(0.25, 0.5, 0.75), n = 40)
  expect_true(all(c("source", "level", "group", "x", "y") %in% names(d)))
  expect_setequal(unique(d$source), "fit")
  expect_setequal(sort(unique(d$level)), sort(c(0.25, 0.5, 0.75) * 1))
})

test_that("isobole_data adds a reference set when reference_fit is given", {
  fit <- list(reference = "CA", deviation = "SA", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5, a = 2))
  ref <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 0.5, 1), C2 = c(0, 0.5, 1),
                   Affected = c(10, 5, 1), Exposed = rep(10, 3))
  d <- isobole_data(fit, df, levels = 0.5, reference_fit = ref, n = 40)
  expect_setequal(unique(d$source), c("fit", "reference"))
})

test_that("isobole_data rejects out-of-range levels", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 1), C2 = c(0, 1), Affected = c(10, 1), Exposed = c(10, 10))
  expect_error(isobole_data(fit, df, levels = c(0.5, 1.5)), "between 0 and 1")
})
