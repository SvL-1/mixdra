test_that("plot_dose_response builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "mixture", conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 2, slope2 = 1, ec501 = 0.5, ec502 = 0.3))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), C2 = 0,
                   Affected = c(10, 6, 3, 1), Exposed = rep(10, 4))
  p <- plot_dose_response(fit, df, chem = 1)
  expect_s3_class(p, "plotly")
})

test_that("plot_dose_response builds a plotly object from raw data (no fit)", {
  skip_if_not_installed("plotly")
  df <- data.frame(C1 = c(0, 0.5, 1, 2), Res = c(10, 8, 5, 2))
  p <- plot_dose_response(NULL, df)
  expect_s3_class(p, "plotly")
})

test_that("plot_dose_response draws controls at the offset on a log axis", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), Res = c(10, 8, 5, 2))
  p <- plotly::plotly_build(plot_dose_response(fit, df, log_x = TRUE))
  ctrl <- Filter(function(tr) isTRUE(tr$name == "control"), p$x$data)
  expect_length(ctrl, 1)
  expect_equal(unique(ctrl[[1]]$x), 0.5 / 10)   # offset = min(positive)/10
})

test_that("plot_dose_response keeps controls at 0 on a linear axis", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), Res = c(10, 8, 5, 2))
  p <- plotly::plotly_build(plot_dose_response(fit, df, log_x = FALSE))
  expect_equal(p$x$layout$xaxis$type, "linear")
  # no separate "control" trace; the control sits at its true x = 0
  obs <- Filter(function(tr) isTRUE(tr$name == "observed"), p$x$data)
  expect_true(0 %in% obs[[1]]$x)
})

test_that("plot_obs_pred builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "mixture", pred = c(0.9, 0.5, 0.2))
  df <- data.frame(C1 = c(0, 1, 2), C2 = 0, Affected = c(9, 5, 2), Exposed = rep(10, 3))
  p <- plot_obs_pred(fit, df)
  expect_s3_class(p, "plotly")
})

test_that("plot_surface builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 0.9, slope1 = 1.5, slope2 = 1.2, ec501 = 0.1, ec502 = 0.2))
  df <- data.frame(C1 = c(0, 0.1, 0.3), C2 = c(0, 0.2, 0.4),
                   Affected = c(9, 5, 2), Exposed = rep(10, 3))
  p <- plot_surface(fit, df, n = 10)
  expect_s3_class(p, "plotly")
})

test_that("plot_isobole builds a plotly object, including a reference overlay", {
  skip_if_not_installed("plotly")
  fit <- list(reference = "CA", deviation = "SA", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5, a = 2))
  ref <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 0.5, 1), C2 = c(0, 0.5, 1),
                   Affected = c(10, 5, 1), Exposed = rep(10, 3))
  p <- plot_isobole(fit, df, levels = c(0.25, 0.5, 0.75), reference_fit = ref, n = 30)
  expect_s3_class(p, "plotly")
})
