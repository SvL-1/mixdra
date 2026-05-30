make_binary_df <- function() {
  # Tiny synthetic 2-chem dataset, decreasing endpoint
  expand <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  truth <- list(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expand$Res <- ca_bi_vec(expand$C1, expand$C2, truth$max, truth$slope1,
                          truth$slope2, truth$ec501, truth$ec502)
  expand
}

test_that("fit_model (continuous CA) recovers parameters and reports SS objective", {
  df <- make_binary_df()
  start <- c(max = 700, slope1 = 2, slope2 = 1, ec501 = 0.1, ec502 = 2)
  fit <- fit_model(df, reference = "CA", deviation = "reference",
                   response = "continuous", start = start)
  expect_equal(unname(fit$par["ec501"]), 0.08, tolerance = 1e-2)
  expect_lt(fit$objective, 1e-2)          # near-perfect fit on noiseless data
  expect_equal(fit$df, 5)                  # 5 free parameters
})

test_that("fit_model respects fixed parameters", {
  df <- make_binary_df()
  start <- c(max = 800, slope1 = 2, slope2 = 1, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   fixed = c("max", "ec501", "ec502"))
  expect_equal(unname(fit$par["max"]), 800)    # untouched
  expect_equal(unname(fit$par["ec501"]), 0.08) # untouched
  expect_equal(fit$df, 2)                       # only slope1, slope2 free
})

test_that("fit_model (binary) minimises deviance and predicts probabilities", {
  df <- make_binary_df()
  df$Exposed <- 10
  # Build affected counts from a known probability surface (max as proportion)
  p <- ca_bi_vec(df$C1, df$C2, 0.95, 4, 1.5, 0.08, 1)
  df$Affected <- round(p * df$Exposed)
  start <- c(max = 0.9, slope1 = 2, slope2 = 1, ec501 = 0.1, ec502 = 2)
  fit <- fit_model(df, "CA", "reference", response = "binary", start = start)
  expect_true(is.finite(fit$objective))
  expect_true(all(fit$pred >= 0 & fit$pred <= 1))
})

test_that("fit_model with all parameters fixed evaluates without optimising", {
  df <- make_binary_df()  # truth: max=800, slope1=4, slope2=1.5, ec501=0.08, ec502=1
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   fixed = names(start))
  expect_equal(fit$df, 0)
  expect_equal(fit$convergence, 0)
  expect_equal(unname(fit$par[["max"]]), 800)     # untouched
  expect_lt(fit$objective, 1e-6)                  # truth params -> ~perfect fit
})
