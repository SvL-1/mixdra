test_that("fit_model with n_starts returns the best of several starts", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  # A deliberately poor single start; multi-start should still find a good fit.
  bad <- c(max = 400, slope1 = 0.5, slope2 = 0.5, ec501 = 1, ec502 = 0.1)
  set.seed(1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = bad, n_starts = 8)
  expect_lt(fit$objective, 1)
})
