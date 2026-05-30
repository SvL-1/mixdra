test_that("ll3_predict matches the log-logistic formula", {
  # Y = max / (1 + (C/EC50)^slope); at C = EC50, Y = max/2
  expect_equal(ll3_predict(c(0, 1, Inf), max = 100, slope = 2, ec50 = 1),
               c(100, 50, 0))
})

test_that("fit_single recovers known parameters from clean data", {
  conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
  truth <- list(max = 800, slope = 2, ec50 = 0.5)
  resp <- ll3_predict(conc, truth$max, truth$slope, truth$ec50)
  fit <- fit_single(conc, resp)
  expect_equal(unname(fit$par["max"]),  800, tolerance = 1e-3)
  expect_equal(unname(fit$par["slope"]), 2,   tolerance = 1e-3)
  expect_equal(unname(fit$par["ec50"]),  0.5, tolerance = 1e-3)
  expect_lt(fit$ssr, 1e-6)
})
