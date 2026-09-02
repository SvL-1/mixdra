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

test_that("fit_single honours a partial upper bound", {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  resp <- ll3_predict(conc, 100, 2, 0.5)   # true slope 2
  fit <- fit_single(conc, resp, upper = c(slope = 1.5))
  expect_lte(fit$par[["slope"]], 1.5 + 1e-6)
})

test_that("fit_single clamps an out-of-bounds start instead of erroring", {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  resp <- ll3_predict(conc, 100, 2, 0.5)
  expect_error(
    fit_single(conc, resp, lower = c(ec50 = 1), start = c(ec50 = 0.001)),
    NA)
})

test_that("eval_single returns the fit shape without optimising", {
  conc <- c(0, 0.5, 1, 2, 4)
  resp <- ll3_predict(conc, 10, 2, 1)
  ev <- eval_single(conc, resp, 10, 2, 1)
  expect_equal(ev$kind, "single")
  expect_equal(unname(ev$par[c("max", "slope", "ec50")]), c(10, 2, 1))
  expect_equal(ev$ssr, 0)
  expect_true(is.na(ev$convergence))
})

test_that("analyse_single forwards bounds to fit_single", {
  d <- data.frame(C1 = c(0, 0.1, 0.3, 1, 3, 10),
                  Res = ll3_predict(c(0, 0.1, 0.3, 1, 3, 10), 100, 2, 0.5))
  fit <- analyse_single(d, upper = c(slope = 1.5))
  expect_lte(fit$par[["slope"]], 1.5 + 1e-6)
})

test_that("fit_single holds a pinned parameter at its start value", {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  d <- data.frame(C1 = conc, Res = ll3_predict(conc, 100, 2, 0.5))
  fit <- analyse_single(d, start = c(ec50 = 0.7), fixed = "ec50")
  expect_equal(unname(fit$par[["ec50"]]), 0.7)
  expect_true(fit$ssr > 0)          # the pin costs fit quality
})

test_that("fit_single with every parameter pinned just evaluates the curve", {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  d <- data.frame(C1 = conc, Res = ll3_predict(conc, 100, 2, 0.5))
  fit <- analyse_single(d, start = c(max = 100, slope = 2, ec50 = 0.5),
                        fixed = c("max", "slope", "ec50"))
  expect_equal(unname(fit$par[c("max", "slope", "ec50")]), c(100, 2, 0.5))
  expect_equal(fit$ssr, 0)
  expect_equal(fit$convergence, 0L)
})
