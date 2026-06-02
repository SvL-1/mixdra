# eval_mixture(): forward-evaluate a mixture model at fixed curve + interaction
# params, returning the fit_model() shape. Uses the noiseless CA surfaces from
# the engine's own vectorised model functions as ground truth.

cp <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)

make_grid <- function() {
  d <- expand.grid(C1 = c(0, 0.04, 0.08, 0.16), C2 = c(0, 0.5, 1, 2))
  d
}

test_that("eval_mixture reproduces the CA reference surface (objective ~ 0)", {
  d <- make_grid()
  d$Res <- ca_bi_vec(d$C1, d$C2, cp[["max"]], cp[["slope1"]], cp[["slope2"]],
                     cp[["ec501"]], cp[["ec502"]])
  fit <- eval_mixture(d, "CA", "reference", "continuous", cp)
  expect_equal(fit$objective, 0, tolerance = 1e-8)
  expect_equal(unname(fit$pred), d$Res, tolerance = 1e-8)
  expect_true(fit$simulated)
  expect_equal(fit$kind, "mixture")
  expect_equal(fit$deviation, "reference")
  expect_setequal(names(fit$par), c("max", "slope1", "slope2", "ec501", "ec502"))
})

test_that("eval_mixture reproduces a known CA + S/A surface from entered a", {
  d <- make_grid()
  d$Res <- ca_sa_bi_vec(d$C1, d$C2, cp[["max"]], cp[["slope1"]], cp[["slope2"]],
                        cp[["ec501"]], cp[["ec502"]], a = 5)
  fit <- eval_mixture(d, "CA", "SA", "continuous", cp, interaction = c(a = 5))
  expect_equal(fit$objective, 0, tolerance = 1e-6)
  expect_equal(fit$par[["a"]], 5)
  expect_true("a" %in% names(fit$par))
  expect_equal(length(fit$pred), nrow(d))
})

test_that("eval_mixture drops interaction params the model does not use", {
  d <- make_grid()
  d$Res <- ca_sa_bi_vec(d$C1, d$C2, cp[["max"]], cp[["slope1"]], cp[["slope2"]],
                        cp[["ec501"]], cp[["ec502"]], a = 5)
  # pass a stray `b`; S/A has no b, so it must be ignored (not error)
  fit <- eval_mixture(d, "CA", "SA", "continuous", cp, interaction = c(a = 5, b = 99))
  expect_false("b" %in% names(fit$par))
  expect_equal(fit$par[["a"]], 5)
})
