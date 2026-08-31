# Tests for the final joint "Optimize all params" engine wrapper.
# Spec: docs/design/specs/2026-06-02-binary-joint-fit-design.md

test_that("joint refine recovers known params from noise-free data", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.5)
  df  <- simulate_mixture(par, "CA", "SA", "continuous", cv = 0)
  # A staged-style starting fit, perturbed away from truth.
  start_fit <- fit_model(df, "CA", "SA", "continuous",
                         start = par * c(1.2, 0.8, 1.3, 0.7, 1.1, 0.6),
                         fixed = c("max", "slope1", "slope2", "ec501", "ec502"),
                         n_starts = 1)
  joint <- refine_joint(start_fit, df, n_starts = 5, time_limit = 60)
  expect_equal(unname(joint$par[["a"]]),   1.5, tolerance = 0.05)
  expect_equal(unname(joint$par[["max"]]), 800, tolerance = 2)
  expect_true(isTRUE(joint$joint))
})

test_that("joint refine never worsens the seed objective", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.2)
  df  <- simulate_mixture(par, "CA", "SA", "continuous", cv = 0)
  staged <- fit_model(df, "CA", "SA", "continuous", start = par,
                      fixed = c("max", "slope1", "slope2", "ec501", "ec502"),
                      n_starts = 1)
  joint <- refine_joint(staged, df, n_starts = 1, time_limit = 60)
  expect_lte(joint$objective, staged$objective + 1e-6)
})

test_that("refine_joint pins a parameter when lower == upper", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.2)
  df  <- simulate_mixture(par, "CA", "SA", "continuous", cv = 0)
  staged <- fit_model(df, "CA", "SA", "continuous", start = par,
                      fixed = c("max", "slope1", "slope2", "ec501", "ec502"),
                      n_starts = 1)
  joint <- refine_joint(staged, df, lower = c(max = 800), upper = c(max = 800),
                        n_starts = 1, time_limit = 60)
  expect_equal(unname(joint$par[["max"]]), 800, tolerance = 1e-8)  # held
  expect_true("max" %in% joint$fixed)
})
