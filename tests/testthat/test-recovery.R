# Recovery oracle: data generated forward from known params must be inverted back
# to those params by the engine. At zero noise the true model's objective -> 0, so
# the LR-test selection picks the generating deviation. These tests prove
# correctness (params in == params out), not just regression-stability.
#
# Coverage is a representative subset of the model grid BY DESIGN (the full
# reference x deviation x response cross-product would be 16+ slow multi-start
# fits): CA reference (continuous), IA reference (binary), CA SA/DR (continuous),
# CA DL (binary), and CA ternary reference. IA deviations, binary SA/DR, and
# ternary deviations are intentionally not covered here; ternary ASA is on hold.

# Mute ONLY the benign base-R "one-dimensional optimization is unreliable"
# warning that stats::optim (Nelder-Mead) emits for single-free-parameter
# sub-fits. Any other warning still surfaces, so a genuine new regression in the
# fitter is not hidden.
mute_1d <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("one-dimensional optimization", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}

test_that("recovery: continuous CA reference round-trips", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  df <- simulate_mixture(par, "CA", "reference", "continuous")
  res <- mute_1d(analyse_mixture(df, reference = "CA", response = "continuous", n_starts = 1))
  ref <- res$fits$reference$par
  expect_equal(unname(ref["max"]),   800,  tolerance = 1e-2)
  expect_equal(unname(ref["slope1"]), 4,   tolerance = 1e-2)
  expect_equal(unname(ref["ec501"]), 0.08, tolerance = 1e-2)
  expect_equal(unname(ref["ec502"]), 1,    tolerance = 1e-2)
  expect_lt(res$fits$reference$objective, 1e-3)
  expect_equal(res$chosen, "reference")
})

test_that("recovery: binary IA reference round-trips", {
  par <- c(max = 0.95, slope1 = 3, slope2 = 2, ec501 = 0.1, ec502 = 0.5)
  df <- simulate_mixture(par, "IA", "reference", "binary")
  res <- mute_1d(analyse_mixture(df, reference = "IA", response = "binary", n_starts = 1))
  ref <- res$fits$reference$par
  expect_equal(unname(ref["max"]),   0.95, tolerance = 1e-2)
  expect_equal(unname(ref["ec501"]), 0.1,  tolerance = 1e-2)
  expect_lt(res$fits$reference$objective, 1e-4)
  expect_equal(res$chosen, "reference")
})

# The deviation tests below fit interaction params (a, b) via a seeded 10-start
# optimiser. Tolerances (0.1-0.15) absorb cross-platform multi-start variance
# (runif streams / BLAS / optimiser convergence can differ across R builds). If
# one of these fails on a new platform, INVESTIGATE THE FIT -- do not blindly
# loosen the tolerance (that guts the oracle), and do not tighten below the
# multi-start noise floor.
test_that("recovery: continuous CA DR recovers a, b and is selected", {
  skip_on_cran()
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1,
           a = 2.5, b = 1.5)
  df <- simulate_mixture(par, "CA", "DR", "continuous")
  set.seed(1)
  res <- mute_1d(analyse_mixture(df, reference = "CA", response = "continuous",
                         n_starts = 10))
  dr <- res$fits$DR$par
  expect_equal(unname(dr["ec501"]), 0.08, tolerance = 1e-2)  # curve from singles
  expect_equal(unname(dr["a"]), 2.5, tolerance = 0.1)
  expect_equal(unname(dr["b"]), 1.5, tolerance = 0.1)
  expect_lt(res$fits$DR$objective, 1e-2)
  expect_equal(res$chosen, "DR")
})

test_that("recovery: continuous CA SA recovers a and is selected", {
  skip_on_cran()
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 3)
  df <- simulate_mixture(par, "CA", "SA", "continuous")
  set.seed(2)
  res <- mute_1d(analyse_mixture(df, reference = "CA", response = "continuous",
                         n_starts = 10))
  expect_equal(unname(res$fits$SA$par["a"]), 3, tolerance = 0.1)
  expect_lt(res$fits$SA$objective, 1e-2)
  expect_equal(res$chosen, "SA")
})

test_that("recovery: binary CA DL recovers a, b and is selected", {
  skip_on_cran()
  par <- c(max = 0.9, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1,
           a = 2, b = 0.5)
  df <- simulate_mixture(par, "CA", "DL", "binary")
  set.seed(3)
  res <- mute_1d(analyse_mixture(df, reference = "CA", response = "binary",
                         n_starts = 10))
  dl <- res$fits$DL$par
  expect_equal(unname(dl["a"]), 2,   tolerance = 0.15)
  expect_equal(unname(dl["b"]), 0.5, tolerance = 0.15)
  expect_equal(res$chosen, "DL")
})

# Ternary: reference-level recovery only. Advanced S/A (ASA) fitting is on hold,
# so its round-trip is intentionally not validated here.
test_that("recovery: continuous CA ternary reference round-trips", {
  par <- c(max = 100, slope1 = 2, slope2 = 1.5, slope3 = 3,
           ec50_1 = 0.1, ec50_2 = 0.5, ec50_3 = 2)
  df <- simulate_mixture(par, "CA", "reference", "continuous")
  res <- mute_1d(analyse_mixture(df, reference = "CA", response = "continuous", n_starts = 1))
  ref <- res$fits$reference$par
  expect_equal(unname(ref["ec50_1"]), 0.1, tolerance = 2e-2)
  expect_equal(unname(ref["ec50_3"]), 2,   tolerance = 2e-2)
  expect_equal(unname(ref["max"]),    100, tolerance = 2e-2)
  expect_lt(res$fits$reference$objective, 1e-2)
})
