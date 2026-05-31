# Tests for explicit per-parameter bounds in fit_model().
# Spec: docs/superpowers/specs/2026-05-30-parameter-constraints-design.md

bounds_cont_df <- function() {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
  g
}

bounds_binary_df <- function() {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Exposed <- 10
  p <- ca_bi_vec(g$C1, g$C2, 0.95, 4, 1.5, 0.08, 1)
  g$Affected <- round(p * g$Exposed)
  g
}

test_that("an explicit upper bound caps a base parameter", {
  df <- bounds_cont_df()
  # True optimum for ec501 is 0.08; cap it below that and it must sit at the cap.
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.04, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   upper = c(ec501 = 0.05))
  expect_lte(unname(fit$par["ec501"]), 0.05 + 1e-6)
})

test_that("an explicit lower bound floors a base parameter", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.25, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   lower = c(ec501 = 0.2))
  expect_gte(unname(fit$par["ec501"]), 0.2 - 1e-6)
})

test_that("an unknown bound name raises an error", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expect_error(
    fit_model(df, "CA", "reference", "continuous", start = start,
              upper = c(nonsense = 1)),
    "base parameter")
})

test_that("naming a deviation parameter in bounds raises an error", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expect_error(
    fit_model(df, "CA", "SA", "continuous", start = start,
              lower = c(a = 0)),
    "base parameter")
})

test_that("binary max defaults to an upper bound of 1", {
  df <- bounds_binary_df()
  start <- c(max = 0.9, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "binary", start = start)
  expect_lte(unname(fit$par["max"]), 1 + 1e-8)
})

test_that("binary max upper bound can be tightened below 1", {
  df <- bounds_binary_df()
  start <- c(max = 0.9, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "binary", start = start,
                   upper = c(max = 0.98))
  expect_lte(unname(fit$par["max"]), 0.98 + 1e-6)
})

test_that("a binary max upper above 1 is clamped to 1 with a warning", {
  df <- bounds_binary_df()
  start <- c(max = 0.9, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expect_warning(
    fit <- fit_model(df, "CA", "reference", "binary", start = start,
                     upper = c(max = 1.5)),
    "capped")
  expect_lte(unname(fit$par["max"]), 1 + 1e-8)
})

test_that("a start outside its bounds is clamped with a warning", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.1, ec502 = 1)
  expect_warning(
    fit_model(df, "CA", "reference", "continuous", start = start,
              upper = c(ec501 = 0.05)),
    "clamped")
})

test_that("deviation parameters stay unconstrained (free to go negative)", {
  # A DL fit on data generated with a negative interaction should let `a`
  # move below 0 -- proof the deviation params are not floored at positivity.
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "DL", "continuous", start = start, n_starts = 1)
  expect_true("a" %in% names(fit$par))
  # `a` is unbounded; the fit must not error and must keep it finite.
  expect_true(is.finite(unname(fit$par["a"])))
})

test_that("an illegal bound name errors even when all params are fixed", {
  # The all-fixed fast path must still validate bound names; otherwise the staged
  # reference fit would silently ignore an illegal bound that the other fits reject.
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expect_error(
    fit_model(df, "CA", "reference", "continuous", start = start,
              fixed = names(start), upper = c(a = 1)),
    "base parameter")
})

test_that("analyse_mixture forwards bounds to every fit", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.04, ec502 = 1)
  res <- analyse_mixture(df, "CA", "continuous", start = start,
                         upper = c(ec501 = 0.05), n_starts = 1)
  expect_lte(unname(res$fits$reference$par["ec501"]), 0.05 + 1e-6)
  expect_lte(unname(res$fits$DL$par["ec501"]), 0.05 + 1e-6)
})
