# Behaviour of the parameter-constraint scheme (see
# docs/parameter-constraints-decisions.md): default bounds are positivity-only
# for the base curve parameters, deviation params are unconstrained, binary
# `max` is capped at 1, and the user may override any bound manually.

make_cap_df <- function() {
  # CA surface whose chem-2 EC50 (8) sits far above a deliberately low seed, so
  # the retired seed-multiplier cap (upper = start * 10) could not reach it.
  g <- expand.grid(C1 = c(0, 0.05, 0.2, 0.5), C2 = c(0, 1, 4, 16))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 8)
  g
}

test_that("default base-parameter bounds are positivity-only (no seed cap)", {
  df <- make_cap_df()
  # Seed EC50_2 (0.5) is far below the truth (8); the old cap (0.5 * 10 = 5)
  # could never reach 8, positivity-only bounds can.
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 0.5)
  set.seed(1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start, n_starts = 6)
  expect_gt(unname(fit$par["ec502"]), 6)
})

test_that("a user-supplied upper bound is respected", {
  df <- make_cap_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 0.5)
  set.seed(1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   upper = c(ec502 = 3), n_starts = 6)
  expect_lte(unname(fit$par["ec502"]), 3 + 1e-6)
})

test_that("a user-supplied lower bound is respected", {
  df <- make_cap_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 0.5)
  set.seed(1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   lower = c(slope1 = 5), n_starts = 6)
  expect_gte(unname(fit$par["slope1"]), 5 - 1e-6)
})

test_that("binary control parameter keeps a default upper bound of 1", {
  df <- make_cap_df()
  df$Exposed <- 10
  df$Affected <- round(ca_bi_vec(df$C1, df$C2, 0.95, 4, 1.5, 0.08, 8) * 10)
  start <- c(max = 0.9, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 0.5)
  set.seed(1)
  fit <- fit_model(df, "CA", "reference", "binary", start = start, n_starts = 6)
  expect_lte(unname(fit$par["max"]), 1 + 1e-9)
  expect_true(all(fit$pred >= 0 & fit$pred <= 1))
})
