# Synthetic binary data with a full single-compound dose series per chemical
# (>= 6 distinct marginal concentrations) so stage 1 can identify each curve,
# plus a few interior mixture points. Noiseless CA-reference surface.
make_staged_binary_df <- function(truth = list(max = 800, slope1 = 4,
                                               slope2 = 1.5, ec501 = 0.08,
                                               ec502 = 1)) {
  c1 <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
  c2 <- c(0.125, 0.25, 0.5, 1, 2, 4)
  singles <- rbind(data.frame(C1 = c(0, c1), C2 = 0),
                   data.frame(C1 = 0,        C2 = c2))
  mix <- expand.grid(C1 = c(0.04, 0.08), C2 = c(0.5, 1))
  d <- unique(rbind(singles, mix))
  d$Res <- ca_bi_vec(d$C1, d$C2, truth$max, truth$slope1, truth$slope2,
                     truth$ec501, truth$ec502)
  d
}

test_that("fit_curve_from_singles recovers curve params from the marginals", {
  d <- make_staged_binary_df()
  base <- fit_curve_from_singles(d, "CA", "continuous", n_starts = 1)
  expect_setequal(names(base), c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_equal(unname(base[["max"]]),   800,  tolerance = 1e-2)
  expect_equal(unname(base[["ec501"]]), 0.08, tolerance = 1e-2)
  expect_equal(unname(base[["ec502"]]), 1,    tolerance = 1e-2)
  expect_equal(unname(base[["slope1"]]), 4,   tolerance = 1e-2)
  expect_equal(unname(base[["slope2"]]), 1.5, tolerance = 1e-2)
})

test_that("analyse_mixture holds the curve params fixed across all four fits", {
  d <- make_staged_binary_df()
  res <- analyse_mixture(d, "CA", "continuous", n_starts = 1)
  base <- c("max", "slope1", "slope2", "ec501", "ec502")
  ref_base <- res$fits$reference$par[base]
  for (m in c("SA", "DR", "DL"))
    expect_equal(res$fits[[m]]$par[base], ref_base)
  # df now counts only the free interaction params.
  expect_equal(res$fits$reference$df, 0)
  expect_equal(res$fits$SA$df, 1)   # a
  expect_equal(res$fits$DR$df, 2)   # a, b
  expect_equal(res$fits$DL$df, 2)   # a, b
})

test_that("analyse_mixture leaves the interaction near 0 when there is none", {
  # On reference-only data, the staged fit must not invent an interaction. We
  # assert the fitted `a` directly rather than the model SELECTION, because on
  # perfectly noiseless data the reference/SA objectives are ~1e-15 and the
  # likelihood-ratio test on such near-zero objectives is numerically degenerate.
  d <- make_staged_binary_df()
  res <- analyse_mixture(d, "CA", "continuous", n_starts = 1)
  expect_equal(unname(res$fits$SA$par[["a"]]), 0, tolerance = 1e-2)
})

test_that("analyse_mixture recovers a known S/A interaction with fixed curves", {
  truth <- list(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d <- make_staged_binary_df(truth)
  # Regenerate the response from a CA + S/A surface with a = 5. On single-compound
  # rows the S/A term is 0, so the marginals are unchanged and stage 1 still
  # recovers the true curves; the interaction lives only in the mixture rows.
  d$Res <- ca_sa_bi_vec(d$C1, d$C2, truth$max, truth$slope1, truth$slope2,
                        truth$ec501, truth$ec502, a = 5)
  set.seed(1)
  d$Res <- d$Res + stats::rnorm(nrow(d), sd = 2)   # small noise -> well-posed LR test
  res <- analyse_mixture(d, "CA", "continuous", n_starts = 1)
  # Curve params come from the (interaction-free) marginals and stay fixed:
  base <- c("max", "slope1", "slope2", "ec501", "ec502")
  expect_equal(res$fits$SA$par[base], res$fits$reference$par[base])
  # The interaction is detected and `a` is recovered near its true value:
  expect_equal(unname(res$fits$SA$par[["a"]]), 5, tolerance = 0.2)
  expect_false(res$chosen == "reference")          # an interaction is selected
})
