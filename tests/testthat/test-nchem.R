test_that("model_spec handles ternary CA selections", {
  ref <- model_spec("CA", "reference", n_chem = 3)
  expect_equal(ref$params,
               c("max", "slope1", "slope2", "slope3",
                 "ec50_1", "ec50_2", "ec50_3"))
  expect_identical(ref$fn, ca_tri_vec)
  expect_equal(model_spec("CA", "SA", 3)$extra, "a")
})

make_staged_ternary_df <- function() {
  # truth: max=800; slopes 4,1.5,1; ec50s 0.08,1,5 (matches ca_tri_vec arg order)
  c1 <- c(0.01, 0.02, 0.04, 0.08, 0.16)
  c2 <- c(0.125, 0.25, 0.5, 1, 2)
  c3 <- c(0.5, 1, 2, 4, 8)
  singles <- rbind(
    data.frame(C1 = c(0, c1), C2 = 0,        C3 = 0),
    data.frame(C1 = 0,        C2 = c(0, c2), C3 = 0),
    data.frame(C1 = 0,        C2 = 0,        C3 = c(0, c3)))
  mix <- expand.grid(C1 = c(0.04, 0.08), C2 = c(0.5, 1), C3 = c(1, 2))
  d <- unique(rbind(singles, mix))
  d$Res <- ca_tri_vec(d$C1, d$C2, d$C3, 800, 4, 1.5, 1, 0.08, 1, 5)
  d
}

test_that("analyse_mixture runs end-to-end (staged) on a ternary dataset", {
  g <- make_staged_ternary_df()
  res <- analyse_mixture(g, reference = "CA", response = "continuous")
  # DR / DL are binary-only: ternary fits only reference + S/A here.
  expect_setequal(names(res$fits), c("reference", "SA"))
  # Curve params are fixed (from the single-compound rows) across both fits:
  base <- c("max", "slope1", "slope2", "slope3", "ec50_1", "ec50_2", "ec50_3")
  expect_equal(res$fits$SA$par[base], res$fits$reference$par[base])
  # df counts only the free interaction params (S/A carries a single `a`).
  expect_equal(unname(res$fits$reference$df), 0)
  expect_equal(unname(res$fits$SA$df), 1)
  # No real interaction in this reference-only data: `a` stays ~0. We assert the
  # parameter rather than model SELECTION because on noiseless data the
  # reference/SA objectives are ~1e-15 and the LR test on them is numerically
  # degenerate (it would otherwise pick a spurious deviation with a ~ -1e-8).
  expect_equal(unname(res$fits$SA$par[["a"]]), 0, tolerance = 1e-2)
})

test_that("ternary DR/DL deviations are rejected (binary-only)", {
  expect_error(model_spec("CA", "DR", 3), "binary-only")
  expect_error(model_spec("IA", "DL", 3), "binary-only")
  # ...but they remain available for binary mixtures.
  expect_equal(model_spec("CA", "DR", 2)$extra, c("a", "b"))
})

test_that("single-chemical analysis returns just the dose-response fit", {
  d <- data.frame(C1 = c(0, 0.01, 0.03, 0.1, 0.3, 1, 3),
                  Res = ll3_predict(c(0, 0.01, 0.03, 0.1, 0.3, 1, 3), 800, 2, 0.5))
  fit <- analyse_single(d)
  expect_equal(unname(fit$par["ec50"]), 0.5, tolerance = 1e-2)
})
