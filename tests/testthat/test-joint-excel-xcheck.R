# Acceptance/smoke check: joint refinement on the real binary CPF/IMI(MPs) data
# must run, never worsen the staged objective, and land at finite, plausible
# params. NOT a bit-exact Excel match (different optimizer; CA bisection is
# non-smooth). Spec: docs/superpowers/specs/2026-06-02-binary-joint-fit-design.md

test_that("joint fit on real binary data runs and never worsens the staged fit", {
  skip_on_cran()
  csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                             "binary_ca_mps_cpf_imi_continuous.csv")
  skip_if_not(file.exists(csv), "binary fixture missing")
  df <- utils::read.csv(csv)

  res <- analyse_mixture(df, reference = "CA", response = "continuous", n_starts = 1)
  staged <- res$fits[["SA"]]
  joint  <- refine_joint(staged, df, n_starts = 4, time_limit = 30)

  expect_true(isTRUE(joint$joint))
  expect_lte(joint$objective, staged$objective + 1e-6)   # never worse than the seed
  expect_true(all(is.finite(joint$par)))                 # all params finite
  expect_gt(unname(joint$par[["max"]]), 0)               # plausible control plateau
  # Surface the fitted values for human comparison against the workbook.
  message(sprintf("joint fit: max=%.3f a=%.4f SSR=%.4f (staged SSR=%.4f)",
                  joint$par[["max"]], joint$par[["a"]],
                  joint$objective, staged$objective))
})
