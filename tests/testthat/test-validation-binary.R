# Validation of the STAGED fitting method on the binary MPs+CPF workbook data.
#
# IMPORTANT: these numbers intentionally differ from the published MixTox/Excel
# workbook. The workbook fits the curve parameters (max, slopes, EC50s) using the
# mixture data (joint Jonker-2005 fitting); this engine fixes the curve
# parameters from the single-compound data ONLY, then fits the interaction
# parameters (a, b) to the mixture data. Constraining the reference curves to the
# marginals raises the reference residual (the workbook instead absorbs
# interaction into re-fitted curves) and lets genuine deviations surface -- e.g.
# the quantal data below now selects a dose-ratio (DR) deviation that the joint
# workbook reported as 'reference'.

test_that("staged engine: binary CA continuous analysis on MPs+CPF", {
  csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi", "binary_ca_mps_cpf_imi_continuous.csv")
  skip_if_not(file.exists(csv), "fixture CSV not generated")
  df <- read.csv(csv)

  set.seed(42)
  res <- analyse_mixture(df, reference = "CA", response = "continuous",
                         n_starts = 1)

  # Staged reference residual SS (workbook joint fit was 1,633,769.68).
  expect_equal(res$fits$reference$objective, 2224215.9, tolerance = 0.02)
  # EC50 of the toxic chemical (CPF/C1) from its single-compound curve.
  expect_equal(unname(res$fits$reference$par["ec501"]), 0.094, tolerance = 0.05)
  # Still selects the dose-level (DL) deviation as most parsimonious.
  expect_equal(res$chosen, "DL")
  # S/A remains significant, though weaker than the joint workbook's chi ~13.
  ca_sa <- res$comparison[res$comparison$model == "SA", ]
  expect_gt(ca_sa$chi, 3)
  expect_lt(ca_sa$p, 0.05)
})

test_that("staged engine: binary CA quantal analysis on MPs+CPF", {
  csv <- testthat::test_path("fixtures", "binary", "survival", "binary_mps_cpf_quantal.csv")
  skip_if_not(file.exists(csv), "quantal fixture CSV not generated")
  df <- read.csv(csv)   # columns C1, C2, Affected, Exposed

  set.seed(42)
  res <- analyse_mixture(df, reference = "CA", response = "binary",
                         n_starts = 1)

  # Staged reference residual deviance (workbook joint fit was 184.424).
  expect_equal(res$fits$reference$objective, 231.418, tolerance = 0.02)
  # control proportion (max) near the workbook value 0.9556.
  expect_equal(unname(res$fits$reference$par["max"]), 0.9642, tolerance = 0.02)
  # Staged fitting detects a dose-ratio-dependent (DR) deviation that the joint
  # workbook (which re-fits the curves on the mixture) reports as 'reference'.
  expect_equal(res$chosen, "DR")
  # The interaction is strongly supported (CA-vs-S/A p well below 0.001).
  ca_sa <- res$comparison[res$comparison$model == "SA", ]
  expect_lt(ca_sa$p, 0.001)
})
