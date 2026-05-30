test_that("engine reproduces the binary workbook CA continuous analysis", {
  csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
  skip_if_not(file.exists(csv), "fixture CSV not generated")
  df <- read.csv(csv)

  set.seed(42)  # multi-start uses runif(); fix the seed for a deterministic test
  res <- analyse_mixture(df, reference = "CA", response = "continuous")

  # Reference fit: residual SS within 2% of the workbook (1,633,769.68).
  expect_equal(res$fits$reference$objective, 1633769.68, tolerance = 0.02)
  # EC50 of the toxic chemical (CPF/C1) near the workbook value 0.082.
  expect_equal(unname(res$fits$reference$par["ec501"]), 0.082, tolerance = 0.05)
  # Workbook concludes DL is the most parsimonious (CA vs DL p ~ 0).
  expect_equal(res$chosen, "DL")
  # CA-vs-S/A LR statistic in the right ballpark (workbook 13.17).
  ca_sa <- res$comparison[res$comparison$model == "SA", ]
  expect_gt(ca_sa$chi, 8)
  expect_lt(ca_sa$p, 0.01)
})

test_that("engine reproduces the binary workbook CA quantal analysis", {
  csv <- testthat::test_path("fixtures", "binary_mps_cpf_quantal.csv")
  skip_if_not(file.exists(csv), "quantal fixture CSV not generated")
  df <- read.csv(csv)   # columns C1, C2, Affected, Exposed

  set.seed(42)
  res <- analyse_mixture(df, reference = "CA", response = "binary")

  # Reference residual deviance within 2% of the workbook (184.424).
  expect_equal(res$fits$reference$objective, 184.424, tolerance = 0.02)
  # control proportion (max) near the workbook value 0.9556.
  expect_equal(unname(res$fits$reference$par["max"]), 0.9556, tolerance = 0.02)
  # Workbook: no deviation is significant (all p > 0.3) -> reference chosen.
  expect_equal(res$chosen, "reference")
})
