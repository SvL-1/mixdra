test_that("fit_model records model descriptors on its result", {
  df <- read.csv(testthat::test_path("fixtures", "binary", "survival", "binary_mps_cpf_quantal.csv"))
  start <- seed_from_singles(df, "binary")
  set.seed(1)
  fit <- fit_model(df, "CA", "SA", "binary", start = start, n_starts = 1)
  expect_equal(fit$reference, "CA")
  expect_equal(fit$deviation, "SA")
  expect_equal(fit$response, "binary")
  expect_equal(fit$conc_cols, c("C1", "C2"))
  expect_equal(fit$n_chem, 2)
  expect_equal(fit$kind, "mixture")
})

test_that("fit_single tags its result kind", {
  f <- fit_single(c(0, 0.5, 1, 2, 4), c(10, 8, 5, 2, 1))
  expect_equal(f$kind, "single")
})
