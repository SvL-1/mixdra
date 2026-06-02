test_that("param_ci returns finite lower/upper bounds for free parameters", {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1) + rep(c(-2, 2), length.out = 9)
  fit <- fit_model(g, "CA", "reference", "continuous",
                   start = c(max = 800, slope1 = 4, slope2 = 1.5,
                             ec501 = 0.08, ec502 = 1))
  ci <- param_ci(fit, g, "CA", "reference", "continuous")
  expect_true(all(c("parameter", "estimate", "lower", "upper") %in% names(ci)))
  expect_true(all(ci$lower <= ci$estimate & ci$estimate <= ci$upper))
})

test_that("result_table assembles one column per fitted model", {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
  res <- analyse_mixture(g, "CA", "continuous")
  tab <- result_table(res)
  expect_true(all(c("reference", "SA", "DR", "DL") %in% colnames(tab)))
  expect_true("max" %in% rownames(tab))
  expect_true("objective" %in% rownames(tab))
})

test_that("result_table folds in each model's LR test (chi-sq + p vs its parent)", {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
  res <- analyse_mixture(g, "CA", "continuous")
  tab <- result_table(res)
  expect_true(all(c("chi-sq (vs parent)", "p (vs parent)") %in% rownames(tab)))
  # the reference has no parent, so no chi-sq / p-value
  expect_true(is.na(tab["chi-sq (vs parent)", "reference"]))
  expect_true(is.na(tab["p (vs parent)", "reference"]))
  # SA/DR/DL each carry the chi-sq and p-value from the comparison table
  for (m in c("SA", "DR", "DL")) {
    expect_equal(tab["chi-sq (vs parent)", m],
                 res$comparison$chi[res$comparison$model == m])
    expect_equal(tab["p (vs parent)", m],
                 res$comparison$p[res$comparison$model == m])
  }
})

test_that("result_table handles a single fitted model with no comparison", {
  # mirrors the Binary tab populating the bottom table from one Autofit, before
  # Find best has run (so there is no $comparison and no chosen model).
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
  res <- analyse_mixture(g, "CA", "continuous")
  one <- result_table(list(fits = res$fits["SA"], comparison = NULL))
  expect_equal(colnames(one), "SA")
  expect_true(is.na(one["p (vs parent)", "SA"]))  # no comparison -> no p-value
  expect_false(is.na(one["objective", "SA"]))
})
