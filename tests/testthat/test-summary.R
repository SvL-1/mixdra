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
