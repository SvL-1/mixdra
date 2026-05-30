test_that("model_spec handles ternary CA selections", {
  ref <- model_spec("CA", "reference", n_chem = 3)
  expect_equal(ref$params,
               c("max", "slope1", "slope2", "slope3",
                 "ec50_1", "ec50_2", "ec50_3"))
  expect_identical(ref$fn, ca_tri_vec)
  expect_equal(model_spec("CA", "SA", 3)$extra, "a")
})

test_that("analyse_mixture runs end-to-end on a ternary dataset", {
  g <- expand.grid(C1 = c(0, 0.05), C2 = c(0, 0.5), C3 = c(0, 2))
  g$Res <- ca_tri_vec(g$C1, g$C2, g$C3, 800, 4, 1.5, 1, 0.08, 1, 5)
  res <- analyse_mixture(g, reference = "CA", response = "continuous")
  expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
  expect_equal(res$chosen, "reference")  # noiseless reference data
})

test_that("single-chemical analysis returns just the dose-response fit", {
  d <- data.frame(C1 = c(0, 0.01, 0.03, 0.1, 0.3, 1, 3),
                  Res = ll3_predict(c(0, 0.01, 0.03, 0.1, 0.3, 1, 3), 800, 2, 0.5))
  fit <- analyse_single(d)
  expect_equal(unname(fit$par["ec50"]), 0.5, tolerance = 1e-2)
})
