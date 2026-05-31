test_that("analyse_mixture returns a fit per deviation and a chosen model", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  res <- analyse_mixture(df, reference = "CA", response = "continuous")
  expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
  expect_true(res$chosen %in% c("reference", "SA", "DR", "DL"))
  # On noiseless reference-generated data the staged fit must not invent an
  # interaction: `a` stays ~0. We assert the parameter rather than model
  # SELECTION because the reference/SA objectives here are ~1e-15, which makes
  # the likelihood-ratio test degenerate (it would otherwise pick a spurious SA).
  expect_equal(unname(res$fits$SA$par[["a"]]), 0, tolerance = 1e-2)
  # Comparison table has one row per non-reference model.
  expect_equal(nrow(res$comparison), 3)
  expect_true(all(c("model", "chi", "df", "p") %in% names(res$comparison)))
})
