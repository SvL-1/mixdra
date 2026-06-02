# Joint model-selection loop: every model is fully optimised (curves free),
# warm-started down the nesting chain. See
# docs/superpowers/specs/2026-06-02-binary-model-selection-loop-ux-design.md

test_that("joint_chain_order is reference->SA->DR->DL for binary, reference->SA for ternary", {
  expect_equal(joint_chain_order(2), c("reference", "SA", "DR", "DL"))
  expect_equal(joint_chain_order(3), c("reference", "SA"))
})

test_that("joint_fit_one frees the curves and (warm-started) never worsens the parent", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2, 0.8), C2 = c(0, 0.5, 5, 20))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  seed <- c(max = 700, slope1 = 3, slope2 = 2, ec501 = 0.1, ec502 = 1.5)  # off-true seed

  ref <- joint_fit_one(df, "CA", "reference", "continuous", seed)
  expect_equal(ref$df, 5L)                       # all five curve params were FREE
  expect_false(isTRUE(ref$simulated))

  sa <- joint_fit_one(df, "CA", "SA", "continuous", seed, parent_fit = ref)
  expect_equal(sa$df, 6L)                         # + a
  expect_true("a" %in% names(sa$par))
  expect_lte(sa$objective, ref$objective + 1e-6)  # warm-start monotonicity

  dr <- joint_fit_one(df, "CA", "DR", "continuous", seed, parent_fit = sa)
  expect_equal(dr$df, 7L)                         # + b
  expect_lte(dr$objective, sa$objective + 1e-6)
})
