test_that("classify_rows labels control/single/binary/ternary", {
  df <- data.frame(C1 = c(0, 1, 0, 1, 1),
                   C2 = c(0, 0, 2, 2, 2),
                   C3 = c(0, 0, 0, 0, 3))
  expect_equal(as.character(classify_rows(df)),
               c("control", "single", "single", "binary", "ternary"))
})

test_that("ternary_ratio_key groups fixed proportions across dose levels", {
  # two ratios: (0.2,0.6,0.2) at 3 doses, (0.6,0.2,0.2) at 2 doses
  df <- data.frame(
    C1 = c(0.2, 0.4, 1.0, 0.6, 1.2),
    C2 = c(0.6, 1.2, 3.0, 0.2, 0.4),
    C3 = c(0.2, 0.4, 1.0, 0.2, 0.4))
  keys <- ternary_ratio_key(df)
  expect_equal(length(unique(keys)), 2)
  expect_equal(keys[1], keys[2]); expect_equal(keys[2], keys[3])
  expect_false(keys[1] == keys[4])
  expect_equal(keys[4], keys[5])
})

test_that("fit_model holds a fixed deviation parameter without error", {
  df <- data.frame(C1 = c(0, 1, 2, 0, 0, 1),
                   C2 = c(0, 0, 0, 1, 2, 1),
                   Res = c(800, 500, 300, 450, 250, 200))
  # `a` is a deviation param; fixing it must not corrupt parscale / perturbation.
  expect_error(
    fit_model(df, "CA", "SA", "continuous",
              start = c(max = 800, slope1 = 1, slope2 = 1,
                        ec501 = 1, ec502 = 1, a = 0),
              fixed = "a", n_starts = 2),
    NA)  # NA => assert NO error is thrown
})
