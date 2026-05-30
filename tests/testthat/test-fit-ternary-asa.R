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
