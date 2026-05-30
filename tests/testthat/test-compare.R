test_that("lr_test (continuous) matches the workbook CA-vs-S/A statistic", {
  # From MixTox Model_binary_MPs_CPF.xlsm, CA continuous: N=145
  out <- lr_test(obj_parent = 1633769.68, obj_child = 1491922.28,
                 df_parent = 5, df_child = 6, n = 145, response = "continuous")
  expect_equal(round(out$chi, 2), 13.17)
  expect_equal(out$df, 1)
  expect_equal(round(out$p, 4), 0.0003)
})

test_that("lr_test (binary) is the deviance difference", {
  # CA vs S/A binary: 184.424 - 183.5394 = 0.8846, df = 1
  out <- lr_test(obj_parent = 184.424, obj_child = 183.5394,
                 df_parent = 5, df_child = 6, n = 145, response = "binary")
  expect_equal(round(out$chi, 4), 0.8846)
  expect_equal(round(out$p, 3), 0.347)
})

test_that("select_parsimonious keeps the simplest model not significantly beaten", {
  # Reference not improved by any deviation -> pick reference
  fits <- list(
    reference = list(objective = 100, df = 5),
    SA        = list(objective = 99.9, df = 6),
    DR        = list(objective = 99.8, df = 7),
    DL        = list(objective = 99.7, df = 7)
  )
  pick <- select_parsimonious(fits, n = 145, response = "continuous", alpha = 0.05)
  expect_equal(pick, "reference")
})
