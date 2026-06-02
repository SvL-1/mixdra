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

test_that("compare_fits builds the LR comparison frame and chooses a model", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  seed <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fits <- lapply(c("reference", "SA", "DR", "DL"), function(d)
    fit_model(df, "CA", d, "continuous", start = seed, fixed = names(seed)))
  names(fits) <- c("reference", "SA", "DR", "DL")

  cmp <- compare_fits(fits, n = nrow(df), response = "continuous", alpha = 0.05)
  expect_true(all(c("model", "parent", "chi", "df", "p") %in% names(cmp$comparison)))
  expect_equal(nrow(cmp$comparison), 3)          # SA, DR, DL each vs their parent
  expect_true(cmp$chosen %in% c("reference", "SA", "DR", "DL"))
})

test_that("compare_fits handles a ternary (reference + SA only) fit set", {
  fits <- list(
    reference = list(objective = 10, df = 7),
    SA        = list(objective = 9,  df = 8))
  cmp <- compare_fits(fits, n = 27, response = "continuous", alpha = 0.05)
  expect_equal(nrow(cmp$comparison), 1)          # only SA vs reference
  expect_equal(cmp$comparison$model, "SA")
  expect_equal(cmp$chosen, "reference")          # p=0.092 > 0.05; SA not significant
})

test_that("compare_fits drops a child whose parent is absent", {
  fits <- list(
    reference = list(objective = 10, df = 5),
    DR        = list(objective = 8,  df = 7))   # DR's parent SA is absent
  cmp <- compare_fits(fits, n = 20, response = "continuous", alpha = 0.05)
  expect_null(cmp$comparison)                    # no comparable child -> NULL
  expect_equal(cmp$chosen, "reference")          # walk stops at reference
})
