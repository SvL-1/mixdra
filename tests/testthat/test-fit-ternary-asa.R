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

test_that("fit_ternary_asa recovers known base, A1-A3, and per-ratio A4", {
  truth <- list(max = 800, slope1 = 3, slope2 = 2, slope3 = 2.5,
                ec50_1 = 1, ec50_2 = 5, ec50_3 = 2,
                A1 = 0.6, A2 = -0.4, A3 = 0.3)
  doses <- c(0.25, 0.5, 1, 2, 4)

  mk <- function(c1, c2, c3) data.frame(C1 = c1, C2 = c2, C3 = c3)
  rows <- list(
    mk(0, 0, 0),
    mk(doses, 0, 0), mk(0, doses, 0), mk(0, 0, doses),
    mk(doses, doses * 5, 0), mk(doses, 0, doses * 2), mk(0, doses * 5, doses * 2),
    mk(doses, doses * 5, doses * 2),
    mk(doses * 2, doses, doses))
  df <- do.call(rbind, rows)

  key <- rep("", nrow(df))
  tern <- classify_rows(df) == "ternary"
  key[tern] <- ternary_ratio_key(df[tern, , drop = FALSE])
  uk <- unique(key[tern]); A4_true <- c(0.8, -0.5); names(A4_true) <- uk

  Res <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    a4 <- if (key[i] == "") 0 else A4_true[[key[i]]]
    Res[i] <- ca_asa_tri(df$C1[i], df$C2[i], df$C3[i], truth$max,
                         truth$slope1, truth$slope2, truth$slope3,
                         truth$ec50_1, truth$ec50_2, truth$ec50_3,
                         truth$A1, truth$A2, truth$A3, a4)
  }
  df$Res <- Res

  set.seed(1)
  fit <- fit_ternary_asa(df, reference = "CA", response = "continuous",
                         n_starts = 4)

  expect_equal(unname(fit$base[["max"]]), 800, tolerance = 0.02)
  expect_equal(unname(fit$base[["ec50_1"]]), 1, tolerance = 0.05)
  expect_equal(unname(fit$pairwise[["A1"]]), 0.6, tolerance = 0.1)
  expect_equal(unname(fit$pairwise[["A3"]]), 0.3, tolerance = 0.1)
  ind <- fit$individual
  expect_equal(nrow(ind), 2)
  expect_gt(ind$A4[ind$ratio == uk[1]], 0.3)   # synergy
  expect_lt(ind$A4[ind$ratio == uk[2]], -0.1)  # antagonism
})
