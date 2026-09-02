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

test_that("fit_model optimises with a fixed deviation parameter (regression)", {
  # `a` is a deviation param. Fixing it must not corrupt parscale / perturbation:
  # pre-fix, indexing the free-keyed `ps`/`theta_i` with `spec$extra` (which still
  # contains the fixed `a`) caused a length mismatch that errored BOTH optimisers,
  # so fit_model returned the unchanged start with convergence == 99. Post-fix at
  # least one optimiser runs and moves the free base params off the start.
  df <- data.frame(C1 = c(0, 0.5, 1, 2, 0, 0, 0, 0.5, 1),
                   C2 = c(0, 0, 0, 0, 0.5, 1, 2, 0.5, 1),
                   Res = c(800, 600, 400, 200, 620, 420, 220, 480, 300))
  set.seed(7)
  # Deliberately-off start so a genuine optimiser run must move the free params.
  start <- c(max = 800, slope1 = 1, slope2 = 1, ec501 = 1, ec502 = 1, a = 0)
  fit <- fit_model(df, "CA", "SA", "continuous", start = start,
                   fixed = "a", n_starts = 2)

  # (a) fixed param honoured: `a` stays at its fixed start value.
  expect_equal(unname(fit$par[["a"]]), 0)

  # (b) the fit genuinely ran rather than falling back to the corrupted-state
  # result. Pre-fix both optimisers errored -> convergence == 99 (all-failed
  # fallback) with params unchanged.
  expect_false(identical(fit$convergence, 99L))

  # ...and the free base params actually moved away from the deliberately-off
  # start (impossible under the pre-fix unchanged-start fallback).
  moved <- any(abs(fit$par[c("slope1", "slope2", "ec501", "ec502")] -
                     start[c("slope1", "slope2", "ec501", "ec502")]) > 1e-6)
  expect_true(moved)
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

test_that("analyse_ternary uses a supplied base verbatim and skips the Stage-1 fit", {
  skip_on_cran()
  df <- utils::read.csv(testthat::test_path(
    "fixtures", "ternary", "fbsa_cpf_imi", "ternary_fbsa_cpf_imi_continuous.csv"))
  base <- c(max = 872.2, slope1 = 4.674, slope2 = 13, slope3 = 3.629,
            ec50_1 = 0.1275, ec50_2 = 5.58, ec50_3 = 0.575)
  res <- analyse_ternary(df, reference = "CA", response = "continuous", base = base)
  # base passed through unchanged (not re-fit from the singles)
  expect_equal(res$base[["max"]], 872.2)
  expect_equal(res$base[["ec50_1"]], 0.1275)
  expect_null(res$fits$base)                 # no Stage-1 fit object when base supplied
  expect_setequal(names(res$pairwise), c("A1", "A2", "A3"))
  expect_true(is.numeric(res$A4_overall))
  expect_gt(nrow(res$individual), 0)
})

test_that("analyse_ternary errors when the supplied base is missing a param", {
  df <- data.frame(C1 = c(0, 1, 1), C2 = c(0, 1, 0), C3 = c(0, 1, 1), Res = c(100, 20, 30))
  expect_error(analyse_ternary(df, base = c(max = 800, slope1 = 4)), "missing")
})

test_that("analyse_ternary returns overall + per-ratio structure", {
  truth <- list(max = 800, slope1 = 3, slope2 = 2, slope3 = 2.5,
                ec50_1 = 1, ec50_2 = 5, ec50_3 = 2, A1 = 0.5, A2 = 0, A3 = 0)
  doses <- c(0.5, 1, 2, 4)
  mk <- function(c1, c2, c3) data.frame(C1 = c1, C2 = c2, C3 = c3)
  df <- do.call(rbind, list(
    mk(0, 0, 0), mk(doses, 0, 0), mk(0, doses, 0), mk(0, 0, doses),
    mk(doses, doses * 5, 0), mk(doses, 0, doses * 2), mk(0, doses * 5, doses * 2),
    mk(doses, doses * 5, doses * 2)))
  df$Res <- ca_asa_tri_vec(df$C1, df$C2, df$C3, truth$max, truth$slope1,
                           truth$slope2, truth$slope3, truth$ec50_1,
                           truth$ec50_2, truth$ec50_3, truth$A1, truth$A2,
                           truth$A3, 0.7)
  res <- analyse_ternary(df, reference = "CA", response = "continuous",
                         n_starts = 3)
  expect_named(res, c("reference", "response", "base", "pairwise",
                      "A4_overall", "individual", "fits"))
  expect_equal(nrow(res$individual), 1)
  expect_true(is.numeric(res$A4_overall))

  eff <- ternary_effect_table(res, df)
  expect_true(all(c("ratio", "C1", "C2", "C3", "pred_CA", "pred_SA",
                    "pred_ASA", "a4_effect") %in% names(eff)))
  expect_equal(eff$a4_effect, eff$pred_ASA - eff$pred_SA, tolerance = 1e-8)
})

test_that("classify_rows handles a two-stressor frame: no ternary class", {
  two <- data.frame(C1 = c(0, 1, 0, 1), C2 = c(0, 0, 1, 1), Res = c(100, 80, 85, 60))
  cls <- classify_rows(two)
  expect_equal(as.character(cls), c("control", "single", "single", "binary"))
  expect_false("ternary" %in% as.character(cls))
})

test_that("supplying pairwise skips Stage 2 and holds A1/A2/A3 fixed", {
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  auto <- fit_ternary_asa(df, n_starts = 1)
  pw   <- auto$pairwise

  frozen <- fit_ternary_asa(df, base = auto$base, pairwise = pw, n_starts = 1)
  expect_equal(unname(frozen$pairwise), unname(pw), tolerance = 1e-10)
  expect_null(frozen$fits$pairwise)          # Stage 2 was not run
  expect_true(is.finite(frozen$A4_overall))
})

test_that("pairwise validates its names", {
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  expect_error(fit_ternary_asa(df, pairwise = c(A1 = 0, A2 = 0)),
               "missing param")
  expect_error(fit_ternary_asa(df, pairwise = list(A1 = 0, A2 = 0, A3 = 0)),
               "named numeric")
})

test_that("LINCHPIN: each pair's S/A a IS its ternary A-term at the same base", {
  # Spec section 3. If this fails, reusing the fitted binaries in the ternary
  # stage is no longer exact and the campaign design needs revisiting.
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  auto <- fit_ternary_asa(df, n_starts = 1)
  base <- auto$base

  # The binary base uses ec501/ec502; map the ternary base onto each pair
  # using the shipped production helper (not a private copy).
  a_of <- function(i, j) {
    b <- pair_base(base, i, j)
    fit_model(pair_df(df, i, j), "CA", "SA", "continuous",
              start = c(b, a = 0), fixed = names(b),
              n_starts = 1, time_limit = 60)$par[["a"]]
  }

  expect_equal(a_of(1, 2), unname(auto$pairwise[["A1"]]), tolerance = 1e-3)
  expect_equal(a_of(1, 3), unname(auto$pairwise[["A2"]]), tolerance = 1e-3)
  expect_equal(a_of(2, 3), unname(auto$pairwise[["A3"]]), tolerance = 1e-3)
})

test_that("tu_shares converts a mass ratio into toxic-unit shares", {
  z <- tu_shares(c(0.5, 0.5), c(0.1, 1))
  expect_equal(z, c(10 / 11, 1 / 11))
  expect_equal(sum(z), 1)
})

test_that("tu_shares leaves an equipotent mixture unchanged", {
  expect_equal(tu_shares(c(0.2, 0.3, 0.5), c(1, 1, 1)), c(0.2, 0.3, 0.5))
})

test_that("ternary_ratio_key groups identically in concentration and TU space", {
  # The point of issue #11: TU is a fixed per-column rescaling, so switching the
  # ratio to toxic units is a relabelling and never a regrouping.
  df <- data.frame(C1 = c(1, 2, 3, 1), C2 = c(2, 4, 6, 1), C3 = c(4, 8, 12, 1))
  ec <- c(0.1, 5, 0.7)
  conc_keys <- ternary_ratio_key(df)
  tu <- as.data.frame(t(apply(as.matrix(df), 1, function(p) p / ec)))
  names(tu) <- c("C1", "C2", "C3")
  expect_equal(as.integer(factor(conc_keys)),
               as.integer(factor(ternary_ratio_key(tu))))
})

test_that("ratio_label formats a proportion vector", {
  expect_equal(ratio_label(c(0.5875, 0.0219, 0.3906)), "0.59:0.02:0.39")
})
