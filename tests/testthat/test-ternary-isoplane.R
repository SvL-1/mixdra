# A minimal analyse_ternary-shaped result for unit testing the pure isoplane API
# (no fitting needed — the functions only read these fields).
mock_res <- function(A1 = 0, A2 = 0, A3 = 0, A4 = 0,
                     ec50 = c(2, 5, 0.5), max = 100,
                     A4_ind = c(0.5, -0.3)) {
  list(
    reference = "CA", response = "continuous",
    base = c(max = max, slope1 = 3, slope2 = 3, slope3 = 3,
             ec50_1 = ec50[1], ec50_2 = ec50[2], ec50_3 = ec50[3]),
    pairwise = c(A1 = A1, A2 = A2, A3 = A3),
    A4_overall = A4,
    individual = data.frame(
      ratio = c("a", "b"), C1 = c(1/3, 0.6), C2 = c(1/3, 0.2),
      C3 = c(1/3, 0.2), A4 = A4_ind, n = c(5L, 5L),
      stringsAsFactors = FALSE))
}

test_that("simplex_grid has (n+1)(n+2)/2 rows summing to 1", {
  g <- simplex_grid(4)
  expect_equal(nrow(g), (4 + 1) * (4 + 2) / 2)   # 15
  expect_equal(rowSums(g[c("z1", "z2", "z3")]), rep(1, nrow(g)), tolerance = 1e-12)
  expect_true(all(g >= 0))
})

test_that("ec50_isoplane with no interaction is the pure-CA reference", {
  res <- mock_res()                       # all A = 0
  d <- ec50_isoplane(res, "SA", n = 10)
  expect_equal(d$sigma_tu, rep(1, nrow(d)), tolerance = 1e-12)
  # C_i = EC50_i * z_i when F4 = 1
  expect_equal(d$C1, res$base[["ec50_1"]] * d$z1, tolerance = 1e-12)
  expect_equal(d$C2, res$base[["ec50_2"]] * d$z2, tolerance = 1e-12)
  expect_equal(d$C3, res$base[["ec50_3"]] * d$z3, tolerance = 1e-12)
})

test_that("ec50_isoplane vertices give single-chemical EC50 points", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 2)
  d <- ec50_isoplane(res, "ASA", n = 6)
  v1 <- d[d$z1 == 1, ]                     # vertex z = (1,0,0)
  expect_equal(nrow(v1), 1)
  expect_equal(v1$sigma_tu, 1, tolerance = 1e-12)        # F4 = exp(0) = 1
  expect_equal(v1$C1, res$base[["ec50_1"]], tolerance = 1e-12)
  expect_equal(c(v1$C2, v1$C3), c(0, 0), tolerance = 1e-12)
})

test_that("ec50_isoplane SA ignores A4, ASA uses A4_overall", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 5)
  sa  <- ec50_isoplane(res, "SA", n = 8)
  asa <- ec50_isoplane(res, "ASA", n = 8)
  # interior point with all z > 0 differs between SA and ASA; edges (a z == 0) match
  interior <- sa$z1 > 0 & sa$z2 > 0 & sa$z3 > 0
  expect_false(isTRUE(all.equal(sa$sigma_tu[interior], asa$sigma_tu[interior])))
  expect_equal(sa$sigma_tu[!interior], asa$sigma_tu[!interior], tolerance = 1e-12)
})

test_that("sigma_tu_curve: no interaction gives sigma_tu == 1 everywhere", {
  res <- mock_res()                       # all A = 0
  d <- sigma_tu_curve(res, "SA", n = 11)
  expect_setequal(unique(d$chem), c("C1", "C2", "C3"))
  expect_equal(d$sigma_tu, rep(1, nrow(d)), tolerance = 1e-12)
})

test_that("sigma_tu_curve: each chemical's vertex (z = 1) has sigma_tu == 1", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 2)
  d <- sigma_tu_curve(res, "ASA", n = 11)
  tip <- d[d$z == 1, ]
  expect_equal(nrow(tip), 3)              # one per chemical
  expect_equal(tip$sigma_tu, rep(1, 3), tolerance = 1e-12)
})

test_that("sigma_tu_curve: C1 at z = 0 is the equal-split FBSA/IMI binary point", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 2)
  d <- sigma_tu_curve(res, "SA", n = 11)
  c1_0 <- d[d$chem == "C1" & d$z == 0, "sigma_tu"]
  # other two equal at z = 0.5 each => F4 = exp(A3 * 0.5 * 0.5)
  expect_equal(c1_0, exp(res$pairwise[["A3"]] * 0.25), tolerance = 1e-12)
})

test_that("sigma_tu_curve SA ignores A4, ASA uses A4_overall; model column written correctly", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 5)
  sa  <- sigma_tu_curve(res, "SA",  n = 11)
  asa <- sigma_tu_curve(res, "ASA", n = 11)
  # model column is populated correctly
  expect_true(all(sa$model  == "SA"))
  expect_true(all(asa$model == "ASA"))
  # interior points (focal z strictly between 0 and 1, so the other two are > 0)
  # => z1*z2*z3 > 0 => three-way term active => SA != ASA
  interior <- sa$z > 0 & sa$z < 1
  expect_false(isTRUE(all.equal(sa$sigma_tu[interior], asa$sigma_tu[interior])))
  # endpoints (z == 0 or z == 1): at least one of the other z's is 0 => term vanishes
  expect_equal(sa$sigma_tu[!interior], asa$sigma_tu[!interior], tolerance = 1e-12)
})

test_that("ec50_markers returns one row per individual ratio", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8)
  m <- ec50_markers(res, df = NULL)
  expect_equal(nrow(m), nrow(res$individual))
  expect_equal(m$ratio, res$individual$ratio)
  expect_true(all(is.finite(m$sigma_tu)))
  expect_true(all(c("z1", "z2", "z3", "C1", "C2", "C3", "sigma_tu") %in% names(m)))
})

test_that("ec50_markers TU-fraction conversion: equal proportions are NOT equal z", {
  # ratio "a" has equal concentration proportions (1/3 each) but unequal EC50s,
  # so its z (TU fractions) must be unequal and weighted by 1/EC50_i.
  res <- mock_res(A1 = 0, A2 = 0, A3 = 0, ec50 = c(2, 5, 0.5), A4_ind = c(0, 0))
  m <- ec50_markers(res, df = NULL)
  ra <- m[m$ratio == "a", ]
  ec <- c(2, 5, 0.5); z_expected <- (1 / ec) / sum(1 / ec)
  expect_equal(c(ra$z1, ra$z2, ra$z3), z_expected, tolerance = 1e-12)
  expect_equal(ra$sigma_tu, 1, tolerance = 1e-12)   # all A = 0 here
})
