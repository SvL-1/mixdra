# A minimal analyse_ternary-shaped result for unit testing the pure isoplane API
# (no fitting needed — the functions only read these fields).
mock_res <- function(A1 = 0, A2 = 0, A3 = 0, A4 = 0,
                     ec50 = c(2, 5, 0.5), max = 100) {
  list(
    reference = "CA", response = "continuous",
    base = c(max = max, slope1 = 3, slope2 = 3, slope3 = 3,
             ec50_1 = ec50[1], ec50_2 = ec50[2], ec50_3 = ec50[3]),
    pairwise = c(A1 = A1, A2 = A2, A3 = A3),
    A4_overall = A4,
    individual = data.frame(
      ratio = c("a", "b"), C1 = c(1/3, 0.6), C2 = c(1/3, 0.2),
      C3 = c(1/3, 0.2), A4 = c(0.5, -0.3), n = c(5L, 5L),
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
