# make_adapter() must reproduce the legacy *_vec predictors exactly (same maths,
# same by-name call contract). The legacy functions are the equivalence oracle in
# tests/testthat/helper-legacy-models.R (auto-sourced by testthat).

test_that("make_adapter reproduces binary CA reference/SA over a grid", {
  m <- 800; sl <- c(6, 0.4); e <- c(0.08, 50); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))

  fn_ref <- make_adapter("CA", "reference", 2)
  got <- fn_ref(c1 = g$c1, c2 = g$c2, max = m,
                slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2])
  want <- ca_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2])
  expect_equal(unname(got), unname(want), tolerance = 1e-4)

  fn_sa <- make_adapter("CA", "SA", 2)
  got <- fn_sa(c1 = g$c1, c2 = g$c2, max = m,
               slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2], a = a)
  want <- ca_sa_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a)
  expect_equal(unname(got), unname(want), tolerance = 1e-4)
})

test_that("make_adapter maps binary DR b -> c(b, 0) and DL b scalar", {
  m <- 800; sl <- c(6, 0.4); e <- c(0.08, 50); a <- 0.5; b <- 0.3
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))

  fn_dr <- make_adapter("IA", "DR", 2)
  got <- fn_dr(c1 = g$c1, c2 = g$c2, max = m,
               slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2], a = a, b = b)
  want <- ia_dr_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, b)
  expect_equal(unname(got), unname(want), tolerance = 1e-8)

  fn_dl <- make_adapter("IA", "DL", 2)
  got <- fn_dl(c1 = g$c1, c2 = g$c2, max = m,
               slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2], a = a, b = 0.2)
  want <- ia_dl_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, 0.2)
  expect_equal(unname(got), unname(want), tolerance = 1e-8)
})

test_that("make_adapter reproduces ternary reference/SA", {
  m <- 800; sl <- c(6, 0.4, 2); e <- c(0.08, 50, 10); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60), c3 = c(0, 5, 30))

  fn <- make_adapter("CA", "SA", 3)
  got <- fn(c1 = g$c1, c2 = g$c2, c3 = g$c3, max = m,
            slope1 = sl[1], slope2 = sl[2], slope3 = sl[3],
            ec50_1 = e[1], ec50_2 = e[2], ec50_3 = e[3], a = a)
  want <- ca_sa_tri_vec(g$c1, g$c2, g$c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3], a)
  expect_equal(unname(got), unname(want), tolerance = 1e-4)
})
