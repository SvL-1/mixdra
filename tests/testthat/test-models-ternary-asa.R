test_that("ca_asa_tri reduces to ca_tri when all A = 0", {
  args <- list(c1 = 0.05, c2 = 10, c3 = 2, max = 800,
               slope1 = 6, slope2 = 0.4, slope3 = 1,
               ec50_1 = 0.08, ec50_2 = 50, ec50_3 = 10)
  ref <- do.call(ca_tri, args)
  got <- do.call(ca_asa_tri, c(args, A1 = 0, A2 = 0, A3 = 0, A4 = 0))
  expect_equal(got, ref, tolerance = 1e-6)
})

test_that("ca_asa_tri binary (C3=0) branch equals ca_sa_tri with a = A1", {
  base <- list(max = 800, slope1 = 6, slope2 = 0.4, slope3 = 1,
               ec50_1 = 0.08, ec50_2 = 50, ec50_3 = 10)
  sa  <- do.call(ca_sa_tri,  c(list(c1 = 0.05, c2 = 10, c3 = 0), base, a = 0.7))
  asa <- do.call(ca_asa_tri, c(list(c1 = 0.05, c2 = 10, c3 = 0), base,
                               A1 = 0.7, A2 = 0, A3 = 0, A4 = 0))
  expect_equal(asa, sa, tolerance = 1e-6)
})

test_that("ca_asa_tri three-way branch with only A4 equals ca_sa_tri a = A4", {
  base <- list(max = 800, slope1 = 6, slope2 = 0.4, slope3 = 1,
               ec50_1 = 0.08, ec50_2 = 50, ec50_3 = 10)
  sa  <- do.call(ca_sa_tri,  c(list(c1 = 0.05, c2 = 10, c3 = 2), base, a = 0.9))
  asa <- do.call(ca_asa_tri, c(list(c1 = 0.05, c2 = 10, c3 = 2), base,
                               A1 = 0, A2 = 0, A3 = 0, A4 = 0.9))
  expect_equal(asa, sa, tolerance = 1e-6)
})

test_that("ca_asa_tri_vec maps over concentration vectors", {
  out <- ca_asa_tri_vec(c(0, 0.05), c(0, 10), c(0, 2), 800, 6, 0.4, 1,
                        0.08, 50, 10, A1 = 0.5, A2 = 0.1, A3 = -0.2, A4 = 0.3)
  expect_length(out, 2)
  expect_equal(out[1], 800)  # control row -> Max
})
