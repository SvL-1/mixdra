test_that("ca_bi reduces to the single-chemical curve when one conc is 0", {
  # C2 = 0 -> Y = max / (1 + (C1/ec501)^slope1)
  expect_equal(ca_bi(c1 = 0.05, c2 = 0, max = 800, slope1 = 6, slope2 = 0.4,
                     ec501 = 0.08, ec502 = 100),
               800 / (1 + (0.05 / 0.08)^6))
  # Both zero, decreasing endpoint (positive slopes) -> control = max
  expect_equal(ca_bi(0, 0, 800, 6, 0.4, 0.08, 100), 800)
})

test_that("ia_bi matches the closed-form independent-action product", {
  f1 <- 1 / (1 + (0.05 / 0.08)^6)
  f2 <- 1 / (1 + (10 / 50)^0.4)
  expect_equal(ia_bi(0.05, 10, 800, 6, 0.4, 0.08, 50), 800 * f1 * f2)
})

test_that("ca_sa_bi reduces to ca_bi when a = 0", {
  expect_equal(ca_sa_bi(0.05, 10, 800, 6, 0.4, 0.08, 50, a = 0),
               ca_bi(0.05, 10, 800, 6, 0.4, 0.08, 50),
               tolerance = 1e-5)
})

test_that("vectorised wrapper maps over concentration vectors", {
  out <- ca_bi_vec(c1 = c(0, 0.05), c2 = c(0, 0),
                   max = 800, slope1 = 6, slope2 = 0.4, ec501 = 0.08, ec502 = 100)
  expect_length(out, 2)
  expect_equal(out[1], 800)
})
