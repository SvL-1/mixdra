test_that("ca_tri reduces to a single-chemical curve when two concs are 0", {
  expect_equal(ca_tri(c1 = 0.05, c2 = 0, c3 = 0, max = 800,
                      slope1 = 6, slope2 = 0.4, slope3 = 1,
                      ec50_1 = 0.08, ec50_2 = 100, ec50_3 = 10),
               800 / (1 + (0.05 / 0.08)^6))
})

test_that("ca_tri reduces to ca_bi when the third conc is 0", {
  expect_equal(
    ca_tri(0.05, 10, 0, 800, 6, 0.4, 1, 0.08, 50, 10),
    ca_bi(0.05, 10, 800, 6, 0.4, 0.08, 50),
    tolerance = 1e-4)
})

test_that("ia_tri matches the triple independent-action product", {
  f1 <- 1 / (1 + (0.05 / 0.08)^6)
  f2 <- 1 / (1 + (10 / 50)^0.4)
  f3 <- 1 / (1 + (2 / 10)^1)
  expect_equal(ia_tri(0.05, 10, 2, 800, 6, 0.4, 1, 0.08, 50, 10),
               800 * f1 * f2 * f3)
})

test_that("vectorised ternary wrapper maps over three concentration vectors", {
  out <- ca_tri_vec(c(0, 0.05), c(0, 0), c(0, 0), 800, 6, 0.4, 1, 0.08, 100, 10)
  expect_length(out, 2)
  expect_equal(out[1], 800)
})
