test_that("obj_ss returns the sum of squared residuals", {
  expect_equal(obj_ss(c(1, 2, 3), c(1, 2, 3)), 0)
  expect_equal(obj_ss(c(1, 2, 3), c(2, 2, 2)), 1 + 0 + 1)
})

test_that("binlik per-row contribution matches the Jonker/VBA formula", {
  # P*log(pi_hat/pi) + (T-P)*log((1-pi_hat)/(1-pi)); pi = P/T
  # T=10, P=8 (pi=0.8), pi_hat=0.9656 -> matches workbook row (Contr.L = -2.0154)
  expect_equal(binlik(exposed = 10, affected = 8, pi_hat = 0.9656),
               8 * log(0.9656 / 0.8) + 2 * log((1 - 0.9656) / (1 - 0.8)),
               tolerance = 1e-9)
  expect_equal(round(binlik(10, 8, 0.9656), 4), -2.0154)
  # Saturated cell (pi = 1) contributes only the positive term
  expect_equal(binlik(10, 10, 0.9656), 10 * log(0.9656 / 1), tolerance = 1e-9)
})

test_that("obj_deviance is -2 * sum of per-row binlik contributions", {
  exposed <- c(10, 10, 10)
  affected <- c(10, 8, 10)
  pi_hat <- c(0.9656, 0.9656, 0.9656)
  expected <- -2 * sum(mapply(binlik, exposed, affected, pi_hat))
  expect_equal(obj_deviance(exposed, affected, pi_hat), expected)
})
