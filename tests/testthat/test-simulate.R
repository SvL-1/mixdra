test_that("mixture_design builds control + single ladders + rays (binary)", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d <- mixture_design(par, "CA", "reference")
  expect_named(d, c("C1", "C2"))
  # exactly one control row (all zero)
  expect_equal(sum(d$C1 == 0 & d$C2 == 0), 1L)
  # single-C1 ladder present: C2 == 0 and 7 distinct positive C1 values
  c1_only <- d[d$C2 == 0 & d$C1 > 0, ]
  expect_equal(nrow(c1_only), 7L)
  # ladder spans 1/8 .. 8 of ec501
  expect_equal(min(c1_only$C1), 0.08 * 2^-3, tolerance = 1e-9)
  expect_equal(max(c1_only$C1), 0.08 * 2^3, tolerance = 1e-9)
  # mixture rays present (both chemicals positive on some rows)
  expect_gt(sum(d$C1 > 0 & d$C2 > 0), 0L)
})

test_that("mixture_design builds a 3-column grid for ternary par", {
  par <- c(max = 100, slope1 = 2, slope2 = 1.5, slope3 = 3,
           ec50_1 = 0.1, ec50_2 = 0.5, ec50_3 = 2)
  d <- mixture_design(par, "CA", "reference")
  expect_named(d, c("C1", "C2", "C3"))
  expect_equal(sum(d$C1 == 0 & d$C2 == 0 & d$C3 == 0), 1L)   # one control
  expect_equal(nrow(d[d$C2 == 0 & d$C3 == 0 & d$C1 > 0, ]), 7L)  # C1 ladder
})

test_that("mixture_design rejects par without 2 or 3 slopes", {
  expect_error(mixture_design(c(max = 1, slope1 = 1, ec501 = 1)),
               "2 or 3 slope")
})
