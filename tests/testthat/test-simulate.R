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

test_that("mixture_design ray concentrations equal dose * w * ec50", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d <- mixture_design(par, "CA", "reference")
  ec50 <- c(0.08, 1)
  # The 1:1 ray: equal toxic-unit fractions w = c(0.5, 0.5). For every such row
  # the toxic-unit ratio C1/ec501 must equal C2/ec502, and the top dose point
  # corresponds to dose = 2^3 in toxic units (split across the two chemicals).
  mix <- d[d$C1 > 0 & d$C2 > 0, ]
  tu1 <- mix$C1 / ec50[1]
  tu2 <- mix$C2 / ec50[2]
  # rows from the 1:1 ray have equal toxic units in each chemical
  on_11 <- abs(tu1 - tu2) < 1e-9
  expect_gt(sum(on_11), 0L)
  # on the 1:1 ray, total toxic units span 2^-3 .. 2^3 (dose series)
  total_tu <- (tu1 + tu2)[on_11]
  expect_equal(min(total_tu), 2^-3, tolerance = 1e-9)
  expect_equal(max(total_tu), 2^3, tolerance = 1e-9)
})

test_that("mixture_design errors when an EC50 is missing from par", {
  expect_error(
    mixture_design(c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08)),
    "must supply every EC50"
  )
})
