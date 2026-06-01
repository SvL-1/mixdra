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

test_that("simulate_mixture continuous at cv=0 returns exact predictions", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d <- simulate_mixture(par, "CA", "reference", "continuous")  # cv = 0 default
  mu <- mixture_predict(d[c("C1", "C2")], par, "CA", "reference")
  expect_equal(d$Res, unname(mu))
})

test_that("simulate_mixture binary at group_size=Inf gives exact proportions", {
  par <- c(max = 0.95, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d <- simulate_mixture(par, "CA", "reference", "binary")  # group_size = Inf
  mu <- mixture_predict(d[c("C1", "C2")], par, "CA", "reference")
  expect_true(all(d$Exposed == 100))
  expect_equal(d$Affected / d$Exposed, unname(mu))   # exact, fractional Affected
})

test_that("simulate_mixture continuous noise scales with cv", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d0 <- simulate_mixture(par, "CA", "reference", "continuous", cv = 0)
  d1 <- simulate_mixture(par, "CA", "reference", "continuous",
                         cv = 0.1, reps = 200, seed = 1)
  mu1 <- mixture_predict(d1[c("C1", "C2")], par, "CA", "reference")
  rel <- (d1$Res - mu1) / abs(mu1)
  rel <- rel[is.finite(rel) & abs(mu1) > 1]   # ignore rows where mu ~ 0
  expect_equal(sd(rel), 0.1, tolerance = 0.03) # empirical CV near 0.1
})

test_that("simulate_mixture rejects binary max > 1", {
  par <- c(max = 5, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expect_error(simulate_mixture(par, "CA", "reference", "binary"),
               "must be <= 1")
})

test_that("simulate_mixture binary binomial sampling honours group_size", {
  par <- c(max = 0.9, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  d <- simulate_mixture(par, "CA", "reference", "binary",
                        group_size = 50, seed = 7)
  expect_true(all(d$Exposed == 50))
  expect_true(all(d$Affected == round(d$Affected)))      # integer counts
  expect_true(all(d$Affected >= 0 & d$Affected <= 50))
})
