test_that("init_start_par defaults a=0, b=1, base=0; start overrides", {
  p <- init_start_par(
    params = c("max", "slope1", "slope2", "ec501", "ec502", "a", "b"),
    extra  = c("a", "b"),
    start  = c(max = 800, ec501 = 0.08))
  expect_equal(unname(p[["a"]]), 0)        # a starts at the solver origin 0
  expect_equal(unname(p[["b"]]), 1)        # b starts at 1
  expect_equal(unname(p[["slope1"]]), 0)   # base param not in start -> 0
  expect_equal(unname(p[["max"]]), 800)    # start overrides
  expect_equal(unname(p[["ec501"]]), 0.08) # start overrides
})

test_that("init_start_par sets ternary b1/b2/b3 to 1 and a to 0", {
  p <- init_start_par(
    params = c("max", "slope1", "slope2", "slope3",
               "ec50_1", "ec50_2", "ec50_3", "a", "b1", "b2", "b3"),
    extra  = c("a", "b1", "b2", "b3"),
    start  = numeric(0))
  expect_equal(unname(p[c("b1", "b2", "b3")]), c(1, 1, 1))
  expect_equal(unname(p[["a"]]), 0)
})
