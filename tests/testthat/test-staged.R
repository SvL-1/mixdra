# Synthetic binary data with a full single-compound dose series per chemical
# (>= 6 distinct marginal concentrations) so stage 1 can identify each curve,
# plus a few interior mixture points. Noiseless CA-reference surface.
make_staged_binary_df <- function(truth = list(max = 800, slope1 = 4,
                                               slope2 = 1.5, ec501 = 0.08,
                                               ec502 = 1)) {
  c1 <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
  c2 <- c(0.125, 0.25, 0.5, 1, 2, 4)
  singles <- rbind(data.frame(C1 = c(0, c1), C2 = 0),
                   data.frame(C1 = 0,        C2 = c2))
  mix <- expand.grid(C1 = c(0.04, 0.08), C2 = c(0.5, 1))
  d <- unique(rbind(singles, mix))
  d$Res <- ca_bi_vec(d$C1, d$C2, truth$max, truth$slope1, truth$slope2,
                     truth$ec501, truth$ec502)
  d
}

test_that("fit_curve_from_singles recovers curve params from the marginals", {
  d <- make_staged_binary_df()
  base <- fit_curve_from_singles(d, "CA", "continuous", n_starts = 1)
  expect_setequal(names(base), c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_equal(unname(base[["max"]]),   800,  tolerance = 1e-2)
  expect_equal(unname(base[["ec501"]]), 0.08, tolerance = 1e-2)
  expect_equal(unname(base[["ec502"]]), 1,    tolerance = 1e-2)
  expect_equal(unname(base[["slope1"]]), 4,   tolerance = 1e-2)
  expect_equal(unname(base[["slope2"]]), 1.5, tolerance = 1e-2)
})
