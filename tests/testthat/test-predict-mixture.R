test_that("mixture_predict matches the verbatim binary predictor", {
  df <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  got <- mixture_predict(df, par, "CA", "reference")
  want <- ca_bi_vec(df$C1, df$C2, 800, 4, 1.5, 0.08, 1)
  expect_equal(unname(got), unname(want))
})

test_that("mixture_predict matches the verbatim ternary predictor", {
  df <- expand.grid(C1 = c(0, 0.1), C2 = c(0, 0.5), C3 = c(0, 2))
  par <- c(max = 100, slope1 = 2, slope2 = 1.5, slope3 = 3,
           ec50_1 = 0.1, ec50_2 = 0.5, ec50_3 = 2)
  got <- mixture_predict(df, par, "CA", "reference")
  want <- ca_tri_vec(df$C1, df$C2, df$C3, 100, 2, 1.5, 3, 0.1, 0.5, 2)
  expect_equal(unname(got), unname(want))
})

test_that("mixture_predict errors on a missing parameter", {
  df <- data.frame(C1 = c(0, 1), C2 = c(0, 1))
  expect_error(
    mixture_predict(df, c(max = 1, slope1 = 1, slope2 = 1, ec501 = 1), "CA", "reference"),
    "missing parameter"
  )
})
