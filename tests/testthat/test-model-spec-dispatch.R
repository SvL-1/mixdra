# model_spec(...)$fn is the single production seam. It must reproduce the legacy
# *_vec predictors for every dispatched (reference, deviation, n_chem) combo and
# leave ASA resolving to ca_asa_tri_vec.

test_that("binary dispatch reproduces legacy predictors", {
  m <- 800; sl <- c(6, 0.4); e <- c(0.08, 50); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))
  call_fn <- function(fn, extra = list())
    do.call(fn, c(list(c1 = g$c1, c2 = g$c2, max = m,
                       slope1 = sl[1], slope2 = sl[2],
                       ec501 = e[1], ec502 = e[2]), extra))

  expect_equal(unname(call_fn(model_spec("CA", "reference", 2)$fn)),
               unname(ca_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2])),
               tolerance = 1e-4)
  expect_equal(unname(call_fn(model_spec("CA", "DR", 2)$fn, list(a = a, b = 0.3))),
               unname(ca_dr_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, 0.3)),
               tolerance = 1e-4)
  expect_equal(unname(call_fn(model_spec("IA", "DL", 2)$fn, list(a = a, b = 0.2))),
               unname(ia_dl_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, 0.2)),
               tolerance = 1e-8)
})

test_that("ternary dispatch reproduces legacy reference/SA", {
  m <- 800; sl <- c(6, 0.4, 2); e <- c(0.08, 50, 10); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60), c3 = c(0, 5, 30))
  fn <- model_spec("IA", "SA", 3)$fn
  got <- fn(c1 = g$c1, c2 = g$c2, c3 = g$c3, max = m,
            slope1 = sl[1], slope2 = sl[2], slope3 = sl[3],
            ec50_1 = e[1], ec50_2 = e[2], ec50_3 = e[3], a = a)
  expect_equal(unname(got),
               unname(ia_sa_tri_vec(g$c1, g$c2, g$c3, m, sl[1], sl[2], sl[3],
                                    e[1], e[2], e[3], a)),
               tolerance = 1e-8)
})

test_that("ASA still resolves to the dedicated ternary predictor", {
  fn <- model_spec("CA", "ASA", 3)$fn
  expect_identical(fn, ca_asa_tri_vec)
})

test_that("dispatched params are unchanged", {
  expect_equal(model_spec("CA", "DR", 2)$params,
               c("max", "slope1", "slope2", "ec501", "ec502", "a", "b"))
  expect_equal(model_spec("CA", "ASA", 3)$extra, c("A1", "A2", "A3", "A4"))
})
