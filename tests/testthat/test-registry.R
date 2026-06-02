test_that("binary CA registry entries expose the right parameters", {
  ref <- model_spec("CA", "reference", n_chem = 2)
  expect_equal(ref$params, c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_null(ref$parent)
  # The 16 reference/SA/DR/DL models now dispatch through make_adapter (a closure
  # over mix_response), not the legacy *_vec functions. The contract is numeric
  # equivalence with the legacy predictor, not function identity.
  expect_equal(
    unname(ref$fn(c1 = c(0, 0.1), c2 = c(0, 20), max = 800,
                  slope1 = 6, slope2 = 0.4, ec501 = 0.08, ec502 = 50)),
    unname(ca_bi_vec(c(0, 0.1), c(0, 20), 800, 6, 0.4, 0.08, 50)),
    tolerance = 1e-4)

  sa <- model_spec("CA", "SA", n_chem = 2)
  expect_equal(sa$params, c("max", "slope1", "slope2", "ec501", "ec502", "a"))
  expect_equal(sa$parent, "reference")
  expect_equal(sa$extra, "a")

  dr <- model_spec("CA", "DR", n_chem = 2)
  expect_equal(dr$extra, c("a", "b"))
  expect_equal(dr$parent, "SA")

  dl <- model_spec("CA", "DL", n_chem = 2)
  expect_equal(dl$extra, c("a", "b"))
  expect_equal(dl$parent, "SA")
})

test_that("binary IA registry resolves to the IA predictor functions", {
  # Now dispatched through make_adapter; assert numeric equivalence with the
  # legacy IA predictors rather than function identity.
  ref <- model_spec("IA", "reference", 2)$fn
  expect_equal(
    unname(ref(c1 = c(0, 0.1), c2 = c(0, 20), max = 800,
               slope1 = 6, slope2 = 0.4, ec501 = 0.08, ec502 = 50)),
    unname(ia_bi_vec(c(0, 0.1), c(0, 20), 800, 6, 0.4, 0.08, 50)),
    tolerance = 1e-8)

  dr <- model_spec("IA", "DR", 2)$fn
  expect_equal(
    unname(dr(c1 = c(0, 0.1), c2 = c(0, 20), max = 800,
              slope1 = 6, slope2 = 0.4, ec501 = 0.08, ec502 = 50,
              a = 0.5, b = 0.3)),
    unname(ia_dr_bi_vec(c(0, 0.1), c(0, 20), 800, 6, 0.4, 0.08, 50, 0.5, 0.3)),
    tolerance = 1e-4)
})

test_that("model_spec resolves ternary Advanced S/A (ASA)", {
  spec <- model_spec("CA", "ASA", 3)
  expect_equal(spec$params,
               c("max", "slope1", "slope2", "slope3",
                 "ec50_1", "ec50_2", "ec50_3", "A1", "A2", "A3", "A4"))
  expect_equal(spec$extra, c("A1", "A2", "A3", "A4"))
  expect_identical(spec$fn, ca_asa_tri_vec)
})

test_that("model_spec rejects ASA for binary", {
  expect_error(model_spec("CA", "ASA", 2), "ASA")
})
