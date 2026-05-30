test_that("binary CA registry entries expose the right parameters", {
  ref <- model_spec("CA", "reference", n_chem = 2)
  expect_equal(ref$params, c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_null(ref$parent)
  expect_identical(ref$fn, ca_bi_vec)

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
  expect_identical(model_spec("IA", "reference", 2)$fn, ia_bi_vec)
  expect_identical(model_spec("IA", "DR", 2)$fn, ia_dr_bi_vec)
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
