test_that("upload_schema returns the fixed columns per stage and response", {
  expect_equal(upload_schema("single", "continuous"), c("Conc", "Res"))
  expect_equal(upload_schema("single", "quantal"), c("Conc", "Affected", "Exposed"))
  expect_equal(upload_schema("binary", "continuous"), c("C1", "C2", "Res"))
  expect_equal(upload_schema("binary", "quantal"), c("C1", "C2", "Affected", "Exposed"))
})

test_that("template_df has exactly the schema columns and at least one example row", {
  d <- template_df("binary", "continuous")
  expect_equal(names(d), c("C1", "C2", "Res"))
  expect_gt(nrow(d), 0)
  expect_true(all(vapply(d, is.numeric, logical(1))))
  # binary template includes single-chemical rows (one chem at 0) for seeding
  expect_true(any(d$C2 == 0 & d$C1 > 0))
  expect_true(any(d$C1 == 0 & d$C2 > 0))
})

test_that("validate_upload returns no errors for a valid file", {
  good <- data.frame(C1 = c(0, 1, 0, 2), C2 = c(0, 0, 1, 2), Res = c(100, 50, 60, 20))
  expect_length(validate_upload(good, "binary", "continuous"), 0)
})

test_that("validate_upload reports missing columns", {
  bad <- data.frame(C1 = c(0, 1), Res = c(100, 50))   # no C2
  errs <- validate_upload(bad, "binary", "continuous")
  expect_match(paste(errs, collapse = " "), "C2")
})

test_that("validate_upload rejects non-numeric, negative conc, and Affected > Exposed", {
  expect_match(paste(validate_upload(
    data.frame(C1 = c("a", "b"), C2 = c(0, 1), Res = c(1, 2)),
    "binary", "continuous"), collapse = " "), "[Nn]on-numeric")
  expect_match(paste(validate_upload(
    data.frame(C1 = c(-1, 1), C2 = c(0, 1), Res = c(1, 2)),
    "binary", "continuous"), collapse = " "), ">= 0|negative|0")
  expect_match(paste(validate_upload(
    data.frame(C1 = c(0, 1), C2 = c(0, 1), Affected = c(2, 12), Exposed = c(10, 10)),
    "binary", "quantal"), collapse = " "), "Affected")
})

test_that("validate_upload needs >= 4 distinct concentrations for a single fit", {
  short <- data.frame(Conc = c(0, 1, 2), Res = c(100, 50, 10))
  expect_match(paste(validate_upload(short, "single", "continuous"), collapse = " "),
               "distinct")
})

test_that("read_upload reads a CSV file into a data frame", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Conc = c(0, 1), Res = c(100, 50)), path, row.names = FALSE)
  d <- read_upload(path)
  expect_equal(names(d), c("Conc", "Res"))
  expect_equal(nrow(d), 2)
})

test_that("to_engine_df renames Conc to C1 for single only", {
  s <- to_engine_df(data.frame(Conc = c(0, 1), Res = c(1, 2)), "single")
  expect_true("C1" %in% names(s))
  expect_false("Conc" %in% names(s))
  b <- to_engine_df(data.frame(C1 = 0, C2 = 1, Res = 3), "binary")
  expect_equal(names(b), c("C1", "C2", "Res"))
})

test_that("collect_bounds keeps only supplied bounds, named by parameter", {
  vals <- list(lo_max = NA, hi_max = 1, lo_slope1 = NA, hi_slope1 = NA,
               lo_slope2 = NA, hi_slope2 = NA, lo_ec501 = 0.01, hi_ec501 = NA,
               lo_ec502 = NA, hi_ec502 = NA)
  b <- collect_bounds(vals)
  expect_equal(b$upper, c(max = 1))
  expect_equal(b$lower, c(ec501 = 0.01))
})

test_that("collect_bounds returns NULL bounds when nothing supplied", {
  vals <- setNames(as.list(rep(NA, 10)),
                   c(paste0("lo_", c("max","slope1","slope2","ec501","ec502")),
                     paste0("hi_", c("max","slope1","slope2","ec501","ec502"))))
  b <- collect_bounds(vals)
  expect_null(b$lower)
  expect_null(b$upper)
})

test_that("collect_bounds reads a custom parameter set (single tab)", {
  vals <- list(lo_max = NA, hi_max = NA, lo_slope = NA, hi_slope = 1.5,
               lo_ec50 = 0.01, hi_ec50 = NA)
  b <- collect_bounds(vals, c("max", "slope", "ec50"))
  expect_equal(b$upper, c(slope = 1.5))
  expect_equal(b$lower, c(ec50 = 0.01))
})

test_that("marginal_df extracts a chemical's single series, renaming its conc to C1", {
  df <- data.frame(C1 = c(0, 1, 2, 0, 0, 3),
                   C2 = c(0, 0, 0, 1, 2, 4),
                   Res = c(100, 60, 40, 70, 50, 10))

  m1 <- marginal_df(df, 1)               # rows where C2 == 0
  expect_equal(m1$C1, c(0, 1, 2))
  expect_equal(m1$Res, c(100, 60, 40))
  expect_false("C2" %in% names(m1))

  m2 <- marginal_df(df, 2)               # rows where C1 == 0; C2 renamed to C1
  expect_equal(m2$C1, c(0, 1, 2))
  expect_equal(m2$Res, c(100, 70, 50))
  expect_false("C2" %in% names(m2))
})

test_that("assemble_curve_params averages max and keeps per-chemical slope/ec50", {
  f1 <- list(par = c(max = 700, slope = 2, ec50 = 1))
  f2 <- list(par = c(max = 600, slope = 1, ec50 = 5))
  p <- assemble_curve_params(f1, f2)
  expect_equal(names(p), c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_equal(unname(p[["max"]]),    650)   # mean(700, 600)
  expect_equal(unname(p[["slope1"]]), 2)
  expect_equal(unname(p[["slope2"]]), 1)
  expect_equal(unname(p[["ec501"]]),  1)
  expect_equal(unname(p[["ec502"]]),  5)
})

test_that("collect_bounds_all reads lower/upper for all seven binary params", {
  vals <- list(olo_max = 100, ohi_max = 1000,
               olo_a = 0,   ohi_a = 5,
               olo_b = NA,  ohi_b = NA)
  b <- collect_bounds_all(vals)
  expect_equal(b$lower[["max"]], 100)
  expect_equal(b$upper[["max"]], 1000)
  expect_equal(b$lower[["a"]], 0)
  expect_equal(b$upper[["a"]], 5)
  expect_false("b" %in% names(b$lower))   # blank dropped
})

test_that("split_fixed_bounds turns equal lower/upper into a fixed param", {
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.5)
  sp <- split_fixed_bounds(lower = c(max = 800, a = 0),
                           upper = c(max = 800, a = 5), start = start)
  expect_true("max" %in% sp$fixed)     # equal bounds -> fixed
  expect_false("a" %in% sp$fixed)      # a is a true range, not pinned
  expect_equal(sp$start[["max"]], 800) # pinned at the equal-bound value
  expect_null(sp$lower[["max"]])       # pinned param dropped from bounds
  expect_equal(sp$lower[["a"]], 0)
  expect_equal(sp$upper[["a"]], 5)
})

test_that("split_fixed_bounds with no bounds returns empty fixed and NULL bounds", {
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  sp <- split_fixed_bounds(lower = numeric(0), upper = numeric(0), start = start)
  expect_length(sp$fixed, 0)
  expect_null(sp$lower)
  expect_null(sp$upper)
  expect_equal(sp$start, start)
})
