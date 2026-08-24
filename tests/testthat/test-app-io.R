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
  expect_false("max" %in% names(sp$lower))   # pinned param dropped from bounds
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

test_that("upload_schema knows the ternary stage", {
  expect_equal(upload_schema("ternary", "continuous"), c("C1", "C2", "C3", "Res"))
  expect_equal(upload_schema("ternary", "quantal"),
               c("C1", "C2", "C3", "Affected", "Exposed"))
})

test_that("ternary template has the schema columns and spans every tier", {
  d <- template_df("ternary", "continuous")
  expect_equal(names(d), c("C1", "C2", "C3", "Res"))
  expect_true(all(vapply(d, is.numeric, logical(1))))
  cls <- as.character(classify_rows(d))
  expect_true(all(c("control", "single", "binary", "ternary") %in% cls))
})

test_that("validate_upload accepts a valid ternary file", {
  good <- data.frame(C1 = c(0, 1, 0, 0, 1), C2 = c(0, 0, 1, 0, 1),
                     C3 = c(0, 0, 0, 1, 1), Res = c(100, 60, 60, 60, 20))
  expect_length(validate_upload(good, "ternary", "continuous"), 0)
})

test_that("validate_upload flags a ternary file with no ternary rows", {
  noternary <- data.frame(C1 = c(0, 1, 0), C2 = c(0, 0, 1),
                          C3 = c(0, 0, 0), Res = c(100, 60, 60))
  expect_match(paste(validate_upload(noternary, "ternary", "continuous"), collapse = " "),
               "ternary rows")
})

test_that("validate_upload catches a negative C3", {
  bad <- data.frame(C1 = c(0, 1), C2 = c(0, 1), C3 = c(-1, 1), Res = c(100, 20))
  expect_match(paste(validate_upload(bad, "ternary", "continuous"), collapse = " "),
               ">= 0|negative|0")
})

test_that("assemble_curve_params3 averages max and uses underscore ec50 names", {
  f1 <- list(par = c(max = 870, slope = 4.7, ec50 = 0.13))
  f2 <- list(par = c(max = 872, slope = 13,  ec50 = 5.58))
  f3 <- list(par = c(max = 874, slope = 3.6, ec50 = 0.57))
  p <- assemble_curve_params3(f1, f2, f3)
  expect_equal(names(p), c("max", "slope1", "slope2", "slope3",
                           "ec50_1", "ec50_2", "ec50_3"))
  expect_equal(unname(p[["max"]]), 872)         # mean(870, 872, 874)
  expect_equal(unname(p[["slope2"]]), 13)
  expect_equal(unname(p[["ec50_3"]]), 0.57)
})

test_that("campaign_chems reports the stressors that are actually dosed", {
  df <- campaign_fixture()
  expect_equal(campaign_chems(df), c(1L, 2L, 3L))
  expect_equal(campaign_n_chem(df), 3L)

  df$C3 <- 0                              # column present but never dosed
  expect_equal(campaign_chems(df), c(1L, 2L))
  expect_equal(campaign_n_chem(df), 2L)

  two <- df[c("C1", "C2", "Res")]         # column absent entirely
  expect_equal(campaign_n_chem(two), 2L)
})

test_that("single_df keeps the control row and renames the concentration to C1", {
  df <- campaign_fixture()

  s2 <- single_df(df, 2)
  expect_equal(names(s2), c("C1", "Res"))
  expect_equal(s2$C1, c(0, 1, 2, 3))      # shared control + chem-2 series
  expect_equal(s2$Res, c(100, 92, 84, 76))

  s1 <- single_df(df, 1)                  # chem 1 needs no rename
  expect_equal(s1$C1, c(0, 1, 2, 3))

  two <- df[df$C3 == 0, c("C1", "C2", "Res")]
  expect_equal(single_df(two, 2)$C1, c(0, 1, 2, 3))   # works with only C1/C2
})

test_that("pair_df keeps that pair's rows and renames to C1/C2", {
  df <- campaign_fixture()

  p23 <- pair_df(df, 2, 3)
  expect_equal(names(p23), c("C1", "C2", "Res"))
  expect_equal(p23$C1, c(0, 1, 2, 3, 0, 0, 0, 1, 2))  # C2 became C1
  expect_equal(p23$C2, c(0, 0, 0, 0, 1, 2, 3, 1, 2))  # C3 became C2
  expect_equal(nrow(p23), 9L)                         # exactly the C1 == 0 rows

  p12 <- pair_df(df, 1, 2)
  expect_equal(names(p12), c("C1", "C2", "Res"))
  expect_equal(nrow(p12), 9)              # control + 3 C1 singles + 3 C2 singles + 2 mixture rows
  expect_true(all(p12$C1 >= 0))
})

test_that("pair_df on a two-stressor frame is a no-op slice", {
  two <- campaign_fixture()[c("C1", "C2", "Res")]
  two <- two[two$C1 > 0 | two$C2 > 0 | seq_len(nrow(two)) == 1, ]
  out <- pair_df(two, 1, 2)
  expect_equal(names(out), c("C1", "C2", "Res"))
  expect_equal(nrow(out), nrow(two))
})

test_that("pair_df rejects a descending pair rather than corrupting the frame", {
  expect_error(pair_df(campaign_fixture(), 2, 1))
})

test_that("pair_base maps a three-stressor campaign base onto one pair's binary names", {
  base3 <- c(max = 872.2065, slope1 = 4.6740, slope2 = 11.2592, slope3 = 3.6291,
            ec50_1 = 0.12747, ec50_2 = 34.6707, ec50_3 = 0.57505)

  p13 <- pair_base(base3, 1, 3)
  expect_equal(names(p13), c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_equal(unname(p13[["max"]]),    872.2065)
  expect_equal(unname(p13[["slope1"]]), 4.6740)
  # stressor 3's slope/EC50 land in slot 2 -- NOT stressor 2's (11.2592 /
  # 34.6707), which is the exact mis-mapping the bug produced. Stressor 2's
  # EC50 (~34.7) and stressor 3's (~0.575) are nearly two orders of magnitude
  # apart, so a wrong mapping is unmissable.
  expect_equal(unname(p13[["slope2"]]), 3.6291)
  expect_equal(unname(p13[["ec501"]]),  0.12747)
  expect_equal(unname(p13[["ec502"]]),  0.57505)

  p23 <- pair_base(base3, 2, 3)
  expect_equal(unname(p23[["slope1"]]), 11.2592)
  expect_equal(unname(p23[["slope2"]]), 3.6291)
  expect_equal(unname(p23[["ec501"]]),  34.6707)
  expect_equal(unname(p23[["ec502"]]),  0.57505)
})

test_that("pair_base passes a two-stressor campaign base through unchanged", {
  base2 <- c(max = 800, slope1 = 2, slope2 = 3, ec501 = 0.1, ec502 = 0.5)
  expect_identical(pair_base(base2, 1, 2), base2)
})

test_that("campaign schema requires C1/C2 and treats C3 as optional", {
  expect_equal(upload_schema("campaign", "continuous"), c("C1", "C2", "Res"))
  expect_equal(upload_schema("campaign", "quantal"),
               c("C1", "C2", "Affected", "Exposed"))
})

test_that("the campaign template is a full three-stressor campaign", {
  tpl <- template_df("campaign", "continuous")
  expect_true(all(c("C1", "C2", "C3", "Res") %in% names(tpl)))
  cls <- classify_rows(tpl)
  expect_true(all(c("control", "single", "binary", "ternary") %in%
                  as.character(cls)))
})

test_that("a well-formed campaign validates clean", {
  expect_equal(validate_upload(campaign_fixture(), "campaign", "continuous"),
               character(0))
})

test_that("a campaign needs at least two dosed stressors", {
  df <- campaign_fixture()
  df$C2 <- 0
  df$C3 <- 0
  errs <- validate_upload(df, "campaign", "continuous")
  expect_true(any(grepl("at least two stressors", errs)))
})

test_that("each single series needs 4 distinct concentrations", {
  df <- campaign_fixture()
  df <- df[!(df$C2 > 0 & df$C1 == 0 & df$C3 == 0 & df$C2 > 2), ]  # thin chem 2
  errs <- validate_upload(df, "campaign", "continuous")
  expect_true(any(grepl("Stressor 2", errs)))
  expect_true(any(grepl("distinct concentrations", errs)))
})

test_that("a pair with no mixture rows is NOT an error (the sub-tab disables instead)", {
  df <- campaign_fixture()
  df <- df[!(df$C2 > 0 & df$C3 > 0), ]     # drop every 2x3 mixture row
  expect_equal(validate_upload(df, "campaign", "continuous"), character(0))
})
