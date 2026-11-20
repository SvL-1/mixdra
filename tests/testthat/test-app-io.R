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

potentiation_shaped <- function(levels2 = c(0, 1)) {
  c1 <- c(0, 0.1, 0.2, 0.4, 0.8, 1.6)
  data.frame(C1  = rep(c1, length(levels2)),
             C2  = rep(levels2, each = length(c1)),
             Res = 100 - seq_len(length(c1) * length(levels2)))
}

test_that("sparse_level_stressors finds a stressor applied at only a few levels", {
  expect_equal(sparse_level_stressors(potentiation_shaped(c(0, 1)), 1:2), c("2" = 1L))
  expect_equal(sparse_level_stressors(potentiation_shaped(c(0, 0.5, 1, 2)), 1:2),
               c("2" = 3L))
  # four levels is a dose-response curve, not a handful of fixed levels
  expect_length(sparse_level_stressors(potentiation_shaped(c(0, 0.5, 1, 2, 4)), 1:2), 0)
})

test_that("sparse_level_stressors stays quiet when no stressor has a full curve", {
  # every stressor at three levels: nothing here is a reference curve to compare
  # against, so this is not the (de)potentiation case
  expect_length(sparse_level_stressors(campaign_fixture(), 1:3), 0)
})

test_that("the bundled example campaign is not flagged as (de)potentiation", {
  df <- to_engine_df(read_upload(system.file(
    "extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")),
    "campaign")
  expect_null(potentiation_design_note(df, c("FBSA", "CPF", "IMI"),
                                       campaign_chems(df)))
})

test_that("potentiation_design_note names the stressor and points to the tab", {
  note <- potentiation_design_note(potentiation_shaped(c(0, 0.5, 1)),
                                   c("CPF (mg/kg)", "Heat (C)"), 1:2)
  txt <- paste(as.character(note), collapse = "")
  expect_match(txt, "Heat (C)", fixed = TRUE)
  expect_match(txt, "only 2 concentrations", fixed = TRUE)
  expect_match(txt, "(De)potentiation", fixed = TRUE)
  expect_false(grepl("CPF", txt, fixed = TRUE))
  one <- potentiation_design_note(potentiation_shaped(c(0, 1)), c("A", "B"), 1:2)
  expect_match(paste(as.character(one), collapse = ""), "only 1 concentration.",
               fixed = TRUE)
})

test_that("a rejected (de)potentiation-shaped campaign upload says where it belongs", {
  errs <- validate_upload(potentiation_shaped(c(0, 1)), "campaign", "continuous")
  expect_true(any(grepl("Stressor 2 has only 2 distinct", errs, fixed = TRUE)))
  hint <- errs[grepl("(De)potentiation tab", errs, fixed = TRUE)]
  expect_length(hint, 1)
  expect_match(hint, "Stressor 2 at only 1 concentration.", fixed = TRUE)
  expect_false(grepl("with caution", hint, fixed = TRUE))
  # a mixture design that is merely short of single-stressor points gets no hint
  thin_errs <- validate_upload(campaign_fixture()[-4, ], "campaign", "continuous")
  expect_true(any(grepl("Stressor 1 has only 3 distinct", thin_errs, fixed = TRUE)))
  expect_false(any(grepl("(De)potentiation", thin_errs, fixed = TRUE)))
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

test_that("panel_constraints turns a ticked Fix box into equal lower/upper", {
  vals <- list(lo_slope = 1, hi_slope = 5,
               pin_ec50 = TRUE, val_ec50 = 0.42,
               pin_max = FALSE, val_max = 900)
  b <- panel_constraints(vals)
  expect_equal(b$lower[["ec50"]], 0.42)
  expect_equal(b$upper[["ec50"]], 0.42)
  expect_equal(b$lower[["slope"]], 1)
  expect_equal(b$upper[["slope"]], 5)
  expect_false("max" %in% names(b$lower))   # unticked Fix contributes nothing
})

test_that("panel_constraints ignores a ticked Fix box with no value", {
  b <- panel_constraints(list(pin_ec50 = TRUE, val_ec50 = NA))
  expect_null(b$lower)
  expect_null(b$upper)
})

test_that("a pinned panel parameter survives split_fixed_bounds as fixed", {
  b  <- panel_constraints(list(pin_ec50 = TRUE, val_ec50 = 0.42))
  sp <- split_fixed_bounds(b$lower, b$upper, c(max = 1, slope = 1, ec50 = 9))
  expect_equal(sp$fixed, "ec50")
  expect_equal(sp$start[["ec50"]], 0.42)
  expect_null(sp$lower)
})

test_that("base_param_names follows the registry's per-arity EC50 spelling", {
  expect_equal(unname(base_param_names(2, 2)), c("slope2", "ec502"))
  expect_equal(unname(base_param_names(2, 3)), c("slope2", "ec50_2"))
})

test_that("campaign_base_bounds renames per-stressor bounds to fit positions", {
  panels <- list("1" = list(lower = c(ec50 = 0.1), upper = c(ec50 = 0.2)),
                 "3" = list(lower = c(slope = 2),  upper = NULL))
  b <- campaign_base_bounds(panels, max_bounds = NULL, chems = c(1L, 2L, 3L))
  expect_equal(b$lower[["ec50_1"]], 0.1)
  expect_equal(b$upper[["ec50_1"]], 0.2)
  expect_equal(b$lower[["slope3"]], 2)
})

test_that("campaign_base_bounds maps a non-contiguous stressor set by position", {
  # Stressors 1 and 3 only: stressor 3 is the SECOND column of the fit.
  panels <- list("3" = list(lower = c(ec50 = 0.5), upper = NULL))
  b <- campaign_base_bounds(panels, max_bounds = NULL, chems = c(1L, 3L))
  expect_equal(b$lower[["ec502"]], 0.5)   # binary names for a 2-stressor fit
})

test_that("campaign_base_bounds takes max from the campaign row, not the panels", {
  panels <- list("1" = list(lower = c(max = 100), upper = c(max = 200)))
  b <- campaign_base_bounds(panels, max_bounds = list(lower = c(max = 850),
                                                      upper = c(max = 850)),
                            chems = c(1L, 2L))
  expect_equal(b$lower[["max"]], 850)
  expect_equal(b$upper[["max"]], 850)
})

test_that("pair_bounds renames a ternary campaign's bounds onto the pair", {
  b <- list(lower = c(max = 850, slope2 = 1, ec50_2 = 0.3, ec50_3 = 0.9),
            upper = c(ec50_3 = 1.1))
  p <- pair_bounds(b, 2, 3)
  expect_equal(p$lower[["slope1"]], 1)      # stressor 2 is the pair's first
  expect_equal(p$lower[["ec501"]], 0.3)
  expect_equal(p$lower[["ec502"]], 0.9)
  expect_equal(p$lower[["max"]], 850)
  expect_equal(p$upper[["ec502"]], 1.1)
})

test_that("pair_bounds drops bounds for a stressor outside the pair", {
  b <- list(lower = c(ec50_1 = 0.1, ec50_2 = 0.3), upper = NULL)
  p <- pair_bounds(b, 2, 3)
  expect_equal(names(p$lower), "ec501")
})

# ---- "potentiation" stage ----------------------------------------------------

pot_df <- function(series, conc = c(0, 0.1, 0.3, 1, 3, 10)) {
  do.call(rbind, lapply(series, function(s)
    data.frame(Series = s, C1 = conc, Res = round(100 / (1 + conc), 1))))
}

test_that("potentiation schema has Series first", {
  expect_equal(upload_schema("potentiation", "continuous"), c("Series", "C1", "Res"))
  expect_equal(upload_schema("potentiation", "quantal"),
               c("Series", "C1", "Affected", "Exposed"))
})

test_that("potentiation template matches schema and validates", {
  for (r in c("continuous", "quantal")) {
    t <- template_df("potentiation", r)
    expect_equal(names(t), upload_schema("potentiation", r))
    expect_equal(unique(t$Series), c("Stressor 1", "Level 1", "Level 2"))
    expect_equal(validate_upload(t, "potentiation", r), character(0))
  }
})

test_that("potentiation: Series is text, C1/Res are still numeric-checked", {
  d <- pot_df(c("Stressor 1", "Level 1"))
  expect_equal(validate_upload(d, "potentiation", "continuous"), character(0))
  d$C1 <- as.character(d$C1)
  d$Res <- as.character(d$Res)
  e <- validate_upload(d, "potentiation", "continuous")
  expect_match(e, "Non-numeric column(s): C1, Res", all = FALSE, fixed = TRUE)
  expect_false(any(grepl("Series", e[grepl("Non-numeric", e)])))
})

test_that("potentiation: missing Series, empty labels, one series", {
  d <- pot_df(c("Stressor 1", "Level 1"))
  expect_match(validate_upload(d[c("C1", "Res")], "potentiation", "continuous"),
               "Missing column(s): Series", all = FALSE, fixed = TRUE)
  d2 <- d; d2$Series[3] <- ""
  expect_match(validate_upload(d2, "potentiation", "continuous"),
               "Series labels must not be empty.", all = FALSE, fixed = TRUE)
  d3 <- d; d3$Series[3] <- NA
  expect_match(validate_upload(d3, "potentiation", "continuous"),
               "Series labels must not be empty.", all = FALSE, fixed = TRUE)
  expect_match(validate_upload(pot_df("Stressor 1"), "potentiation", "continuous"),
               "at least two series", all = FALSE)
})

test_that("potentiation: a series needs 3 distinct positive concentrations", {
  d <- rbind(pot_df("Stressor 1"), pot_df("Level 1", conc = c(0, 0.1, 1, 1)))
  expect_match(validate_upload(d, "potentiation", "continuous"),
               "Series 'Level 1' has only 2 positive concentrations",
               all = FALSE, fixed = TRUE)
})

test_that("potentiation: numeric Series labels validate and analyse", {
  d <- pot_df(c(1, 2))
  expect_equal(validate_upload(d, "potentiation", "continuous"), character(0))
  r <- analyse_potentiation(d)
  expect_equal(r$series, c("1", "2"))
})

test_that("potentiation: labels differing by whitespace are one series", {
  # untrimmed this would be 3 series; trimmed it is 2, so no error
  d <- rbind(pot_df("Stressor 1"), pot_df("Level 1 "), pot_df("Level 1"))
  expect_equal(validate_upload(d, "potentiation", "continuous"), character(0))
  # trimmed, "Stressor 1" and " Stressor 1" are ONE series, so only one here
  d1 <- rbind(pot_df("Stressor 1"), pot_df(" Stressor 1"))
  expect_match(validate_upload(d1, "potentiation", "continuous"),
               "at least two series", all = FALSE)
})

test_that("potentiation: rows with missing values are reported by row number", {
  d <- pot_df(c("Stressor 1", "Level 1", "Level 2"))   # 18 rows
  d$Res[12] <- NA; d$C1[4] <- NA
  expect_match(validate_upload(d, "potentiation", "continuous"),
               "Row(s) 4, 12 have missing values; remove or complete them.",
               all = FALSE, fixed = TRUE)
  t <- template_df("potentiation", "quantal")
  t$Affected[2] <- NA; t$Exposed[7] <- NA
  expect_match(validate_upload(t, "potentiation", "quantal"),
               "Row(s) 2, 7 have missing values", all = FALSE, fixed = TRUE)
  d3 <- pot_df(c("Stressor 1", "Level 1")); d3$Series[3] <- NA
  expect_match(validate_upload(d3, "potentiation", "continuous"),
               "Row(s) 3 have missing values", all = FALSE, fixed = TRUE)
  # a long list is cut short
  d4 <- pot_df(c("Stressor 1", "Level 1", "Level 2")); d4$Res <- NA_real_
  expect_match(validate_upload(d4, "potentiation", "continuous"),
               "Row(s) 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 and 8 more have missing values",
               all = FALSE, fixed = TRUE)
  expect_false(any(grepl("missing values",
    validate_upload(pot_df(c("Stressor 1", "Level 1")), "potentiation", "continuous"))))
})

test_that("potentiation quantal template declines with dose, as the engine fits it", {
  t <- template_df("potentiation", "quantal")
  for (s in unique(t$Series)) {
    p <- with(t[t$Series == s, ], Affected / Exposed)
    expect_equal(p[[1]], 1, label = s)
    expect_true(all(diff(p) <= 0), label = s)
  }
  r <- analyse_potentiation(t, response = "quantal")
  expect_equal(r$curves[["Stressor 1"]]$par[["ec50"]], 3, tolerance = 0.3)
  expect_true(all(r$comparisons$ratio < 1))
})

test_that("potentiation quantal: Exposed must be > 0", {
  t <- template_df("potentiation", "quantal")
  t$Exposed[2] <- 0; t$Affected[2] <- 0
  expect_match(validate_upload(t, "potentiation", "quantal"),
               "Exposed must be greater than 0 in every row.", all = FALSE, fixed = TRUE)
})
