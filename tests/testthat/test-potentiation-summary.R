ex <- function() read.csv(test_path("fixtures", "potentiation",
                                    "potentiation_example_continuous.csv"))

# Minimal hand-built analysis result: one row per level, in file order.
fake_res <- function(sig, dir, p = rep(0.01, length(sig)),
                     extrapolated = rep(FALSE, length(sig)),
                     ref_extrapolated = FALSE) {
  k <- length(sig)
  lv <- paste("Level", seq_len(k))
  list(
    series = c("Stressor 1", lv), reference = "Stressor 1",
    curves = list(`Stressor 1` = list(ec50_extrapolated = ref_extrapolated)),
    comparisons = data.frame(
      series = lv, p = p, ec50_reference = 10,
      ec50_series = ifelse(dir == "more toxic", 5, 20),
      significant = sig, direction = ifelse(sig, dir, "no significant influence"),
      ec50_extrapolated = extrapolated, stringsAsFactors = FALSE))
}

test_that("workbook example reads as the spec's example", {
  s <- analyse_potentiation(ex())$summary
  expect_match(s, "less toxic with Level 1 (EC50", fixed = TRUE)
  expect_match(s, "Level 3 (EC50", fixed = TRUE)
  expect_match(s, "extrapolated", fixed = TRUE)
  expect_match(s, "no significant influence with Level 2", fixed = TRUE)
  expect_match(s, "The effect does not increase consistently with the level.$")
  expect_false(grepl(" from ", s, fixed = TRUE))
})

test_that("onset named when every later level is significant", {
  s <- potentiation_summary(fake_res(c(FALSE, TRUE, TRUE),
                                     rep("less toxic", 3)))
  expect_match(s, "from Level 2 on", fixed = TRUE)
  expect_match(s, "Level 2 (EC50 10 → 20, p = 0.01)", fixed = TRUE)
})

test_that("mixed directions list both and give no onset", {
  s <- potentiation_summary(fake_res(c(TRUE, TRUE),
                                     c("more toxic", "less toxic")))
  expect_match(s, "more toxic with Level 1", fixed = TRUE)
  expect_match(s, "less toxic with Level 2", fixed = TRUE)
  expect_false(grepl(" from ", s, fixed = TRUE))
})

test_that("single comparison", {
  expect_equal(
    potentiation_summary(fake_res(FALSE, "more toxic", p = 0.42)),
    "No significant influence of Level 1 on the toxicity of Stressor 1 (p = 0.42).")
  expect_equal(
    potentiation_summary(fake_res(TRUE, "more toxic", p = 0.003,
                                  extrapolated = TRUE)),
    paste0("Stressor 1 is significantly more toxic with Level 1 ",
           "(EC50 10 → 5, p = 0.003, extrapolated)."))
})

test_that("none significant", {
  expect_equal(
    potentiation_summary(fake_res(c(FALSE, FALSE), rep("less toxic", 2))),
    "No significant influence on the toxicity of Stressor 1 at any level.")
})

test_that("stressor name is used", {
  s <- potentiation_summary(fake_res(c(TRUE, TRUE), rep("less toxic", 2)),
                            stressor1 = "Copper")
  expect_match(s, "^Copper ")
})

test_that("an extrapolated reference EC50 is called out", {
  s <- potentiation_summary(fake_res(c(TRUE, TRUE), rep("less toxic", 2),
                                     ref_extrapolated = TRUE))
  expect_match(s, paste0("The EC50 of Stressor 1 lies outside its tested ",
                         "concentrations, so every shift is extrapolated.$"))
  expect_false(grepl("outside its tested",
                     potentiation_summary(fake_res(c(TRUE, TRUE),
                                                   rep("less toxic", 2)))))
})

test_that("format_p rounds and floors", {
  expect_equal(format_p(0.000798), "p = 0.0008")
  expect_equal(format_p(0.00743), "p = 0.0074")
  expect_equal(format_p(5e-6), "p < 0.0001")
})
