# interaction_help() is a pure helper: (reference, deviation) -> Shiny markup.
# These tests render the tagList to an HTML string and assert the formula
# fragment + a/b wording for each reference x deviation combination.

html_of <- function(reference, deviation) {
  as.character(interaction_help(reference, deviation))
}

test_that("reference deviation states there is no interaction term", {
  h <- html_of("CA", "reference")
  expect_match(h, "No interaction", ignore.case = TRUE)
  expect_match(h, "CA", fixed = TRUE)
})

test_that("reference deviation names the IA baseline when reference is IA", {
  expect_match(html_of("IA", "reference"), "IA", fixed = TRUE)
})

test_that("SA shows the F = a.product formula and the a meaning", {
  h <- html_of("CA", "SA")
  expect_match(h, "F = a&middot;&prod;z", fixed = TRUE)
  expect_match(h, "strength and direction", fixed = TRUE)
})

test_that("DR shows a summation formula and the mixture-ratio meaning", {
  h <- html_of("CA", "DR")
  expect_match(h, "&Sigma;", fixed = TRUE)
  expect_match(h, "mixture ratio", fixed = TRUE)
})

test_that("DL formula differs between CA (SigmaTU) and IA (P)", {
  ca <- html_of("CA", "DL")
  ia <- html_of("IA", "DL")
  expect_match(ca, "&Sigma;TU", fixed = TRUE)
  expect_match(ia, "b&middot;P)", fixed = TRUE)
  expect_match(ia, "IA-predicted", fixed = TRUE)
})

test_that("the sign convention is stated for models with an interaction term", {
  for (dev in c("SA", "DR", "DL")) {
    h <- html_of("CA", dev)
    expect_match(h, "antagonism", fixed = TRUE)
    expect_match(h, "synergism", fixed = TRUE)
  }
})

test_that("the sign convention is omitted for the reference (no a term)", {
  h <- html_of("CA", "reference")
  expect_false(grepl("antagonism", h, fixed = TRUE))
})
