# Validation of the (de)potentiation analysis against the draft workbook the
# example data was extracted from. Point values live in inst/validation/oracles.csv;
# what remains asserted here are the claims a value-plus-tolerance cannot express.

test_that("(de)potentiation: three-level comparison reproduces the workbook", {
  expect_matches_reference("potentiation_example_continuous")

  res <- validation_fit("potentiation_example_continuous")
  # Level 2 is the no-difference case: the workbook's X-squared is ~0 (-1.3e-5,
  # a numerical artefact of its optimiser); the engine floors at 0.
  chi <- res$comparisons$chi[res$comparisons$series == "Level 2"]
  expect_gte(chi, 0)
  expect_lt(chi, 1e-3)
})
