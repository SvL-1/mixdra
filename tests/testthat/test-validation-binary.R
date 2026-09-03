# Validation of the STAGED fitting method on the binary MPs+CPF workbook data.
#
# The per-quantity targets now live in inst/validation/oracles.csv, which is also
# what the validation report renders -- so an asserted number and a reported
# number cannot drift apart. Rows flagged `divergent` there record where this
# engine intentionally departs from the workbook, with the reason: the workbook
# fits the curve parameters (max, slopes, EC50s) jointly on the mixture data,
# while this engine fixes them from the single-compound data only and then fits
# the interaction parameters. Constraining the reference curves to the marginals
# raises the reference residual (the workbook instead absorbs interaction into
# re-fitted curves) and lets genuine deviations surface -- e.g. the quantal data
# below selects a dose-ratio (DR) deviation that the joint workbook reported as
# 'reference'.
#
# What remains asserted here directly are the claims the oracle table cannot
# express as a value-plus-tolerance: significance inequalities.

test_that("staged engine: binary CA continuous analysis on MPs+CPF", {
  expect_matches_workbook("binary_mps_cpf_ca_continuous")

  res <- validation_fit("binary_mps_cpf_ca_continuous")
  # S/A remains significant, though weaker than the joint workbook's chi ~13.
  ca_sa <- res$comparison[res$comparison$model == "SA", ]
  expect_gt(ca_sa$chi, 3)
  expect_lt(ca_sa$p, 0.05)
})

test_that("staged engine: binary CA quantal analysis on MPs+CPF", {
  expect_matches_workbook("binary_mps_cpf_quantal")

  res <- validation_fit("binary_mps_cpf_quantal")
  # The interaction is strongly supported (CA-vs-S/A p well below 0.001).
  ca_sa <- res$comparison[res$comparison$model == "SA", ]
  expect_lt(ca_sa$p, 0.001)
})
