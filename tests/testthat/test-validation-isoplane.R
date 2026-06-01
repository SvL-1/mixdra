# Validation gate: the model-driven isoplane / Sigma-TU computation must land on
# the fitted model's TRUE EC50 surface. We assert SELF-CONSISTENCY against the
# actual ca_asa_tri predictor: every isoplane / marker / curve point, fed back
# through the model, must return max/2 (50% effect). This is exact for the CA
# Advanced-S/A model at ANY slopes -- at Y = max/2 the ((max-Y)/Y)^(1/slope)
# factor is 1 for every slope, so sigma_tu = F4(z) holds exactly.
#
# We deliberately do NOT validate against Isoplains_Mixtox.xlsx. That Excel tool
# (Model4Isoplane.xlsm) constructs its isoplane by a different, interior-
# approximate method: it agrees with the exact model only near the simplex edges
# (~0.5%) and diverges in the interior (~9%). Our closed form is exact for the
# fitted model, so the Excel is not a comparable reference. (Decision 2026-06-01.)
#
# SLOW: analyse_ternary on 419 rows ~ 2 min. set.seed for stability.

test_that("isoplane / markers / sigma-TU land on the fitted model's max/2 surface", {
  fx <- testthat::test_path("fixtures", "ternary_fbsa_cpf_imi_continuous.csv")
  skip_if_not(file.exists(fx), "ternary fixture CSV missing")

  df <- read.csv(fx)
  set.seed(42)
  res <- analyse_ternary(df, "CA", "continuous", n_starts = 10, time_limit = 120,
                         lower = c(ec50_2 = 5.579), upper = c(ec50_2 = 5.581))
  b <- res$base
  half <- b[["max"]] / 2
  ec <- c(b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]])
  pred <- function(C1, C2, C3, A4)
    ca_asa_tri(C1, C2, C3, b[["max"]], b[["slope1"]], b[["slope2"]], b[["slope3"]],
               b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]],
               res$pairwise[["A1"]], res$pairwise[["A2"]], res$pairwise[["A3"]], A4)

  ## SA isoplane (A4 = 0): every grid point predicts max/2
  sa <- ec50_isoplane(res, "SA", n = 12)
  ysa <- vapply(seq_len(nrow(sa)), function(i)
    pred(sa$C1[i], sa$C2[i], sa$C3[i], 0), numeric(1))
  expect_equal(ysa, rep(half, nrow(sa)), tolerance = 1e-4)

  ## ASA isoplane (A4 = A4_overall): every grid point predicts max/2
  asa <- ec50_isoplane(res, "ASA", n = 12)
  yasa <- vapply(seq_len(nrow(asa)), function(i)
    pred(asa$C1[i], asa$C2[i], asa$C3[i], res$A4_overall), numeric(1))
  expect_equal(yasa, rep(half, nrow(asa)), tolerance = 1e-4)

  ## EC50 markers: each ratio's marker (its individual A4) predicts max/2
  mk <- ec50_markers(res, df)
  a4 <- res$individual$A4[match(mk$ratio, res$individual$ratio)]
  ymk <- vapply(seq_len(nrow(mk)), function(i)
    pred(mk$C1[i], mk$C2[i], mk$C3[i], a4[i]), numeric(1))
  expect_equal(ymk, rep(half, nrow(mk)), tolerance = 1e-4)

  ## sigma_tu_curve: reconstruct concentrations from z and sigma_tu
  ## (focal = z, others = (1-z)/2; C_i = ec_i * zc_i * sigma_tu) -> predicts max/2
  chem_idx <- c(C1 = 1, C2 = 2, C3 = 3)
  cur <- sigma_tu_curve(res, "SA", n = 11)
  ycur <- vapply(seq_len(nrow(cur)), function(i) {
    zc <- rep((1 - cur$z[i]) / 2, 3)
    zc[chem_idx[[cur$chem[i]]]] <- cur$z[i]
    C <- ec * zc * cur$sigma_tu[i]
    pred(C[1], C[2], C[3], 0)
  }, numeric(1))
  expect_equal(ycur, rep(half, nrow(cur)), tolerance = 1e-4)
})
