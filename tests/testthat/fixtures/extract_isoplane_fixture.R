# One-off: extract the precomputed FBSA/CPF/IMI isoplane reference data from
# Skylar's Isoplains_Mixtox.xlsx into CSV fixtures, for validating the
# model-driven ec50_isoplane / sigma_tu_curve / ec50_markers functions.
# Run from the repo root:
#   R -q -e 'source("tests/testthat/fixtures/extract_isoplane_fixture.R")'
# Sheet layouts verified 2026-06-01 (see plan Task 1).
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required to regenerate this fixture.")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
wb <- "recieved/R scripts & excel/R scripts & excel/Isoplains_Mixtox.xlsx"

# Isoplane point cloud: one row per (point, model). Concentrations map
# CPF->C1, FBSA->C2, IMI->C3 (matches the C1/C2/C3 convention of the fixture
# used by analyse_ternary, where CPF=C1, FBSA=C2, IMI=C3).
pts <- as.data.frame(readxl::read_excel(wb, sheet = "FBSA_CPF_IMI"))
pts <- pts[, c("CPF", "FBSA", "IMI", "Mod")]
names(pts) <- c("C1", "C2", "C3", "model")
pts <- pts[stats::complete.cases(pts), , drop = FALSE]
write.csv(pts, "tests/testthat/fixtures/isoplane_points.csv", row.names = FALSE)
cat("wrote", nrow(pts), "isoplane points (",
    paste(unique(pts$model), collapse = ", "), ")\n")

# z-value curves for this mixture.
zv <- as.data.frame(readxl::read_excel(wb, sheet = "TU-zValues"))
zv <- zv[zv$Exp == "Mixtox1", c("z", "TU", "Chem", "Mod"), drop = FALSE]
zv <- zv[stats::complete.cases(zv), , drop = FALSE]
write.csv(zv, "tests/testthat/fixtures/isoplane_zvalues.csv", row.names = FALSE)
cat("wrote", nrow(zv), "z-value rows\n")

# EC50 marker points for this mixture.
ec <- as.data.frame(readxl::read_excel(wb, sheet = "TER_EC50"))
ec <- ec[ec$Exp == "Mixtox1", c("z", "TU", "Chem", "Ratio"), drop = FALSE]
ec <- ec[stats::complete.cases(ec), , drop = FALSE]
write.csv(ec, "tests/testthat/fixtures/isoplane_ec50.csv", row.names = FALSE)
cat("wrote", nrow(ec), "EC50 marker rows\n")
