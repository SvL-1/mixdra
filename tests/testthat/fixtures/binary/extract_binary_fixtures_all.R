# One-off generator: extract ALL binary-mixture data from Sam's workbook
#   "binary mixtox calculations CA continuous refitted to singles and
#    binary residuals_correct.xlsm"
# into the CSV fixture format used in tests/testthat/fixtures.
# The source workbook lives (gitignored, local-only) next to this script in
# binary/; needs `readxl`. Run from the repo root:
#   R -q -e 'source("tests/testthat/fixtures/binary/extract_binary_fixtures_all.R")'
#
# WHAT THIS WORKBOOK IS (verified, see memory note + Jelle's question):
#   The 4 chemicals (CPF / MPs / IMI / FBSA) were tested in FOUR ternary
#   mixture experiments ("campaigns"), each combining 3 of the 4 chemicals.
#   Each campaign is re-expressed as its 3 binary-pair analysis sheets named
#   "CA - <Chem1> <Chem2> (<third chem>)", where the parenthetical names the
#   third chemical co-tested in that campaign. Hence 11 sheets, not 6: most
#   pairs recur once per campaign they appeared in.
#   *Within* a campaign the 3 sheets SHARE rows: the identical (0,0) control
#   block (20 reps) and each chemical's single-dose arm are reused verbatim
#   (verified: campaign {CPF,MPs,IMI} controls identical across its 3 sheets;
#   CPF singles identical in its two CPF sheets). *Across* campaigns the same
#   pair has different replicate data (different physical run). So the 11 are
#   NOT independent datasets -- they are correlated within a campaign. The
#   per-pair duplication of controls/singles is intentional: each binary CA/IA
#   fit needs the controls + both single arms + that pair's mixture rows.
#
# Layout (all sheets): header row 24 = "<Chem1> <Chem2> Data ...".
#   Continuous (reproduction): col A=Chem1 dose -> C1, B=Chem2 dose -> C2,
#     C=Data -> Res. (Some FBSA-campaign rows hold non-integer "refitted to
#     singles" control/single values -- kept verbatim, that IS the data.)
#   Quantal (survival): A=Chem1 -> C1, B=Chem2 -> C2, C=Survivors -> Affected,
#     D=Exposed. (Matches existing binary_mps_cpf_quantal.csv: the "Affected"
#     column holds the survivor count; observed pi = Affected / Exposed.)
#   A blank separator row sits inside some blocks, so we read whole columns and
#   keep only rows where every needed cell is numeric (robust, like the
#   original extract_binary_fixture.R).
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to regenerate these fixtures.")
}
wb <- "tests/testthat/fixtures/binary/binary mixtox calculations CA continuous refitted to singles and binary residuals_correct.xlsm"
outdir <- "tests/testthat/fixtures/binary"

numeric_rows <- function(df, cols) {
  vals <- lapply(cols, function(j) suppressWarnings(as.numeric(df[[j]])))
  keep <- Reduce(`&`, lapply(vals, function(v) !is.na(v)))
  out <- as.data.frame(vals, col.names = names(df)[cols])
  out[keep, , drop = FALSE]
}

extract <- function(sheet, subdir, file, col_names) {
  rng <- readxl::cell_cols(if (length(col_names) == 3) "A:C" else "A:D")
  raw <- readxl::read_excel(wb, sheet = sheet, range = rng,
                            col_names = col_names)
  dat <- numeric_rows(as.data.frame(raw), seq_along(col_names))
  dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
  write.csv(dat, file.path(outdir, subdir, file), row.names = FALSE)
  cat(sprintf("wrote %3d rows  %-16s/%-40s <- %s\n", nrow(dat), subdir, file, sheet))
}

# Fixtures are grouped by the 3-chemical EXPERIMENT they come from (the binary
# sheets are pairwise projections of those ternary experiments -- verified: the
# (IMI)-campaign controls are byte-identical to the ternary MPs CPF IMI
# controls). Subdir = experiment chem-triple; survival/ holds the quantal
# (different endpoint).
# --- continuous (reproduction): C1, C2, Res ---------------------------------
cont <- list(
  c("CA - MPs CPF (IMI)",  "cpf_mps_imi",  "binary_ca_mps_cpf_imi_continuous.csv"),
  c("CA - MPs IMI (CPF)",  "cpf_mps_imi",  "binary_ca_mps_imi_cpf_continuous.csv"),
  c("CA -IMI CPF (MPs)",   "cpf_mps_imi",  "binary_ca_imi_cpf_mps_continuous.csv"),
  c("CA - CPF FBSA (IMI)", "fbsa_cpf_imi", "binary_ca_cpf_fbsa_imi_continuous.csv"),
  c("CA - FBSA IMI (CPF)", "fbsa_cpf_imi", "binary_ca_fbsa_imi_cpf_continuous.csv"),
  c("CA - CPF IMI (FBSA)", "fbsa_cpf_imi", "binary_ca_cpf_imi_fbsa_continuous.csv"),
  c("IA - continous data model (new)", "fbsa_cpf_imi", "binary_ia_cpf_imi_continuous.csv"),
  c("CA - MPs FBSA (IMI)", "mps_fbsa_imi", "binary_ca_mps_fbsa_imi_continuous.csv"),
  c("CA - MPs IMI (FBSA)", "mps_fbsa_imi", "binary_ca_mps_imi_fbsa_continuous.csv"),
  c("CA - FBSA IMI (MPs)", "mps_fbsa_imi", "binary_ca_fbsa_imi_mps_continuous.csv"),
  c("CA - MPs CPF (FBSA)", "mps_cpf_fbsa", "binary_ca_mps_cpf_fbsa_continuous.csv"),
  c("CA - CPF FBSA (MPs)", "mps_cpf_fbsa", "binary_ca_cpf_fbsa_mps_continuous.csv")
)
for (m in cont) extract(m[1], m[2], m[3], c("C1", "C2", "Res"))

# --- quantal (survival): C1, C2, Affected (= survivors), Exposed ------------
quant <- list(
  c("CA - binary data model (new)", "survival", "binary_ca_cpf_imi_quantal.csv"),
  c("IA - binary data model (new)", "survival", "binary_ia_cpf_imi_quantal.csv")
)
for (m in quant) extract(m[1], m[2], m[3], c("C1", "C2", "Affected", "Exposed"))
