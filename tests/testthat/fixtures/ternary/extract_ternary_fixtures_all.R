# One-off generator: extract ALL ternary raw data into CSV fixtures from Sam's
# "CA Juveniles A4 fitting DR ..._simplified.xls" workbooks (received 2026-06-02,
# under recieved/20260602/). Requires `readxl`. Run from the repo root:
#   R -q -e 'source("tests/testthat/fixtures/ternary/extract_ternary_fixtures_all.R")'
#
# These are the THREE ternary mixture experiments done with the 4 chemicals
# (CPF / MPs / IMI / FBSA): FBSA+CPF+IMI, MPs+CPF+IMI, MPs+FBSA+CPF. (The 4th
# possible triple, MPs+FBSA+IMI, was not run.) Each workbook's primary sheet
# "CA <chems>" holds the full design: shared controls (0,0,0) + every single-
# chemical arm + all binary/ternary mixture rays, pooled in one table -- the
# ratio-specific tabs ("...33% FBSA", "CA 20CPF 60FBSA 20IMI", etc.) are just
# analysis views of subsets, so the primary sheet is the complete raw dataset.
#
# Layout (verified, all three): row 31 = "Chem 1|Chem 2|Chem 3", row 32 = chem
# names + "RGR" (the reproduction response, col D) + "CA" (a model col, F).
# Data from row 33; a blank separator row sits inside the block. Columns A/B/C
# are the three chemical doses in the order printed in row 32 (which differs per
# workbook, see `order` below) -> C1/C2/C3 positionally; D=RGR -> Res. We read
# whole columns A:D and keep only fully-numeric rows (drops header/param/blank
# rows). Some control/single cells hold non-integer "refitted to singles" values
# (e.g. 860.942476056172) -- kept verbatim, that IS the data, consistent with
# the binary workbook (see binary-workbook-structure memory).
#
# Output filenames follow the existing convention (ternary_fbsa_cpf_imi_...):
# the chemical order in the SOURCE FILE NAME, not the column order.
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required to regenerate these fixtures.")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
# Each experiment's source .xls and its CSV are co-located in a per-experiment
# subfolder (named to match the binary/ experiment folders). The source .xls
# files are gitignored (local only) -- the committed artifact is the CSV.
base <- "tests/testthat/fixtures/ternary"

numeric_rows <- function(df, cols) {
  vals <- lapply(cols, function(j) suppressWarnings(as.numeric(df[[j]])))
  keep <- Reduce(`&`, lapply(vals, function(v) !is.na(v)))
  out <- as.data.frame(vals, col.names = c("C1", "C2", "C3", "Res"))
  out[keep, , drop = FALSE]
}

extract <- function(subdir, file, sheet, out, order) {
  raw <- readxl::read_excel(file.path(base, subdir, file), sheet = sheet,
                            range = readxl::cell_cols("A:D"),
                            col_names = c("C1", "C2", "C3", "Res"))
  dat <- numeric_rows(as.data.frame(raw), 1:4)
  write.csv(dat, file.path(base, subdir, out), row.names = FALSE)
  cat(sprintf("wrote %3d rows  %-16s/%-38s <- %-16s [A,B,C = %s]\n",
              nrow(dat), subdir, out, sheet, order))
}

extract("fbsa_cpf_imi", "CA Juveniles A4 fitting DR FBSA CPF IMI_simplified.xls",
        "CA CPF FBSA IMI", "ternary_fbsa_cpf_imi_continuous.csv", "CPF,FBSA,IMI")
extract("cpf_mps_imi",  "CA Juveniles A4 fitting DR MPs CPF IMI_simplified.xls",
        "CA MPs CPF IMI",  "ternary_mps_cpf_imi_continuous.csv",  "CPF,MPs,IMI")
extract("mps_cpf_fbsa", "CA Juveniles A4 fitting DR MPs FBSA CPF_simplified.xls",
        "CA MPs FBSA CPF", "ternary_mps_fbsa_cpf_continuous.csv",  "MPs,FBSA,CPF")
