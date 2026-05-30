# One-off generator: extract the binary MPs+CPF/Imi workbook data into CSV
# fixtures checked into the repo. Requires the workbook at the repo root and the
# `readxl` package (a Suggested dependency). Run from the repo root:
#   R -q -e 'source("tests/testthat/fixtures/extract_binary_fixture.R")'
#
# Layout (verified against the workbook):
#   "CA - continous data model (new)": header row "[CPF] [MPs] Data",
#       cols A=[CPF]->C1, B=[MPs]->C2, C=Data->Res. 145 data rows.
#   "CA - binary data model (new)":    header row "[MPs] [Imi] Survived Exposed"
#       (a blank row follows the header), cols A=[MPs]->C1, B=[Imi]->C2,
#       C=Survived->Affected, D=Exposed. 145 data rows.
#   The endpoint is survival (decreasing): the response count is the number of
#   SURVIVORS, so observed proportion pi = Survived / Exposed.
#
# Rather than hardcode fragile cell ranges, we read whole columns and keep only
# the rows where every needed cell is numeric — robust to header position and
# the blank separator row.
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to regenerate these fixtures.")
}
wb <- "MixTox Model_binary_MPs_CPF.xlsm"

numeric_rows <- function(df, cols) {
  vals <- lapply(cols, function(j) suppressWarnings(as.numeric(df[[j]])))
  keep <- Reduce(`&`, lapply(vals, function(v) !is.na(v)))
  out <- as.data.frame(vals, col.names = names(df)[cols])
  out[keep, , drop = FALSE]
}

cont_raw <- readxl::read_excel(wb, sheet = "CA - continous data model (new)",
                               range = readxl::cell_cols("A:C"),
                               col_names = c("C1", "C2", "Res"))
cont <- numeric_rows(as.data.frame(cont_raw), 1:3)
write.csv(cont, "tests/testthat/fixtures/binary_mps_cpf_continuous.csv",
          row.names = FALSE)
cat("wrote", nrow(cont), "continuous rows\n")

bin_raw <- readxl::read_excel(wb, sheet = "CA - binary data model (new)",
                              range = readxl::cell_cols("A:D"),
                              col_names = c("C1", "C2", "Affected", "Exposed"))
bin <- numeric_rows(as.data.frame(bin_raw), 1:4)
write.csv(bin, "tests/testthat/fixtures/binary_mps_cpf_quantal.csv",
          row.names = FALSE)
cat("wrote", nrow(bin), "quantal rows\n")
