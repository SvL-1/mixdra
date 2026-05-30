# One-off: extract the ternary FBSA+CPF+IMI raw data into a CSV fixture.
# Run from the repo root:
#   R -q -e 'source("tests/testthat/fixtures/extract_ternary_fixture.R")'
# Layout (verified): workbook "FBSA CPF IMI ternary -simplified_correct.xls",
# sheet "CA CPF FBSA IMI", row-32 header CPF | FBSA | IMI | RGR. Columns:
# A=CPF->C1, B=FBSA->C2, C=IMI->C3, D=RGR(Reproduction)->Res. Data from row 33.
# We read whole columns A:D and keep only rows where all four cells are numeric,
# which robustly drops the header/parameter rows and any blank separator rows.
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required to regenerate this fixture.")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
wb <- "FBSA CPF IMI ternary -simplified_correct.xls"

numeric_rows <- function(df, cols) {
  vals <- lapply(cols, function(j) suppressWarnings(as.numeric(df[[j]])))
  keep <- Reduce(`&`, lapply(vals, function(v) !is.na(v)))
  out <- as.data.frame(vals, col.names = c("C1", "C2", "C3", "Res"))
  out[keep, , drop = FALSE]
}

raw <- readxl::read_excel(wb, sheet = "CA CPF FBSA IMI",
                          range = readxl::cell_cols("A:D"),
                          col_names = c("C1", "C2", "C3", "Res"))
dat <- numeric_rows(as.data.frame(raw), 1:4)
write.csv(dat, "tests/testthat/fixtures/ternary_fbsa_cpf_imi_continuous.csv",
          row.names = FALSE)
cat("wrote", nrow(dat), "ternary rows\n")
