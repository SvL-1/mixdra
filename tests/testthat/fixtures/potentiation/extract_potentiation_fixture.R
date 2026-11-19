# One-off generator: extract the (de)potentiation example fixture from Sam's
# draft workbook "DRCs different ratios FBSA with IMI and CPF.xlsx" (sheet
# "75-25 drc CPF on FBSA IMI"). The workbook is NOT committed; pass its path.
# Run from the repo root (an optional second argument overrides the output):
#   Rscript tests/testthat/fixtures/potentiation/extract_potentiation_fixture.R \
#     "C:/path/to/DRCs different ratios FBSA with IMI and CPF.xlsx"
# Requires `readxl`. Writes potentiation_example_continuous.csv next to this
# script (columns Series, C1, Res); copy it to inst/extdata/ unchanged.
#
# The sheet holds three side-by-side comparison blocks. Each block lists the
# stressor-1-alone reference curve (rows 30-74: 20 controls + 25 dosed rows)
# and one "with stressor 2" curve (rows 81-125: the same 20 controls repeated,
# then 25 dosed rows). Concentration / response columns:
#   reference  M30:N74   (identical copies in W30:X74 and AG30:AH74)
#   Level 1    M101:N125 (block "25%", columns L-T)
#   Level 2    AG101:AH125 (block "50%", columns AF-AN)
#   Level 3    W101:X125 (block "75%", columns V-AD)
# The series' control rows (M81:N100, AG81:AH100, W81:X100) are the reference
# controls again, so they are not written: the engine borrows the reference
# controls for a series that has none, which reproduces the workbook.
#
# Labels are deliberately generic ("Stressor 1", "Level 1".."Level 3"): the
# workbook's own stressor and block names are unreliable (e.g. a block label
# names a stressor that is not in the experiment).
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required to regenerate this fixture.")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: extract_potentiation_fixture.R <path to the workbook .xlsx>")
wb <- args[[1]]
sheet <- "75-25 drc CPF on FBSA IMI"

cells <- function(range) {
  x <- readxl::read_excel(wb, sheet = sheet, range = range,
                          col_names = c("C1", "Res"))
  x <- as.data.frame(x)
  stopifnot(is.numeric(x$C1), is.numeric(x$Res), !anyNA(x))
  x
}

ref <- cells("M30:N74")
stopifnot(nrow(ref) == 45, sum(ref$C1 == 0) == 20, all(ref$C1[1:20] == 0))

# (a) every block repeats the same reference rows
stopifnot(identical(cells("W30:X74"), ref), identical(cells("AG30:AH74"), ref))

# (b) every series' control rows are the reference controls (N30:N49)
for (r in c("M81:N100", "AG81:AH100", "W81:X100")) {
  ctl <- cells(r)
  stopifnot(all(ctl$C1 == 0), identical(ctl$Res, ref$Res[1:20]))
}

series <- list("Stressor 1" = ref,
               "Level 1"    = cells("M101:N125"),
               "Level 2"    = cells("AG101:AH125"),
               "Level 3"    = cells("W101:X125"))
stopifnot(all(vapply(series[-1], function(s) all(s$C1 > 0), logical(1))))

out <- do.call(rbind, Map(function(lab, s) data.frame(Series = lab, s),
                          names(series), series))
rownames(out) <- NULL
stopifnot(nrow(out) == 120)

dest <- if (length(args) >= 2) args[[2]] else
  "tests/testthat/fixtures/potentiation/potentiation_example_continuous.csv"
write.csv(out, dest, row.names = FALSE)
cat(sprintf("wrote %d rows -> %s\n", nrow(out), dest))
print(table(Series = out$Series, control = out$C1 == 0))
