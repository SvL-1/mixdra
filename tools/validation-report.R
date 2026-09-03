# Regenerate the validation report.
#
#   Rscript tools/validation-report.R [html-output]
#
# Re-fits every registered validation dataset and writes the comparison against
# the published workbook values three ways: VALIDATION.md (committed, renders on
# GitHub), the badge block in README.md, and a standalone HTML page (not
# committed; CI attaches it to the run). Slow -- the ternary fit alone is a few
# minutes on 419 rows with 10 starts.
#
# CI runs this and then diffs VALIDATION.md and README.md, so the committed
# numbers cannot drift away from what the engine fits.

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
html <- if (length(args)) args[[1]] else "validation-report.html"

options(mixdra.validation_fixtures = normalizePath("tests/testthat/fixtures",
                                                   winslash = "/"))
set.seed(42)
rep <- validation_report()

# A failed row must not quietly ship inside a document that reads as evidence.
bad <- rep[rep$verdict %in% c("fail", "missing"), , drop = FALSE]
if (nrow(bad)) {
  message("Validation rows that did not pass:")
  print(bad[, c("dataset", "label", "workbook", "observed", "verdict")],
        row.names = FALSE)
}

validation_report_md(rep, "VALIDATION.md")
validation_badges(rep, "README.md")
validation_report_html(rep, html)

message(sprintf("%d quantities, %d passing | wrote VALIDATION.md, README badges, %s",
                nrow(rep), sum(rep$verdict == "pass"), html))

# The job summary makes the table readable straight from the Actions run page.
if (nzchar(Sys.getenv("GITHUB_STEP_SUMMARY")))
  writeLines(validation_report_md(rep, NULL),
             Sys.getenv("GITHUB_STEP_SUMMARY"))

if (nrow(bad)) quit(status = 1)
