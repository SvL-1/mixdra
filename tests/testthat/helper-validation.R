# Shared, memoised fits for the workbook-validation tests. The ternary fit alone
# is ~2 min on 419 rows with 10 starts, and three test files need it; caching it
# per R session means CI pays for each dataset once instead of once per file.
# validation_run() seeds each fit, so the multi-start optimiser is
# reproducible and the tests fit exactly what the report reports.

.validation_cache <- new.env(parent = emptyenv())

validation_fit <- function(id) {
  if (!is.null(.validation_cache[[id]])) return(.validation_cache[[id]])
  spec <- validation_datasets()[[id]]
  if (is.null(spec)) stop("unregistered validation dataset: ", id)
  csv <- testthat::test_path("fixtures", spec$path)
  if (!file.exists(csv)) return(NULL)
  .validation_cache[[id]] <- validation_run(spec, utils::read.csv(csv))
  .validation_cache[[id]]
}

# Run the oracle table for one dataset against its cached fit, and assert every
# row that makes a point-value claim. Rows of kind "info"/"missing" are reported
# by validation_report() but not asserted here -- their claims are inequalities
# or sign conditions, asserted explicitly in the calling test.
expect_matches_reference <- function(id) {
  fit <- validation_fit(id)
  testthat::skip_if(is.null(fit), paste0("fixture for ", id, " missing"))
  rep <- validation_report(datasets = id, fits = stats::setNames(list(fit), id))
  for (i in seq_len(nrow(rep))) {
    if (!rep$verdict[i] %in% c("pass", "fail")) next
    testthat::expect_equal(
      rep$verdict[i], "pass",
      info = sprintf("%s / %s: reference %s, expected %s, got %s (rel. diff %.3g)",
                     id, rep$label[i], rep$reference[i], rep$expected[i],
                     rep$observed[i], rep$rel_diff[i]))
  }
  invisible(rep)
}
