# Table-driven validation of the engine against the published MixTox/Excel
# workbooks. `inst/validation/oracles.csv` is the single source of truth: it
# carries, per fitted quantity, the workbook value, the value this engine must
# produce, the tolerance, and where in the workbook the number came from. Both
# the validation tests and the validation report read that one table, so a
# reported number and an asserted number can never drift apart.
#
# Adding a newly published dataset is a data-entry job: drop the CSV fixture in
# tests/testthat/fixtures/, register how to fit it in `validation_datasets()`,
# and add its rows to oracles.csv. No new test code.

#' Validation oracle table
#'
#' The reference targets that [validation_report()] checks the engine against.
#' Columns: `dataset`, `stage`, `parameter`, `label`, `reference` (the
#' externally published value), `expected` (what this engine must produce),
#' `tolerance` (relative, as in [testthat::expect_equal()]), `kind`,
#' `source_type`, `source`, `note`.
#'
#' `source_type` says where the reference value came from: `workbook` for a cell
#' in the MixTox/Excel workbook the fixture was extracted from, and
#' `publication` for a number printed in the peer-reviewed paper. Publication
#' values are the stronger evidence — they are citable and independent of the
#' spreadsheet — but only some quantities are comparable across the two, because
#' the published single-compound curves are fitted per compound with their own
#' asymptote while the mixture model shares one. Where that makes a quantity
#' incomparable, the row says so in its note.
#'
#' `kind` says what sort of claim the row makes:
#' \describe{
#'   \item{`match`}{The engine must reproduce the workbook value.}
#'   \item{`divergent`}{The engine deliberately differs; `note` says why.}
#'   \item{`pinned`}{Held fixed in both; agreement is by construction.}
#'   \item{`unidentified`}{The data do not pin this parameter down; only a
#'     qualitative bound is claimed.}
#'   \item{`bounded`}{Only an inequality against the workbook is claimed.}
#' }
#' @return A data frame.
#' @export
validation_oracles <- function() {
  path <- system.file("validation", "oracles.csv", package = "mixdra")
  if (!nzchar(path)) stop("validation oracle table not found", call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Provenance of the experiments behind the validation datasets
#'
#' One row per laboratory experiment: the organism, life stage, endpoint,
#' exposure duration, test guideline, soil, units and mixture ratios, plus the
#' publication and open data repository the experiment was reported in, and any
#' caveats a reader needs in order to interpret its reference values (separated
#' by `|`). Several fixtures come from the same experiment, so the biology is
#' recorded once here rather than repeated per dataset.
#' @return A data frame, one row per experiment.
#' @export
validation_experiments <- function() {
  path <- system.file("validation", "datasets.csv", package = "mixdra")
  if (!nzchar(path)) stop("validation dataset table not found", call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' How each validation dataset is fitted
#'
#' One entry per `dataset` id used in [validation_oracles()]: the fixture path
#' (relative to the test fixtures root), the experiment it came from (a row of
#' [validation_experiments()]), the measured endpoint, and the call that
#' produces the fit. Kept next to the oracle table so a new dataset is
#' registered in one place.
#' @return A named list of `list(path, experiment, endpoint, title, fit)` entries.
#' @keywords internal
validation_datasets <- function() {
  list(
    binary_mps_cpf_ca_continuous = list(
      path  = file.path("binary", "cpf_mps_imi",
                        "binary_ca_mps_cpf_imi_continuous.csv"),
      experiment = "experiment_2",
      endpoint   = "reproduction (number of juveniles)",
      title = "Binary CPF + MPs, concentration addition, continuous response",
      fit   = function(df) analyse_mixture(df, reference = "CA",
                                           response = "continuous", n_starts = 1)),
    binary_mps_cpf_quantal = list(
      path  = file.path("binary", "survival", "binary_mps_cpf_quantal.csv"),
      experiment = "experiment_2",
      # Survival was recorded but not analysed in the publication, so this
      # dataset has workbook reference values only.
      endpoint   = "adult survival (quantal)",
      title = "Binary CPF + MPs, concentration addition, quantal (survival)",
      fit   = function(df) analyse_mixture(df, reference = "CA",
                                           response = "binary", n_starts = 1)),
    ternary_fbsa_cpf_imi = list(
      path  = file.path("ternary", "fbsa_cpf_imi",
                        "ternary_fbsa_cpf_imi_continuous.csv"),
      experiment = "experiment_1",
      endpoint   = "reproduction (number of juveniles)",
      title = "Ternary CPF + FBSA + IMI, Advanced S/A, continuous response",
      fit   = function(df) analyse_ternary(df, "CA", "continuous", n_starts = 10,
                                           time_limit = 120,
                                           lower = c(ec50_2 = 5.579),
                                           upper = c(ec50_2 = 5.581)))
  )
}

#' Pull one fitted quantity out of a result object
#'
#' Maps an oracle row's (`stage`, `parameter`) pair onto the corresponding slot
#' of an [analyse_mixture()] or [analyse_ternary()] result. Returns `NA` rather
#' than erroring when a slot is absent, so one missing quantity does not sink a
#' whole report.
#' @keywords internal
validation_extract <- function(res, stage, parameter) {
  pick <- function(x, nm) if (!is.null(x) && nm %in% names(x)) unname(x[[nm]]) else NA
  switch(stage,
    reference = if (identical(parameter, "objective"))
                  res$fits$reference$objective %||% NA
                else pick(as.list(res$fits$reference$par), parameter),
    selection = if (identical(parameter, "chosen")) res$chosen %||% NA else NA,
    base      = pick(as.list(res$base), parameter),
    pairwise  = pick(as.list(res$pairwise), parameter),
    overall   = if (identical(parameter, "objective"))
                  res$fits$overall$objective %||% NA else NA,
    NA)
}

#' Where the validation fixture CSVs live
#'
#' The fixtures are test data, so they live under `tests/testthat/fixtures/`
#' rather than in the installed package. Resolution order: the
#' `mixdra.validation_fixtures` option, an installed copy under
#' `inst/validation/fixtures/`, then the source tree relative to the working
#' directory.
#' @return A path, or `""` when no fixture root can be found.
#' @export
validation_fixtures <- function() {
  opt <- getOption("mixdra.validation_fixtures")
  if (!is.null(opt) && dir.exists(opt)) return(opt)
  inst <- system.file("validation", "fixtures", package = "mixdra")
  if (nzchar(inst) && dir.exists(inst)) return(inst)
  for (p in c("tests/testthat/fixtures", "../../tests/testthat/fixtures",
              "../tests/testthat/fixtures"))
    if (dir.exists(p)) return(normalizePath(p, winslash = "/"))
  ""
}

#' Compare fitted results against the workbook oracle table
#'
#' Runs (or accepts) the fits for each registered validation dataset and joins
#' them onto [validation_oracles()], adding the observed value, the relative
#' difference, and a verdict. This is what both the validation tests and the
#' validation report consume.
#'
#' @param datasets Character vector of dataset ids to include; default all.
#' @param fits Optional named list of already-computed fit results, keyed by
#'   dataset id. Supply this to avoid re-running the slow fits (the ternary fit
#'   alone takes ~2 minutes).
#' @param fixtures Root of the test fixtures tree; see [validation_fixtures()].
#' @return The oracle table plus `observed`, `rel_diff`, `rel_diff_reference`
#'   and `verdict` columns.
#'   `verdict` is `"pass"`, `"fail"`, or `"info"` for rows that make no
#'   point-value claim. `rel_diff` is the distance from `expected` (the
#'   assertion target); `rel_diff_reference` the distance from the published
#'   value, which differs for `divergent` rows.
#' @export
validation_report <- function(datasets = NULL, fits = NULL,
                              fixtures = validation_fixtures()) {
  orc <- validation_oracles()
  reg <- validation_datasets()
  if (is.null(datasets)) datasets <- unique(orc$dataset)
  orc <- orc[orc$dataset %in% datasets, , drop = FALSE]

  for (id in datasets) {
    if (!is.null(fits[[id]])) next
    spec <- reg[[id]]
    if (is.null(spec)) stop("unregistered validation dataset: ", id, call. = FALSE)
    csv <- file.path(fixtures, spec$path)
    if (!file.exists(csv)) {
      fits[[id]] <- NULL
      next
    }
    fits[[id]] <- spec$fit(utils::read.csv(csv))
  }

  obs <- vapply(seq_len(nrow(orc)), function(i) {
    res <- fits[[orc$dataset[i]]]
    v <- if (is.null(res)) NA else
      validation_extract(res, orc$stage[i], orc$parameter[i])
    if (length(v) != 1) NA_character_ else as.character(v)
  }, character(1))

  orc$observed <- obs
  num_obs <- suppressWarnings(as.numeric(orc$observed))
  num_exp <- suppressWarnings(as.numeric(orc$expected))
  num_wb  <- suppressWarnings(as.numeric(orc$reference))
  orc$rel_diff <- ifelse(is.na(num_obs) | is.na(num_exp) | num_exp == 0,
                         NA_real_, abs(num_obs - num_exp) / abs(num_exp))
  # Distance from the *published* value, which for `divergent` rows is not the
  # assertion target: it is how far the deliberate methodological difference
  # actually moves the number, and belongs in the report.
  orc$rel_diff_reference <- ifelse(is.na(num_obs) | is.na(num_wb) | num_wb == 0,
                                  NA_real_, abs(num_obs - num_wb) / abs(num_wb))

  orc$verdict <- vapply(seq_len(nrow(orc)), function(i) {
    if (is.na(orc$observed[i])) return("missing")
    if (orc$kind[i] %in% c("unidentified", "bounded")) return("info")
    if (is.na(orc$expected[i])) return("info")
    if (is.na(num_exp[i]))                       # character comparison
      return(if (identical(orc$observed[i], orc$expected[i])) "pass" else "fail")
    if (is.na(orc$tolerance[i])) return("info")
    if (is.na(orc$rel_diff[i])) return("info")
    if (orc$rel_diff[i] <= as.numeric(orc$tolerance[i])) "pass" else "fail"
  }, character(1))

  orc
}
