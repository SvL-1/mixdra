# Pure helpers for the Shiny app's data layer. No Shiny, no plotly here, so these
# are fully unit-testable and run without the UI stack installed.

#' Fixed column schema for a stage and response type
#' @param stage "single" or "binary".
#' @param response "continuous" or "quantal".
#' @return Character vector of required column names.
#' @keywords internal
upload_schema <- function(stage, response) {
  stage <- match.arg(stage, c("single", "binary"))
  response <- match.arg(response, c("continuous", "quantal"))
  if (stage == "single") {
    if (response == "continuous") c("Conc", "Res") else c("Conc", "Affected", "Exposed")
  } else {
    if (response == "continuous") c("C1", "C2", "Res") else c("C1", "C2", "Affected", "Exposed")
  }
}

#' Example template data frame for a stage and response type
#'
#' Returns a small, illustrative dataset with exactly the schema columns. The
#' binary template includes single-chemical rows (one chemical at 0) because
#' [mixdra::analyse_mixture()] seeds itself from them.
#' @inheritParams upload_schema
#' @return A data frame the user can download, fill in, and re-upload.
#' @keywords internal
template_df <- function(stage, response) {
  cols <- upload_schema(stage, response)
  if (stage == "single") {
    conc <- c(0, 0.1, 0.3, 1, 3, 10)
    if (response == "continuous") {
      data.frame(Conc = conc, Res = c(100, 96, 82, 50, 18, 4))
    } else {
      data.frame(Conc = conc, Affected = c(0, 1, 2, 5, 8, 10), Exposed = rep(10, 6))
    }
  } else {
    # single-chemical series for each chemical + a few mixture rows
    c1 <- c(0, 0.1, 0.3, 1, 0, 0, 0, 0.1, 0.3, 1)
    c2 <- c(0, 0,   0,   0, 0.1, 0.3, 1, 0.1, 0.3, 1)
    if (response == "continuous") {
      data.frame(C1 = c1, C2 = c2,
                 Res = c(100, 80, 55, 20, 88, 70, 35, 72, 45, 12))
    } else {
      data.frame(C1 = c1, C2 = c2,
                 Affected = c(0, 2, 4, 8, 1, 3, 7, 3, 6, 9), Exposed = rep(10, 10))
    }
  }
}

#' Validate an uploaded data frame against the fixed schema
#'
#' @inheritParams upload_schema
#' @param df The uploaded data frame.
#' @return Character vector of human-readable error messages; empty if valid.
#' @keywords internal
validate_upload <- function(df, stage, response) {
  errs <- character(0)
  req_cols <- upload_schema(stage, response)

  missing <- setdiff(req_cols, names(df))
  if (length(missing))
    errs <- c(errs, paste0("Missing column(s): ", paste(missing, collapse = ", "),
                           ". Expected exactly: ", paste(req_cols, collapse = ", "), "."))

  present <- intersect(req_cols, names(df))
  non_num <- present[!vapply(df[present], is.numeric, logical(1))]
  if (length(non_num))
    errs <- c(errs, paste0("Non-numeric column(s): ", paste(non_num, collapse = ", "), "."))

  # Range checks only on numeric columns that are present.
  conc_cols <- intersect(c("Conc", "C1", "C2"), present)
  conc_ok <- conc_cols[vapply(df[conc_cols], is.numeric, logical(1))]
  if (length(conc_ok) && any(unlist(df[conc_ok]) < 0, na.rm = TRUE))
    errs <- c(errs, "Concentrations must be >= 0.")

  if (response == "quantal" && all(c("Affected", "Exposed") %in% present) &&
      is.numeric(df$Affected) && is.numeric(df$Exposed)) {
    if (any(df$Affected < 0 | df$Exposed < 0, na.rm = TRUE))
      errs <- c(errs, "Affected and Exposed must be non-negative.")
    if (any(df$Affected > df$Exposed, na.rm = TRUE))
      errs <- c(errs, "Affected must be <= Exposed.")
  }

  if (stage == "single" && "Conc" %in% present && is.numeric(df$Conc)) {
    n_distinct <- length(unique(df$Conc[!is.na(df$Conc)]))
    if (n_distinct < 4)
      errs <- c(errs, "Need at least 4 distinct concentrations to fit a single-chemical curve.")
  }

  errs
}

#' Read an uploaded CSV file
#' @param path File path.
#' @return A data frame.
#' @keywords internal
read_upload <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
}

#' Map an uploaded data frame to the engine's column convention
#'
#' The single-chemical template uses `Conc`; the engine expects `C1`.
#' @param df Uploaded data frame.
#' @param stage "single" or "binary".
#' @return The data frame with engine-ready column names.
#' @keywords internal
to_engine_df <- function(df, stage) {
  if (stage == "single") names(df)[names(df) == "Conc"] <- "C1"
  df
}

#' One chemical's single-compound series from a binary frame
#'
#' Keeps the rows where the *other* chemical's concentration is 0 (so the shared
#' control row is included), drops the other concentration column, and renames
#' this chemical's concentration column to `C1`. The result has the shape a
#' single-chemical fitter expects (`C1` + response columns).
#' @param df Binary engine data frame (`C1`, `C2`, response columns).
#' @param chem 1 or 2 — which chemical's marginal series to extract.
#' @return A data frame with `C1` and the response columns.
#' @keywords internal
marginal_df <- function(df, chem) {
  this  <- paste0("C", chem)
  other <- paste0("C", if (chem == 1) 2 else 1)
  out <- df[df[[other]] == 0, , drop = FALSE]
  out[[other]] <- NULL
  names(out)[names(out) == this] <- "C1"
  rownames(out) <- NULL
  out
}

#' Frozen curve-parameter vector from two single-chemical fits
#'
#' Builds the named vector `analyse_mixture(start = …)` holds fixed: a shared
#' `max` (the average of the two per-chemical fits, matching the engine's
#' [seed_from_singles()] behaviour) plus per-chemical `slope1/slope2` and
#' `ec501/ec502`. Names match the binary registry's base parameters.
#' @param fit1,fit2 Single-fit results (each a list with `par = c(max, slope, ec50)`).
#' @return A named numeric vector: `max`, `slope1`, `slope2`, `ec501`, `ec502`.
#' @keywords internal
assemble_curve_params <- function(fit1, fit2) {
  c(max    = mean(c(fit1$par[["max"]], fit2$par[["max"]])),
    slope1 = fit1$par[["slope"]],
    slope2 = fit2$par[["slope"]],
    ec501  = fit1$par[["ec50"]],
    ec502  = fit2$par[["ec50"]])
}

#' Assemble lower/upper bound vectors from Advanced-panel inputs
#'
#' Reads `lo_<param>` / `hi_<param>` values for `params`; blank/NA entries are
#' dropped. To FIX a parameter, set its lower and upper to the same value.
#' Callers pass the parameter set they need (the curve-fit panel passes
#' `max`/`slope`/`ec50`); the default is the binary base set.
#' @param values Named list (e.g. a Shiny `input`) holding `lo_*`/`hi_*` numbers.
#' @param params Character vector of parameter names to read.
#' @return A list with `lower` and `upper` named numeric vectors (or NULL).
#' @keywords internal
collect_bounds <- function(values,
                           params = c("max", "slope1", "slope2", "ec501", "ec502")) {
  pick <- function(prefix) {
    v <- vapply(params, function(p) {
      x <- values[[paste0(prefix, p)]]
      if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else as.numeric(x)
    }, numeric(1))
    v <- v[!is.na(v)]
    if (length(v)) v else NULL
  }
  list(lower = pick("lo_"), upper = pick("hi_"))
}
