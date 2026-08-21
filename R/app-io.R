# Pure helpers for the Shiny app's data layer. No Shiny, no plotly here, so these
# are fully unit-testable and run without the UI stack installed.

#' Fixed column schema for a stage and response type
#' @param stage "single", "binary", "ternary", or "campaign".
#' @param response "continuous" or "quantal".
#' @return Character vector of required column names.
#' @keywords internal
upload_schema <- function(stage, response) {
  stage <- match.arg(stage, c("single", "binary", "ternary", "campaign"))
  response <- match.arg(response, c("continuous", "quantal"))
  if (stage == "single") {
    if (response == "continuous") c("Conc", "Res") else c("Conc", "Affected", "Exposed")
  } else if (stage == "binary" || stage == "campaign") {
    # A campaign REQUIRES C1/C2; C3 is optional and makes it a three-stressor
    # campaign. Optionality is enforced in validate_upload(), not the schema.
    if (response == "continuous") c("C1", "C2", "Res") else c("C1", "C2", "Affected", "Exposed")
  } else {
    if (response == "continuous") c("C1", "C2", "C3", "Res")
    else c("C1", "C2", "C3", "Affected", "Exposed")
  }
}

#' Example template data frame for a stage and response type
#'
#' Returns a small, illustrative dataset with exactly the schema columns. The
#' binary template includes single-chemical rows (one chemical at 0) because
#' [mixdra::analyse_mixture()] seeds itself from them. The campaign template is
#' the full ternary template with all four strata.
#' @inheritParams upload_schema
#' @return A data frame the user can download, fill in, and re-upload.
#' @keywords internal
template_df <- function(stage, response) {
  # A campaign template IS the ternary template: all four strata in one frame.
  if (stage == "campaign") stage <- "ternary"
  cols <- upload_schema(stage, response)
  if (stage == "single") {
    conc <- c(0, 0.1, 0.3, 1, 3, 10)
    if (response == "continuous") {
      data.frame(Conc = conc, Res = c(100, 96, 82, 50, 18, 4))
    } else {
      data.frame(Conc = conc, Affected = c(0, 1, 2, 5, 8, 10), Exposed = rep(10, 6))
    }
  } else if (stage == "binary") {
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
  } else {  # ternary: all tiers so the staged fit has data at every stage
    rows <- rbind(
      data.frame(C1 = 0,           C2 = 0,           C3 = 0),            # control
      data.frame(C1 = c(0.1, 0.3, 1), C2 = 0,        C3 = 0),           # chem1 single
      data.frame(C1 = 0,           C2 = c(0.1, 0.3, 1), C3 = 0),        # chem2 single
      data.frame(C1 = 0,           C2 = 0,           C3 = c(0.1, 0.3, 1)), # chem3 single
      data.frame(C1 = c(0.5, 0.5, 0), C2 = c(0.5, 0, 0.5), C3 = c(0, 0.5, 0.5)), # 3 binaries
      data.frame(C1 = c(0.2, 0.4, 0.6), C2 = c(0.2, 0.4, 0.6), C3 = c(0.2, 0.4, 0.6)) # 1:1:1 ternary
    )
    tot <- rows$C1 + rows$C2 + rows$C3
    if (response == "continuous") {
      cbind(rows, Res = round(100 / (1 + tot), 1))          # illustrative decline
    } else {
      cbind(rows, Affected = pmin(round(10 * tot / (1 + tot)), 10), Exposed = 10)
    }
  }
}

#' Validate an uploaded data frame against the fixed schema
#'
#' @param stage "single", "binary", "ternary", or "campaign".
#' @param response "continuous" or "quantal".
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
  conc_cols <- intersect(c("Conc", "C1", "C2", "C3"), present)
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

  if (stage == "ternary" && all(c("C1", "C2", "C3") %in% present) &&
      all(vapply(df[c("C1", "C2", "C3")], is.numeric, logical(1)))) {
    nz <- rowSums(as.matrix(df[c("C1", "C2", "C3")]) > 0)
    if (!any(nz == 3, na.rm = TRUE))
      errs <- c(errs, paste0("No ternary rows (all of C1, C2, C3 > 0); the ",
                             "per-ratio A4 step needs at least one ternary mixture."))
  }

  if (stage == "campaign") {
    cols <- intersect(c("C1", "C2", "C3"), present)
    ok   <- cols[vapply(df[cols], is.numeric, logical(1))]
    if (length(ok) >= 2) {
      chems <- campaign_chems(df[ok])
      if (length(chems) < 2) {
        errs <- c(errs, paste0(
          "A campaign needs at least two stressors with a positive ",
          "concentration; found ", length(chems), "."))
      } else {
        for (k in chems) {
          nd <- length(unique(stats::na.omit(single_df(df, k)$C1)))
          if (nd < 4)
            errs <- c(errs, paste0(
              "Stressor ", k, " has only ", nd, " distinct concentrations in ",
              "its single-stressor series; at least 4 are needed to fit a curve."))
        }
      }
    }
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
#' @param stage "single", "binary", or "ternary".
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

#' One chemical's single-compound series from a ternary frame
#'
#' Keeps rows where the OTHER TWO chemicals are 0 (so the shared control row is
#' included), drops their columns, and renames this chemical's concentration
#' column to `C1` — the shape a single-chemical fitter expects. The ternary
#' analogue of [marginal_df()].
#' @param df Ternary engine data frame (`C1`, `C2`, `C3`, response columns).
#' @param chem 1, 2 or 3 — which chemical's marginal series to extract.
#' @return A data frame with `C1` and the response columns.
#' @keywords internal
marginal_df3 <- function(df, chem) {
  cols  <- paste0("C", 1:3)
  this  <- paste0("C", chem)
  other <- setdiff(cols, this)
  keep  <- df[[other[1]]] == 0 & df[[other[2]]] == 0
  out   <- df[keep, , drop = FALSE]
  out[other] <- NULL
  names(out)[names(out) == this] <- "C1"
  rownames(out) <- NULL
  out
}

#' Frozen ternary base-parameter vector from three single-chemical fits
#'
#' Builds the named vector [analyse_ternary()] holds fixed as its `base`: a
#' shared `max` (mean of the three per-chemical fits, matching the engine's
#' seeding) plus per-chemical `slope1/2/3` and `ec50_1/2/3`. Names match the
#' ternary registry's base parameters (underscore `ec50_i`, unlike binary's
#' `ec501`).
#' @param fit1,fit2,fit3 Single-fit results (each `par = c(max, slope, ec50)`).
#' @return A named numeric vector: `max`, `slope1-3`, `ec50_1-3`.
#' @keywords internal
assemble_curve_params3 <- function(fit1, fit2, fit3) {
  c(max    = mean(c(fit1$par[["max"]], fit2$par[["max"]], fit3$par[["max"]])),
    slope1 = fit1$par[["slope"]],
    slope2 = fit2$par[["slope"]],
    slope3 = fit3$par[["slope"]],
    ec50_1 = fit1$par[["ec50"]],
    ec50_2 = fit2$par[["ec50"]],
    ec50_3 = fit3$par[["ec50"]])
}

#' Which stressors a campaign frame actually doses
#'
#' A concentration column that is present but never positive (e.g. a `C3` of
#' zeros pasted in by mistake) does not count — the campaign is then a
#' two-stressor one.
#' @param df Campaign data frame with `C1`, `C2` and optionally `C3`.
#' @return Integer vector of stressor indices, e.g. `c(1L, 2L)`.
#' @keywords internal
campaign_chems <- function(df) {
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  act  <- cols[vapply(df[cols], function(x) any(x > 0, na.rm = TRUE), logical(1))]
  as.integer(sub("^C", "", act))
}

#' Number of stressors in a campaign frame
#' @inheritParams campaign_chems
#' @return 2 or 3 (or fewer, which [validate_upload()] rejects).
#' @keywords internal
campaign_n_chem <- function(df) length(campaign_chems(df))

#' Rows where every stressor outside `keep` is absent
#'
#' The shared core of [single_df()] and [pair_df()]: keeps the rows where every
#' concentration column NOT in `keep` is 0 (so the shared control row always
#' comes along), then drops those columns. Renaming is left to the caller,
#' which knows whether it wants `C1` or `C1`/`C2`.
#' @param df Campaign data frame.
#' @param keep Character vector of concentration columns to retain.
#' @return A data frame with the non-kept concentration columns removed.
#' @keywords internal
slice_to_chems <- function(df, keep) {
  cols  <- intersect(c("C1", "C2", "C3"), names(df))
  other <- setdiff(cols, keep)
  rows  <- if (length(other))
    Reduce(`&`, lapply(other, function(k) df[[k]] == 0)) else rep(TRUE, nrow(df))
  out <- df[rows, , drop = FALSE]
  out[other] <- NULL
  rownames(out) <- NULL
  out
}

#' One stressor's single-stressor series from a campaign frame
#'
#' Keeps the rows where every OTHER stressor is 0 (so the shared control row is
#' included), drops their columns, and renames this stressor's concentration
#' column to `C1` — the shape a single-stressor fitter expects. Generalises
#' [marginal_df()] / [marginal_df3()] to 2- or 3-stressor frames.
#' @inheritParams campaign_chems
#' @param chem 1, 2 or 3 — which stressor's series to extract.
#' @return A data frame with `C1` and the response columns.
#' @keywords internal
single_df <- function(df, chem) {
  this <- paste0("C", chem)
  out <- slice_to_chems(df, this)
  names(out)[names(out) == this] <- "C1"
  out
}

#' One pair's rows from a campaign frame
#'
#' Keeps the rows where the third stressor is 0 (so the singles and the shared
#' control come along, exactly as the binary fitter expects), drops its column,
#' and renames the pair's concentration columns to `C1`/`C2`.
#' @inheritParams campaign_chems
#' @param i,j Stressor indices of the pair, `i < j`.
#' @return A data frame with `C1`, `C2` and the response columns.
#' @keywords internal
pair_df <- function(df, i, j) {
  stopifnot(i < j)
  this <- paste0("C", c(i, j))
  out <- slice_to_chems(df, this)
  names(out)[names(out) == this[1]] <- "C1"
  names(out)[names(out) == this[2]] <- "C2"
  out
}

#' Map a campaign base onto one pair's binary parameter names
#'
#' A three-stressor campaign base is named for the ternary registry
#' (`slope1..3`, `ec50_1..3`), but each pair is fitted with the two-chemical
#' model, whose registry uses `slope1`, `slope2`, `ec501`, `ec502`. This selects
#' the pair's two stressors and renames them into that binary shape, so the
#' frozen curves are actually held fixed. A two-stressor campaign base is
#' already in binary shape and passes through unchanged.
#' @param base Campaign base parameters, from [campaign_base()].
#' @param i,j Stressor indices of the pair, `i < j`.
#' @return A named numeric vector: `max`, `slope1`, `slope2`, `ec501`, `ec502`.
#' @keywords internal
pair_base <- function(base, i, j) {
  stopifnot(i < j)
  if (all(c("ec501", "ec502") %in% names(base))) return(base)   # already binary
  c(max    = unname(base[["max"]]),
    slope1 = unname(base[[paste0("slope", i)]]),
    slope2 = unname(base[[paste0("slope", j)]]),
    ec501  = unname(base[[paste0("ec50_", i)]]),
    ec502  = unname(base[[paste0("ec50_", j)]]))
}

#' Assemble lower/upper bound vectors from Advanced-panel inputs
#'
#' Reads `lo_<param>` / `hi_<param>` values for `params`; blank/NA entries are
#' dropped. To FIX a parameter, set its lower and upper to the same value.
#' Callers pass the parameter set they need (the curve-fit panel passes
#' `max`/`slope`/`ec50`); the default is the binary base set.
#' @param values Named list (e.g. a Shiny `input`) holding `lo_*`/`hi_*` numbers.
#' @param params Character vector of parameter names to read.
#' @param lo_prefix,hi_prefix Input-id prefix for lower/upper bounds
#'   (default `"lo_"` / `"hi_"`).
#' @return A list with `lower` and `upper` named numeric vectors (or NULL).
#' @keywords internal
collect_bounds <- function(values,
                           params = c("max", "slope1", "slope2", "ec501", "ec502"),
                           lo_prefix = "lo_", hi_prefix = "hi_") {
  pick <- function(prefix) {
    v <- vapply(params, function(p) {
      x <- values[[paste0(prefix, p)]]
      if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else as.numeric(x)
    }, numeric(1))
    v <- v[!is.na(v)]
    if (length(v)) v else NULL
  }
  list(lower = pick(lo_prefix), upper = pick(hi_prefix))
}

#' Read the 7-parameter lower/upper grid from the Optimize-all panel
#'
#' Thin wrapper over [collect_bounds()] for the full binary parameter set
#' (curve params + interaction `a`/`b`), reading the Optimize-all panel's
#' `olo_*` / `ohi_*` inputs.
#' @param values Named list (e.g. a Shiny `input`) holding `olo_*`/`ohi_*`.
#' @param params Parameter names to read (default: the binary 7-set).
#' @return A list with `lower` and `upper` named numeric vectors (or NULL).
#' @keywords internal
collect_bounds_all <- function(values,
                               params = c("max", "slope1", "slope2",
                                          "ec501", "ec502", "a", "b")) {
  collect_bounds(values, params, lo_prefix = "olo_", hi_prefix = "ohi_")
}

#' Split a lower/upper bound set into pinned (fixed) params and free bounds
#'
#' A parameter whose lower and upper bounds are both present and equal is treated
#' as PINNED: it goes into `fixed` at that value (and `start` is set to it). The
#' remaining bounds pass through as true ranges. This routes "fix via Lower =
#' Upper" through [fit_model()]'s `fixed` argument, avoiding the `lower == upper`
#' error that L-BFGS-B would otherwise raise.
#' @param lower,upper Named numeric vectors of bounds (may be empty or NULL).
#' @param start Named numeric starting vector (all model params).
#' @return A list: `fixed` (character), `start` (with pinned values applied),
#'   `lower`, `upper` (named numerics with pinned params removed, or NULL).
#' @keywords internal
split_fixed_bounds <- function(lower, upper, start) {
  if (is.null(lower)) lower <- numeric(0)
  if (is.null(upper)) upper <- numeric(0)
  common <- intersect(names(lower), names(upper))
  eq <- common[is.finite(lower[common]) & is.finite(upper[common]) &
                 abs(lower[common] - upper[common]) <= 1e-12]
  start[eq] <- lower[eq]
  drop_lo <- lower[setdiff(names(lower), eq)]
  drop_hi <- upper[setdiff(names(upper), eq)]
  list(fixed = eq,
       start = start,
       lower = if (length(drop_lo)) drop_lo else NULL,
       upper = if (length(drop_hi)) drop_hi else NULL)
}
