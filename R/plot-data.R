# Pure data builders for the plotting layer. No plotly here, so these are fully
# unit-testable. The renderers in R/plot.R consume their output.

#' Observed response vector from a mixture/single data frame
#'
#' Continuous data carry `Res`; binary (quantal) data carry `Affected`/`Exposed`
#' and the modelled response is the proportion `Affected / Exposed`.
#' @keywords internal
obs_response <- function(df) {
  if ("Res" %in% names(df)) df$Res else df$Affected / df$Exposed
}

#' Build the data for one chemical's (marginal) dose-response curve
#'
#' When every other chemical is at 0, both CA and IA reduce to the single
#' three-parameter log-logistic, so the marginal curve is just [ll3_predict()]
#' with that chemical's `slope`/`ec50` and the shared `max`.
#' @param fit An enriched fit (from [fit_model()] or [fit_single()]).
#' @param df The data frame the fit was built from.
#' @param chem Index of the chemical whose marginal to draw (mixtures only).
#' @return A list: `curve` (data frame conc/response), `observed` (data frame
#'   conc/response), `chem` (the concentration column name).
#' @keywords internal
dr_curve_data <- function(fit, df, chem = 1) {
  if (isTRUE(fit$kind == "single")) {
    conc_col <- "C1"
    mx <- fit$par[["max"]]; sl <- fit$par[["slope"]]; ec <- fit$par[["ec50"]]
    keep <- rep(TRUE, nrow(df))
  } else {
    conc_cols <- fit$conc_cols
    if (chem < 1 || chem > length(conc_cols))
      stop("`chem` must be between 1 and ", length(conc_cols), call. = FALSE)
    conc_col <- conc_cols[chem]
    others <- setdiff(conc_cols, conc_col)
    keep <- if (length(others) == 0) rep(TRUE, nrow(df))
            else rowSums(df[others] == 0) == length(others)
    mx <- fit$par[["max"]]
    sl <- fit$par[[paste0("slope", chem)]]
    ec <- fit$par[[paste0("ec50", chem)]]
  }
  conc <- df[[conc_col]][keep]
  resp <- obs_response(df)[keep]
  pos  <- conc[conc > 0]
  grid <- seq(min(pos), max(conc), length.out = 200)
  list(curve    = data.frame(conc = grid, response = ll3_predict(grid, mx, sl, ec)),
       observed = data.frame(conc = conc, response = resp),
       chem     = conc_col)
}
