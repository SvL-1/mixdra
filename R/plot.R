# Thin plotly renderers. Each calls a data builder from R/plot-data.R and
# assembles an interactive plotly object. `plotly` is a Suggested dependency, so
# every renderer first checks it is installed.

#' Error if plotly is not available
#' @keywords internal
require_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("The 'plotly' package is required for plotting. ",
         "Install it with install.packages(\"plotly\").", call. = FALSE)
}

#' Plot a chemical's dose-response curve
#'
#' Draws the fitted log-logistic marginal curve with the observed points
#' overlaid, as an interactive plotly object.
#' @param fit An enriched fit from [fit_model()] or [fit_single()].
#' @param df The data frame the fit was built from.
#' @param chem Index of the chemical whose marginal to draw (mixtures only).
#' @param log_x Use a log10 concentration axis (default `TRUE`). Control points
#'   at concentration 0 are not shown on a log axis.
#' @return A plotly object.
#' @export
plot_dose_response <- function(fit, df, chem = 1, log_x = TRUE) {
  require_plotly()
  d <- dr_curve_data(fit, df, chem)
  p <- plotly::plot_ly()
  p <- plotly::add_markers(p, x = d$observed$conc, y = d$observed$response,
                           name = "observed", marker = list(color = "black", size = 6))
  p <- plotly::add_lines(p, x = d$curve$conc, y = d$curve$response,
                         name = "fitted", line = list(color = "steelblue"))
  plotly::layout(p,
    xaxis = list(title = d$chem, type = if (log_x) "log" else "linear"),
    yaxis = list(title = "response"))
}
