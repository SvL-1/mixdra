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

#' Plot observed vs predicted response
#'
#' Scatter of observed against predicted values with a 1:1 reference line, as a
#' plotly object. Works for any fit.
#' @param fit An enriched fit from [fit_model()] or [fit_single()].
#' @param df The data frame the fit was built from.
#' @return A plotly object.
#' @export
plot_obs_pred <- function(fit, df) {
  require_plotly()
  d <- obs_pred_data(fit, df)
  lim <- range(c(d$observed, d$predicted), na.rm = TRUE)
  p <- plotly::plot_ly()
  p <- plotly::add_markers(p, x = d$observed, y = d$predicted, name = "points",
                           marker = list(color = "black", size = 6))
  p <- plotly::add_lines(p, x = lim, y = lim, name = "1:1",
                         line = list(color = "grey", dash = "dash"))
  plotly::layout(p, xaxis = list(title = "observed"),
                 yaxis = list(title = "predicted"))
}

#' Plot the fitted 3-D response surface of a binary mixture
#'
#' Observed points (`scatter3d`) overlaid on the fitted response surface, as an
#' interactive plotly object. Binary fits only.
#' @param fit An enriched binary fit from [fit_model()].
#' @param df The data frame the fit was built from.
#' @param n Grid resolution per axis (default 100).
#' @return A plotly object.
#' @export
plot_surface <- function(fit, df, n = 100) {
  require_plotly()
  g <- surface_grid_data(fit, df, n = n)
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, x = g$observed$x, y = g$observed$y, z = g$observed$z,
                         type = "scatter3d", mode = "markers",
                         marker = list(size = 3, color = "blue"), name = "observed")
  p <- plotly::add_surface(p, x = g$x_vals, y = g$y_vals, z = g$z,
                           opacity = 0.8, showscale = FALSE)
  plotly::layout(p, scene = list(
    xaxis = list(title = g$labels$x),
    yaxis = list(title = g$labels$y),
    zaxis = list(title = "response")))
}

#' Plot 2-D isoboles (equal-response contours) of a binary mixture
#'
#' Draws contour lines of equal response in the (C1, C2) plane at the requested
#' effect levels (solid, black). If `reference_fit` is supplied, its isoboles are
#' overlaid (dashed, red) so departures from additivity are visible. Binary fits
#' only.
#' @param fit An enriched binary fit from [fit_model()].
#' @param df The data frame the fit was built from.
#' @param levels Effect levels as fractions of `max` (default
#'   `c(0.1, 0.25, 0.5, 0.75, 0.9)`).
#' @param reference_fit Optional enriched reference fit to overlay (dashed).
#' @param n Grid resolution per axis.
#' @return A plotly object.
#' @export
plot_isobole <- function(fit, df, levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
                         reference_fit = NULL, n = 100) {
  require_plotly()
  d <- isobole_data(fit, df, levels = levels, reference_fit = reference_fit, n = n)
  labs <- attr(d, "labels")
  p <- plotly::plot_ly()
  for (grp in unique(d$group)) {
    seg <- d[d$group == grp, ]
    is_ref <- seg$source[1] == "reference"
    p <- plotly::add_lines(p, x = seg$x, y = seg$y, showlegend = FALSE,
                           name = paste0(seg$source[1], " ", signif(seg$level[1], 3)),
                           line = list(color = if (is_ref) "red" else "black",
                                       dash  = if (is_ref) "dash" else "solid"))
  }
  plotly::layout(p, xaxis = list(title = labs$x), yaxis = list(title = labs$y))
}
