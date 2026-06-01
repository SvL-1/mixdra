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

#' Plot the EC50 isoplane of a ternary mixture
#'
#' 3-D scatter of the CA+S/A and CA+S/A+S/A EC50 isoplane point-clouds with the
#' per-ratio EC50 markers overlaid, as an interactive plotly object. Mirrors the
#' MixTox isoplane figures, sourced from the fitted model.
#' @param res An [analyse_ternary()] result.
#' @param df The data frame the fit was built from (passed to [ec50_markers()]).
#' @param n Isoplane grid resolution per simplex edge (default 30).
#' @return A plotly object.
#' @export
plot_isoplane <- function(res, df = NULL, n = 30) {
  require_plotly()
  d <- isoplane_plot_data(res, df, n = n)
  cols <- c("CA+S/A" = "orange", "CA+S/A+S/A" = "blue", "EC50" = "red")
  sizes <- c("CA+S/A" = 3, "CA+S/A+S/A" = 3, "EC50" = 7)
  p <- plotly::plot_ly()
  for (s in c("CA+S/A", "CA+S/A+S/A", "EC50")) {
    seg <- d[d$series == s, ]
    p <- plotly::add_trace(p, x = seg$C1, y = seg$C2, z = seg$C3,
                           type = "scatter3d", mode = "markers", name = s,
                           marker = list(size = sizes[[s]], color = cols[[s]]))
  }
  plotly::layout(p, scene = list(xaxis = list(title = "C1"),
                                 yaxis = list(title = "C2"),
                                 zaxis = list(title = "C3")))
}

#' Plot Sigma-TU vs z for a ternary mixture
#'
#' Per-chemical ΣTU-vs-z curves under CA+S/A (solid) and CA+S/A+S/A (dashed),
#' with the additivity reference line at ΣTU = 1, as an interactive plotly
#' object. Deviation from 1 is the interaction (< 1 synergy, > 1 antagonism).
#' @param res An [analyse_ternary()] result.
#' @param n Number of z points per chemical (default 21).
#' @return A plotly object.
#' @export
plot_sigma_tu <- function(res, n = 21) {
  require_plotly()
  chem_cols <- c(C1 = "#E66100", C2 = "#FFC20A", C3 = "#5D3A9B")
  d <- sigma_tu_plot_data(res, n = n)
  p <- plotly::plot_ly()
  for (ch in c("C1", "C2", "C3")) {
    for (s in c("CA+S/A", "CA+S/A+S/A")) {
      seg <- d[d$chem == ch & d$series == s, ]
      p <- plotly::add_lines(p, x = seg$z, y = seg$sigma_tu,
                             name = paste(ch, s),
                             line = list(color = chem_cols[[ch]],
                                         dash = if (s == "CA+S/A") "solid" else "dash"))
    }
  }
  p <- plotly::add_lines(p, x = c(0, 1), y = c(1, 1), name = "ΣTU = 1",
                         line = list(color = "black", width = 1))
  plotly::layout(p, xaxis = list(title = "z (chemical TU fraction)"),
                 yaxis = list(title = "ΣTU"))
}
