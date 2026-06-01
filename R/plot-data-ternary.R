# Pure builders for the ternary plotting layer. No plotly here; they reshape the
# R/ternary-isoplane.R outputs into single long-format frames with a `series`
# label so the renderers in R/plot.R stay thin.

#' Combine isoplane clouds + EC50 markers into one plotting frame
#'
#' Stacks the CA+S/A isoplane, the CA+S/A+S/A isoplane, and the per-ratio EC50
#' markers, tagging each with a `series` label.
#' @param res An [analyse_ternary()] result.
#' @param df Passed through to [ec50_markers()] (unused there).
#' @param n Isoplane grid resolution.
#' @return A data frame: `C1, C2, C3, series`.
#' @keywords internal
isoplane_plot_data <- function(res, df = NULL, n = 30) {
  sa  <- ec50_isoplane(res, "SA", n)
  asa <- ec50_isoplane(res, "ASA", n)
  mk  <- ec50_markers(res, df)
  rbind(
    data.frame(C1 = sa$C1,  C2 = sa$C2,  C3 = sa$C3,  series = "CA+S/A",
               stringsAsFactors = FALSE),
    data.frame(C1 = asa$C1, C2 = asa$C2, C3 = asa$C3, series = "CA+S/A+S/A",
               stringsAsFactors = FALSE),
    data.frame(C1 = mk$C1,  C2 = mk$C2,  C3 = mk$C3,  series = "EC50",
               stringsAsFactors = FALSE))
}

#' Combine SA and ASA Sigma-TU curves into one plotting frame
#' @param res An [analyse_ternary()] result.
#' @param n Number of z points per chemical.
#' @return A data frame: `chem, z, sigma_tu, series`.
#' @keywords internal
sigma_tu_plot_data <- function(res, n = 21) {
  sa  <- sigma_tu_curve(res, "SA", n)
  asa <- sigma_tu_curve(res, "ASA", n)
  rbind(
    data.frame(chem = sa$chem,  z = sa$z,  sigma_tu = sa$sigma_tu,
               series = "CA+S/A", stringsAsFactors = FALSE),
    data.frame(chem = asa$chem, z = asa$z, sigma_tu = asa$sigma_tu,
               series = "CA+S/A+S/A", stringsAsFactors = FALSE))
}
