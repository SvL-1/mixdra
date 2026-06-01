# Model-driven EC50 isoplane and Sigma-TU / z computation for the ternary
# Advanced-S/A fit. Pure (no plotly); these are scientific outputs in their own
# right, consumed by the plot builders in R/plot-data-ternary.R and reusable by
# the app. At the EC50 surface the implicit CA equation collapses to a closed
# form (see plan/spec): on the isoplane, sigma_tu = F4(z) and C_i = EC50_i*z_i*F4.

#' Triangular grid over the 2-simplex
#'
#' All `(i, j, k)` with `i + j + k = n` and `i, j, k >= 0`, scaled to
#' proportions summing to 1. Vertices (e.g. `(1,0,0)`) and edges (one coordinate
#' 0) are included.
#' @param n Grid resolution (points per simplex edge).
#' @return A data frame with columns `z1`, `z2`, `z3`.
#' @keywords internal
simplex_grid <- function(n) {
  rows <- list()
  for (i in 0:n) for (j in 0:(n - i)) rows[[length(rows) + 1L]] <- c(i, j, n - i - j)
  m <- do.call(rbind, rows) / n
  data.frame(z1 = m[, 1], z2 = m[, 2], z3 = m[, 3])
}

#' F4 deviation factor for a ternary Advanced-S/A model
#' @keywords internal
.f4 <- function(z1, z2, z3, A1, A2, A3, A4) {
  exp(A1 * z1 * z2 + A2 * z1 * z3 + A3 * z2 * z3 + A4 * z1 * z2 * z3)
}

#' EC50 isoplane points for a ternary Advanced-S/A fit
#'
#' Returns the EC50 surface in concentration space over a triangular grid of
#' TU-directions `z` on the 2-simplex. On the EC50 surface `sigma_tu = F4(z)` and
#' `C_i = EC50_i * z_i * F4(z)`.
#' @param res An [analyse_ternary()] result.
#' @param model `"SA"` (pairwise only, A4 = 0) or `"ASA"` (adds the overall
#'   three-way `res$A4_overall`).
#' @param n Grid resolution per simplex edge (default 30).
#' @return A data frame: `z1, z2, z3, C1, C2, C3, sigma_tu, model`.
#' @export
ec50_isoplane <- function(res, model = c("SA", "ASA"), n = 30) {
  model <- match.arg(model)
  b <- res$base
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]
  A4 <- if (model == "ASA") res$A4_overall else 0
  g <- simplex_grid(n)
  F4 <- .f4(g$z1, g$z2, g$z3, A1, A2, A3, A4)
  data.frame(z1 = g$z1, z2 = g$z2, z3 = g$z3,
             C1 = b[["ec50_1"]] * g$z1 * F4,
             C2 = b[["ec50_2"]] * g$z2 * F4,
             C3 = b[["ec50_3"]] * g$z3 * F4,
             sigma_tu = F4, model = model, stringsAsFactors = FALSE)
}

#' Sigma-TU vs z curves for a ternary Advanced-S/A fit
#'
#' For each chemical *x*, varies its TU-fraction `z_x` from 0 to 1 while holding
#' the other two chemicals at equal z (`(1 - z_x)/2` each), and reports `ΣTU =
#' F4(z)` along that path. This is the curve underlying the z-plot: deviation of
#' `sigma_tu` from 1 is the interaction (`< 1` synergy, `> 1` antagonism).
#' @param res An [analyse_ternary()] result.
#' @param model `"SA"` (pairwise only) or `"ASA"` (adds `res$A4_overall`).
#' @param n Number of z points per chemical (default 21, i.e. steps of 0.05).
#' @return A data frame: `chem` ("C1"/"C2"/"C3"), `z`, `sigma_tu`, `model`.
#' @export
sigma_tu_curve <- function(res, model = c("SA", "ASA"), n = 21) {
  model <- match.arg(model)
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]
  A4 <- if (model == "ASA") res$A4_overall else 0
  zx <- seq(0, 1, length.out = n)
  other <- (1 - zx) / 2
  chems <- c("C1", "C2", "C3")
  do.call(rbind, lapply(seq_along(chems), function(ix) {
    z <- matrix(other, nrow = n, ncol = 3)
    z[, ix] <- zx
    data.frame(chem = chems[ix], z = zx,
               sigma_tu = .f4(z[, 1], z[, 2], z[, 3], A1, A2, A3, A4),
               model = model, stringsAsFactors = FALSE)
  }))
}

#' EC50 marker points for each tested ternary ratio
#'
#' For every ratio in `res$individual`, converts its concentration proportions to
#' TU-fractions `z` and evaluates the EC50 isoplane point using that ratio's own
#' individual `A4`. Plotting these against the overall isoplane/curve is the
#' visual form of the per-ratio-vs-overall A4 "averaging-out" comparison.
#' @param res An [analyse_ternary()] result.
#' @param df Unused (kept for signature symmetry with the renderers); the marker
#'   geometry comes entirely from `res$individual` and `res$base`.
#' @return A data frame: `ratio, z1, z2, z3, C1, C2, C3, sigma_tu`.
#' @export
ec50_markers <- function(res, df = NULL) {
  b <- res$base
  ec <- c(b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]])
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]
  ind <- res$individual
  if (nrow(ind) == 0)
    return(data.frame(ratio = character(0), z1 = numeric(0), z2 = numeric(0),
                      z3 = numeric(0), C1 = numeric(0), C2 = numeric(0),
                      C3 = numeric(0), sigma_tu = numeric(0),
                      stringsAsFactors = FALSE))
  do.call(rbind, lapply(seq_len(nrow(ind)), function(i) {
    p <- c(ind$C1[i], ind$C2[i], ind$C3[i])    # concentration proportions
    tu <- p / ec
    z <- tu / sum(tu)
    F4 <- .f4(z[1], z[2], z[3], A1, A2, A3, ind$A4[i])
    data.frame(ratio = ind$ratio[i], z1 = z[1], z2 = z[2], z3 = z[3],
               C1 = ec[1] * z[1] * F4, C2 = ec[2] * z[2] * F4,
               C3 = ec[3] * z[3] * F4, sigma_tu = F4, stringsAsFactors = FALSE)
  }))
}
