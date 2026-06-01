# Unified mixture predictor.
#
# Replaces the 16 verbatim-ported per-(reference x deviation x n_chem) functions
# in models-binary.R / models-ternary.R with a single n-agnostic predictor. The
# whole family differs in exactly two places:
#   (a) the interaction scalar `Fdev` (reference / SA / DR / DL), and
#   (b) the combination rule (CA bisection vs IA probit/product).
# Everything else -- degenerate inputs, toxic units, dose fractions, the
# increasing/decreasing slope-sign duality -- is shared.
#
# Equivalence to the originals is pinned by tests/testthat/test-mixture-predict.R.

#' Predicted response for a single mixture measurement point
#'
#' @param concs Numeric vector of the n chemical concentrations at this point.
#' @param max Control response (response at zero dose for positive slopes).
#' @param slopes,ec50s Per-chemical curve parameters (length n).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param a Scalar interaction parameter (SA/DR/DL).
#' @param b For "DR" a per-chemical numeric vector (length n); for "DL" a scalar.
#' @return Predicted response (scalar).
#' @keywords internal
mix_response <- function(concs, max, slopes, ec50s,
                         reference = c("CA", "IA"),
                         deviation = c("reference", "SA", "DR", "DL"),
                         a = 0, b = 0) {
  reference <- match.arg(reference)
  deviation <- match.arg(deviation)

  active <- which(concs > 0)
  pos_all <- all(slopes > 0)               # positive slopes => response decreases with dose

  # --- degenerate cases, shared by every model ---
  if (length(active) == 0) return(if (pos_all) max else 0)
  if (length(active) == 1) {
    i <- active
    return(max / (1 + (concs[i] / ec50s[i])^slopes[i]))
  }

  ca <- concs[active]; sl <- slopes[active]; ec <- ec50s[active]
  pos <- all(sl > 0)                       # sign regime over the active chemicals
  TU <- ca / ec
  z  <- TU / sum(TU)

  # (a) interaction scalar
  core <- switch(deviation,
    reference = 0,
    SA = a * prod(z),
    DR = (a + sum(b[active] * z)) * prod(z),
    DL = {
      if (reference == "CA") {
        a * (1 - b * sum(TU)) * prod(z)
      } else {
        f <- 1 / (1 + (ca / ec)^sl)
        p50 <- if (pos) 1 - prod(f) else 1 - prod(1 - f)
        a * (1 - b * p50) * prod(z)
      }
    })

  # (b) combination rule
  if (reference == "CA") {
    rhs <- if (deviation == "reference") 1 else exp(core)
    L <- 0; U <- max; Y <- (L + U) / 2
    for (i in seq_len(200)) {
      Y <- (L + U) / 2
      eci <- ec * ((max - Y) / Y)^(1 / sl)
      G <- sum(ca / eci) - rhs
      if (pos) { if (G < 0) L <- Y else U <- Y }
      else     { if (G > 0) L <- Y else U <- Y }
      if (U - L < 1e-10) break
    }
    Y
  } else {
    f <- 1 / (1 + (ca / ec)^sl)
    base <- if (pos) prod(f) else 1 - prod(1 - f)
    shift <- if (pos) core else -core
    max * stats::pnorm(stats::qnorm(base) + shift)
  }
}
