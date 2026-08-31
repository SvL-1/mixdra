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

# --- production dispatch adapter ----------------------------------------------
# model_spec() returns make_adapter(...) as `$fn` for the reference/SA/DR/DL
# family. It preserves the historical by-name call contract used by
# `.mixture_eval()` / `param_ci()` -- do.call(fn, c(conc_list, as.list(par))) --
# and repacks each measurement row into mix_response()'s vector signature.

# Map the flat `b` parameter onto mix_response()'s expected shape.
#   DR binary : scalar b -> c(b, 0)   (b weights chemical 1's dose fraction)
#   DR ternary: b1,b2,b3 -> c(b1,b2,b3)   (not reached via model_spec today, but
#               kept faithful to mix_response()'s contract / the equivalence test)
#   DL        : scalar b passes through
#   reference / SA: b unused -> 0
.adapter_bmap <- function(deviation, n_chem, A) {
  if (deviation == "DR") {
    if (n_chem == 2) c(A$b, 0) else c(A$b1, A$b2, A$b3)
  } else if (deviation == "DL") {
    A$b
  } else {
    0
  }
}

#' Build a vectorised, by-name predictor over `mix_response()`
#'
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param n_chem 2 or 3.
#' @return A function accepting named concentration vectors `c1` and `c2` (and `c3` for ternary)
#'   and named scalar parameters (`max`, `slope1..`, `ec50..`, `a`, `b`/`b1..b3`),
#'   returning one prediction per row.
#' @keywords internal
make_adapter <- function(reference, deviation, n_chem) {
  slope_names <- if (n_chem == 2) c("slope1", "slope2")
                 else            c("slope1", "slope2", "slope3")
  ec50_names  <- if (n_chem == 2) c("ec501", "ec502")
                 else            c("ec50_1", "ec50_2", "ec50_3")
  cc_names    <- paste0("c", seq_len(n_chem))

  function(...) {
    A      <- list(...)
    concs  <- lapply(cc_names, function(k) A[[k]])
    slopes <- vapply(slope_names, function(k) A[[k]], numeric(1))
    ec50s  <- vapply(ec50_names,  function(k) A[[k]], numeric(1))
    a      <- if (is.null(A$a)) 0 else A$a
    b      <- .adapter_bmap(deviation, n_chem, A)
    n_pts  <- length(concs[[1]])
    vapply(seq_len(n_pts), function(i) {
      mix_response(vapply(concs, function(v) v[[i]], numeric(1)),
                   A$max, slopes, ec50s, reference, deviation, a = a, b = b)
    }, numeric(1))
  }
}
