#' Three-parameter log-logistic dose-response prediction
#'
#' `Y = max / (1 + (C / ec50)^slope)`. With `slope > 0` the response
#' decreases with dose (control = `max`).
#' @param conc Numeric vector of concentrations (>= 0).
#' @param max,slope,ec50 Curve parameters.
#' @return Numeric vector of predicted responses.
#' @export
ll3_predict <- function(conc, max, slope, ec50) {
  max / (1 + (conc / ec50)^slope)
}

#' Fit a three-parameter log-logistic curve to single-chemical data
#'
#' Minimises the residual sum of squares. `lower`, `upper`, and `start` are
#' optional **named, partial** overrides of the built-in defaults
#' (`lower = max 1e-8 / slope 1e-3 / ec50 1e-8`, `upper = max Inf / slope 50 /
#' ec50 Inf`, and a data-driven start). The effective start is clamped into
#' `[lower, upper]` so the optimiser never receives an out-of-bounds seed.
#' @param conc Numeric vector of concentrations.
#' @param resp Numeric vector of responses (same length as `conc`).
#' @param lower,upper,start Optional named numeric vectors (`max`/`slope`/`ec50`)
#'   overriding the corresponding defaults.
#' @param fixed Character vector of parameter names to hold FIXED at their
#'   `start` value. Pinning has to go through here rather than through equal
#'   lower/upper bounds: L-BFGS-B cannot take a finite difference inside a
#'   zero-width box and fails with "non-finite finite-difference value".
#' @return A list with `par` (named: max, slope, ec50), `ssr`, `convergence`, and
#'   `kind = "single"`.
#' @export
fit_single <- function(conc, resp, lower = NULL, upper = NULL, start = NULL,
                       fixed = character(0)) {
  stopifnot(length(conc) == length(resp), length(conc) > 3)
  pos <- conc[conc > 0]
  def_start <- c(max = max(resp, na.rm = TRUE), slope = 1, ec50 = stats::median(pos))
  def_lower <- c(max = 1e-8, slope = 1e-3, ec50 = 1e-8)
  def_upper <- c(max = Inf,  slope = 50,   ec50 = Inf)

  lo <- def_lower; if (!is.null(lower)) lo[names(lower)] <- lower
  hi <- def_upper; if (!is.null(upper)) hi[names(upper)] <- upper
  st <- def_start; if (!is.null(start)) st[names(start)] <- start
  st <- pmin(pmax(st, lo), hi)   # keep the seed inside the box

  ssr_of <- function(p)
    sum((resp - ll3_predict(conc, p[["max"]], p[["slope"]], p[["ec50"]]))^2)

  free <- setdiff(names(st), intersect(fixed, names(st)))
  # Every parameter pinned: nothing to optimise, just evaluate the curve.
  if (!length(free))
    return(list(par = st, ssr = ssr_of(st), convergence = 0L, kind = "single"))

  obj <- function(p) { full <- st; full[free] <- p; ssr_of(full) }
  # parscale normalises the three very differently-scaled parameters (max ~ 1e2,
  # slope ~ 1, ec50 ~ 1e-1) so L-BFGS-B's relative tolerance bites uniformly;
  # without it the optimiser stops well short of the true minimum.
  res <- stats::optim(st[free], obj, method = "L-BFGS-B",
                      lower = lo[free], upper = hi[free],
                      control = list(parscale = pmax(abs(st[free]), 1e-8),
                                     factr = 1e-9, maxit = 1000))
  par <- st; par[free] <- res$par
  list(par = par, ssr = res$value, convergence = res$convergence,
       kind = "single")
}

#' Evaluate the log-logistic curve at given parameters (no optimisation)
#'
#' Forward evaluation used by the Shiny app's "Simulate" action: computes the
#' predicted response and residual sum of squares for caller-supplied parameters,
#' returning the same shape as [fit_single()] so the plotting layer consumes it
#' unchanged.
#' @param conc,resp Concentration and observed-response vectors (same length).
#' @param max,slope,ec50 Curve parameters to evaluate.
#' @return A list with `par` (named max/slope/ec50), `ssr`, `convergence = NA`,
#'   and `kind = "single"`.
#' @keywords internal
eval_single <- function(conc, resp, max, slope, ec50) {
  pred <- ll3_predict(conc, max, slope, ec50)
  list(par = c(max = max, slope = slope, ec50 = ec50),
       ssr = sum((resp - pred)^2),
       convergence = NA_integer_, kind = "single")
}
