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
#' Minimises the residual sum of squares.
#' @param conc Numeric vector of concentrations.
#' @param resp Numeric vector of responses (same length as `conc`).
#' @return A list with `par` (named: max, slope, ec50), `ssr`, and `convergence`.
#' @export
fit_single <- function(conc, resp) {
  stopifnot(length(conc) == length(resp), length(conc) > 3)
  pos <- conc[conc > 0]
  start <- c(max = max(resp, na.rm = TRUE),
             slope = 1,
             ec50 = stats::median(pos))
  obj <- function(p) {
    pred <- ll3_predict(conc, p[["max"]], p[["slope"]], p[["ec50"]])
    sum((resp - pred)^2)
  }
  lower <- c(max = 1e-8, slope = 1e-3, ec50 = 1e-8)
  upper <- c(max = Inf,  slope = 50,   ec50 = Inf)
  # parscale normalises the three very differently-scaled parameters (max ~ 1e2,
  # slope ~ 1, ec50 ~ 1e-1) so L-BFGS-B's relative tolerance bites uniformly;
  # without it the optimiser stops well short of the true minimum.
  res <- stats::optim(start, obj, method = "L-BFGS-B",
                      lower = lower, upper = upper,
                      control = list(parscale = pmax(abs(start), 1e-8),
                                     factr = 1e-9, maxit = 1000))
  list(par = res$par, ssr = res$value, convergence = res$convergence)
}
