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
  res <- stats::optim(start, obj, method = "L-BFGS-B",
                      lower = lower, upper = upper)
  list(par = res$par, ssr = res$value, convergence = res$convergence)
}
