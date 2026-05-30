#' Sum-of-squares objective for continuous responses
#' @param obs,pred Numeric vectors of observed and predicted responses.
#' @return The residual sum of squares.
#' @export
obj_ss <- function(obs, pred) {
  sum((obs - pred)^2)
}

#' Per-row binomial log-likelihood contribution (Jonker Eq. 16 / VBA BinLik)
#'
#' Returns `P*log(pi_hat/pi) + (T-P)*log((1-pi_hat)/(1-pi))`, where `pi = P/T`.
#' Saturated terms (pi == 0 or pi == 1) drop the undefined component.
#' @param exposed Number exposed (T).
#' @param affected Number responding (P).
#' @param pi_hat Model-predicted probability in (0, 1).
#' @return A scalar (<= 0) log-likelihood contribution.
#' @export
binlik <- function(exposed, affected, pi_hat) {
  pi <- affected / exposed
  pos <- if (pi > 0) affected * log(pi_hat / pi) else 0
  neg <- if (pi < 1) (exposed - affected) * log((1 - pi_hat) / (1 - pi)) else 0
  pos + neg
}

#' Binomial deviance objective for binary responses
#'
#' `-2 * sum(binlik)`; minimising this maximises the binomial likelihood.
#' @param exposed,affected Numeric vectors of counts.
#' @param pi_hat Numeric vector of predicted probabilities.
#' @return The residual deviance (scalar).
#' @export
obj_deviance <- function(exposed, affected, pi_hat) {
  -2 * sum(mapply(binlik, exposed, affected, pi_hat))
}
