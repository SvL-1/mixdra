#' Fit one mixture model to a dataset
#'
#' @param df Data frame with `C1`, `C2` (and `C3` for ternary) plus either
#'   `Res` (continuous) or `Exposed`/`Affected` (binary).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param start Named numeric vector of starting values for the base parameters
#'   (max, slope*, ec50*). Deviation parameters default to 0.
#' @param fixed Character vector of parameter names to hold fixed at `start`.
#' @param lower_frac,upper_mult Box-constraint multipliers applied to positive
#'   base parameters (deviation params are unconstrained).
#' @param n_starts Number of optimisation starts. The first uses `start`; each
#'   additional start log-uniformly perturbs `start` (clamped to the bounds) to
#'   escape the local minima of the non-smooth CA bisection surface. The best
#'   finite objective is kept.
#' @return A list with `par` (named fitted parameters), `objective`
#'   (residual SS or deviance), `pred` (fitted values), `residuals`, `df`
#'   (number of free parameters), `n`, and `convergence`.
#' @export
fit_model <- function(df, reference, deviation = "reference",
                      response = c("continuous", "binary"),
                      start, fixed = character(0),
                      lower_frac = 0.1, upper_mult = 10, n_starts = 1) {
  response <- match.arg(response)
  conc_cols <- intersect(c("C1", "C2", "C3"), names(df))
  n_chem <- length(conc_cols)
  spec <- model_spec(reference, deviation, n_chem)

  # Assemble the full starting vector. Zero-initialising means any deviation
  # parameter (a, b, b1, b2, b3) not supplied in `start` defaults to 0.
  par <- stats::setNames(numeric(length(spec$params)), spec$params)
  par[names(start)] <- start[names(start)]

  free <- setdiff(spec$params, fixed)
  conc <- stats::setNames(as.list(df[conc_cols]), paste0("c", seq_along(conc_cols)))

  predict_with <- function(p_full) {
    args <- c(conc, as.list(p_full))
    names(args) <- c(names(conc), names(p_full))
    do.call(spec$fn, args)
  }
  objective_of <- function(p_full) {
    pred <- predict_with(p_full)
    if (response == "continuous") {
      obj_ss(df$Res, pred)
    } else {
      obj_deviance(df$Exposed, df$Affected, pmin(pmax(pred, 1e-8), 1 - 1e-8))
    }
  }

  obj_free <- function(theta) {
    p_full <- par
    p_full[free] <- theta
    val <- objective_of(p_full)
    if (!is.finite(val)) 1e12 else val
  }

  theta0 <- par[free]
  base <- setdiff(free, spec$extra)   # bounded curve params; deviation params unconstrained
  lower <- stats::setNames(rep(-Inf, length(free)), free)
  upper <- stats::setNames(rep(Inf, length(free)), free)
  lower[base] <- pmax(1e-8, theta0[base] * lower_frac)
  upper[base] <- pmax(lower[base] * 1.01, theta0[base] * upper_mult)

  # One attempt: L-BFGS-B (parscale normalises the very differently-scaled
  # parameters), with a Nelder-Mead fallback if it fails to converge — robust to
  # the non-smooth CA bisection inner solver.
  run_from <- function(theta_init) {
    ps <- pmax(abs(theta_init), 1e-8)
    res <- tryCatch(
      stats::optim(theta_init, obj_free, method = "L-BFGS-B",
                   lower = lower[free], upper = upper[free],
                   control = list(parscale = ps, factr = 1e-9, maxit = 1000)),
      error = function(e) NULL)
    if (is.null(res) || res$convergence != 0) {
      res2 <- tryCatch(
        stats::optim(theta_init, obj_free, method = "Nelder-Mead",
                     control = list(parscale = ps, maxit = 2000)),
        error = function(e) NULL)
      if (!is.null(res2) && (is.null(res) || res2$value < res$value)) res <- res2
    }
    res
  }

  best <- NULL
  for (i in seq_len(n_starts)) {
    theta_i <- theta0
    if (i > 1) {
      theta_i <- theta0 * exp(stats::runif(length(theta0), -log(3), log(3)))
      theta_i[base] <- pmin(pmax(theta_i[base], lower[base]), upper[base])
    }
    res <- run_from(theta_i)
    if (!is.null(res) && is.finite(res$value) &&
        (is.null(best) || res$value < best$value)) best <- res
  }

  par[free] <- best$par
  pred <- predict_with(par)
  obs <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  list(par = par, objective = best$value, pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = best$convergence)
}
