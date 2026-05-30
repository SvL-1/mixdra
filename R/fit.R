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
#' @param lower,upper Optional named numeric vectors of manual box constraints,
#'   one entry per parameter to constrain (others keep their default). Defaults
#'   encode only domain meaning, not the seed: base curve parameters (`max`,
#'   `slope*`, `ec50*`) get positivity (`lower = 1e-8`, `upper = Inf`); for
#'   binary data `max` additionally gets `upper = 1` (it is a probability);
#'   deviation parameters (`a`, `b`, `b1`, `b2`, `b3`) are left fully
#'   unconstrained so downstream interaction analysis (e.g. the dose at which an
#'   interaction switches between synergism and antagonism) is undistorted.
#' @param n_starts Number of optimisation starts. The first uses `start`; each
#'   additional start perturbs `start` to escape the local minima of the
#'   non-smooth CA bisection surface. The best finite objective is kept.
#' @return A list with `par` (named fitted parameters), `objective`
#'   (residual SS or deviance), `pred` (fitted values), `residuals`, `df`
#'   (number of free parameters), `n`, and `convergence`.
#' @export
fit_model <- function(df, reference, deviation = "reference",
                      response = c("continuous", "binary"),
                      start, fixed = character(0),
                      lower = NULL, upper = NULL, n_starts = 1) {
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
  base <- setdiff(free, spec$extra)   # curve params; deviation params unconstrained

  # Default bounds encode domain meaning only (not the seed): base curve params
  # get positivity; deviation params stay unconstrained. The user may override
  # any individual bound via `lower`/`upper`.
  lo <- stats::setNames(rep(-Inf, length(free)), free)
  up <- stats::setNames(rep(Inf, length(free)), free)
  lo[base] <- 1e-8
  # Binary `max` is a control probability and cannot exceed 1.
  if (response == "binary" && "max" %in% base) up[["max"]] <- 1
  if (!is.null(lower)) lo[names(lower)] <- lower
  if (!is.null(upper)) up[names(upper)] <- upper

  # One attempt: L-BFGS-B (parscale normalises the very differently-scaled
  # parameters), with a Nelder-Mead fallback if it fails to converge — robust to
  # the non-smooth CA bisection inner solver.
  run_from <- function(theta_init) {
    # parscale sets the finite-difference step per parameter. Deviation params
    # (a, b, ...) start at 0; with a 1e-8 floor their FD step is far below the CA
    # bisection solver's noise floor, so L-BFGS-B sees a zero gradient and never
    # moves them off 0. Floor them at 1 so their effect is actually detected.
    ps <- pmax(abs(theta_init), 1e-8)
    if (length(spec$extra)) ps[spec$extra] <- pmax(ps[spec$extra], 1)
    res <- tryCatch(
      stats::optim(theta_init, obj_free, method = "L-BFGS-B",
                   lower = lo[free], upper = up[free],
                   control = list(parscale = ps, factr = 1e-9, maxit = 100)),
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
      # Base (positive curve) params: multiplicative jitter, clamped to bounds.
      theta_i[base] <- theta0[base] * exp(stats::runif(length(base), -log(3), log(3)))
      theta_i[base] <- pmin(pmax(theta_i[base], lo[base]), up[base])
      # Deviation params start at 0, so multiplicative jitter leaves them at 0;
      # perturb additively over a broad symmetric range so the optimiser explores
      # interaction (a, b, ...) away from the reference model.
      if (length(spec$extra))
        theta_i[spec$extra] <- theta0[spec$extra] +
          stats::runif(length(spec$extra), -10, 10)
    }
    res <- run_from(theta_i)
    if (!is.null(res) && is.finite(res$value) &&
        (is.null(best) || res$value < best$value)) best <- res
  }

  # Clamp the optimum into the box: the Nelder-Mead fallback is unbounded, so
  # this guarantees user-supplied (and default) constraints are always honoured.
  par[free] <- pmin(pmax(best$par, lo[free]), up[free])
  pred <- predict_with(par)
  obs <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  list(par = par, objective = objective_of(par), pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = best$convergence)
}
