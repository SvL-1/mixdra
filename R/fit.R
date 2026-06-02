#' Assemble the full named starting vector for a model fit
#'
#' Base curve parameters default to 0 (always supplied via `start`); interaction
#' parameters default to their solver origin: `a` at 0 and every `b`-family
#' parameter (`b`, `b1`, `b2`, `b3`) at 1. Values present in `start` override the
#' defaults. Any other deviation parameter (e.g. the Advanced S/A `A*` terms) is
#' left at 0.
#' @keywords internal
init_start_par <- function(params, extra, start) {
  par <- stats::setNames(numeric(length(params)), params)  # all 0
  bfam <- intersect(extra, c("b", "b1", "b2", "b3"))
  if (length(bfam)) par[bfam] <- 1
  par[names(start)] <- start[names(start)]
  par
}

#' Fit one mixture model to a dataset
#'
#' @param df Data frame with `C1`, `C2` (and `C3` for ternary) plus either
#'   `Res` (continuous) or `Exposed`/`Affected` (binary).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param start Named numeric vector of starting values for the base parameters
#'   (max, slope*, ec50*). Deviation parameters default to their solver origin:
#'   `a` at 0 and `b`/`b1`/`b2`/`b3` at 1; any other deviation parameter at 0.
#' @param fixed Character vector of parameter names to hold fixed at `start`.
#' @param lower,upper Optional named numeric vectors of hard bounds. Any model
#'   parameter may be named here, including interaction parameters (`a`, `b`,
#'   `b1/b2/b3`). Parameters not named fall back to defaults: base curve params
#'   (`max`, `slope*`, `ec50*`) default to positivity (`lower = 1e-8`,
#'   `upper = Inf`); binary `max` additionally defaults to `upper = 1` (a value
#'   above 1 is capped, with a warning, since it is a probability); interaction
#'   params default to unconstrained (`-Inf`/`Inf`) when not named. A start
#'   value outside its bounds is clamped into range with a warning.
#' @param n_starts Number of optimisation starts. The first uses `start`; each
#'   additional start log-uniformly perturbs `start` (clamped to the bounds) to
#'   escape the local minima of the non-smooth CA bisection surface. The best
#'   finite objective is kept.
#' @param time_limit Wall-clock budget in seconds for the whole fit (across all
#'   starts and the Nelder-Mead fallback). When exceeded, the running optimiser
#'   is aborted and the best parameters seen so far are returned with
#'   `convergence = 99` and a warning. `NULL` disables the limit. This is a
#'   debugging guard against pathological non-convergence, not a substitute for
#'   sensible starts/bounds.
#' @return A list with `par` (named fitted parameters), `objective`
#'   (residual SS or deviance), `pred` (fitted values), `residuals`, `df`
#'   (number of free parameters), `n`, `convergence`, and `fixed` (character
#'   vector of parameter names held fixed at their `start` values).
#' @export
fit_model <- function(df, reference, deviation = "reference",
                      response = c("continuous", "binary"),
                      start, fixed = character(0),
                      lower = NULL, upper = NULL, n_starts = 1,
                      time_limit = 30) {
  response <- match.arg(response)
  conc_cols <- intersect(c("C1", "C2", "C3"), names(df))
  n_chem <- length(conc_cols)
  spec <- model_spec(reference, deviation, n_chem)

  # Assemble the full starting vector: interaction params default to their solver
  # origin (a = 0, b-family = 1); curve params come from `start`.
  par <- init_start_par(spec$params, spec$extra, start)

  free <- setdiff(spec$params, fixed)
  conc <- stats::setNames(as.list(df[conc_cols]), paste0("c", seq_along(conc_cols)))

  predict_with <- function(p_full) .mixture_eval(spec$fn, conc, p_full)
  objective_of <- function(p_full) {
    pred <- predict_with(p_full)
    if (response == "continuous") {
      obj_ss(df$Res, pred)
    } else {
      obj_deviance(df$Exposed, df$Affected, pmin(pmax(pred, 1e-8), 1 - 1e-8))
    }
  }

  # Any model parameter may be bounded, including the interaction params (a, b),
  # which the final joint "Optimize all" stage constrains/pins. Validate names
  # here -- before the all-fixed fast path -- so an illegal bound name is rejected
  # consistently regardless of how many params are free (otherwise the staged
  # reference fit would silently ignore it).
  bad <- setdiff(c(names(lower), names(upper)), spec$params)
  if (length(bad))
    stop("`lower`/`upper` may only name a model parameter (",
         paste(spec$params, collapse = ", "), "); got: ",
         paste(unique(bad), collapse = ", "))

  # When every parameter is fixed (the staged reference fit: curve params held at
  # their single-compound values, no interaction params to free) there is nothing
  # to optimise -- just evaluate the model. `stats::optim` cannot run on a
  # zero-length parameter vector, so short-circuit here.
  if (length(free) == 0) {
    pred <- predict_with(par)
    obs <- if (response == "continuous") df$Res else df$Affected / df$Exposed
    return(list(par = par, objective = objective_of(par), pred = pred,
                residuals = obs - pred, df = 0L, n = nrow(df),
                convergence = 0L, reference = reference, deviation = deviation,
                response = response, conc_cols = conc_cols, n_chem = n_chem,
                fixed = fixed, kind = "mixture"))
  }

  theta0 <- par[free]
  base <- setdiff(free, spec$extra)             # positively-bounded curve params

  # Defaults: base params get positivity; deviation params stay +/-Inf unless
  # the caller bounds them.
  lo <- stats::setNames(rep(-Inf, length(free)), free)
  hi <- stats::setNames(rep(Inf, length(free)), free)
  lo[base] <- 1e-8
  # For binary data `max` is the control response *probability*, so it cannot
  # exceed 1; without this cap the optimiser can push it above 1 and the model
  # returns fitted probabilities > 1 (`ca_bi = max / (1 + ...) <= max`).
  bin_max <- response == "binary" && "max" %in% base
  if (bin_max) hi[["max"]] <- 1

  # Apply user-supplied bounds to any free parameter (incl. interaction a/b).
  for (p in intersect(names(lower), free)) lo[[p]] <- lower[[p]]
  for (p in intersect(names(upper), free)) hi[[p]] <- upper[[p]]
  if (bin_max && hi[["max"]] > 1) {
    warning("binary `max` upper bound capped at 1 (requested ", upper[["max"]], ")")
    hi[["max"]] <- 1
  }

  if (length(free) && any(lo[free] >= hi[free]))
    stop("each parameter's lower bound must be below its upper bound; check: ",
         paste(free[lo[free] >= hi[free]], collapse = ", "))

  # Keep the optimiser's starting point feasible.
  if (length(free)) {
    clamped <- pmin(pmax(theta0[free], lo[free]), hi[free])
    if (any(clamped != theta0[free])) {
      off <- free[clamped != theta0[free]]
      warning("start value(s) outside bounds, clamped: ",
              paste(off, collapse = ", "))
      theta0[free] <- clamped
    }
  }
  lower <- lo
  upper <- hi

  # Wall-clock budget. `optim` has no time limit, so we enforce one from inside
  # the objective: track the best feasible point seen and, once the deadline
  # passes, abort the running optimiser with an error (caught in `run_from`).
  # `best_seen` then provides a usable result even when nothing converged.
  t0 <- Sys.time()
  timed_out <- FALSE
  best_seen <- list(value = Inf, par = theta0)

  obj_free <- function(theta) {
    if (!is.null(time_limit) && !timed_out &&
        as.numeric(difftime(Sys.time(), t0, units = "secs")) > time_limit)
      timed_out <<- TRUE
    if (timed_out) stop("mixdra_time_limit")   # abort the active optim run
    # Enforce the box constraints for *every* optimiser. L-BFGS-B respects
    # `lower`/`upper` natively, but the Nelder-Mead fallback does not, so reject
    # infeasible points here to make the bounds bite regardless of method.
    if (any(theta < lower[free]) || any(theta > upper[free])) return(1e12)
    p_full <- par
    p_full[free] <- theta
    val <- objective_of(p_full)
    if (!is.finite(val)) return(1e12)
    if (val < best_seen$value) best_seen <<- list(value = val, par = theta)
    val
  }

  # One attempt: L-BFGS-B (parscale normalises the very differently-scaled
  # parameters), with a Nelder-Mead fallback if it fails to converge — robust to
  # the non-smooth CA bisection inner solver.
  run_from <- function(theta_init) {
    # parscale sets the finite-difference step per parameter. Deviation params
    # (a, b, ...) start at 0; with a 1e-8 floor their FD step is far below the CA
    # bisection solver's noise floor, so L-BFGS-B sees a zero gradient and never
    # moves them off 0. Floor them at 1 so their effect is actually detected.
    ps <- pmax(abs(theta_init), 1e-8)
    ex <- intersect(spec$extra, free)          # only the FREE deviation params
    if (length(ex)) ps[ex] <- pmax(ps[ex], 1)
    res <- tryCatch(
      stats::optim(theta_init, obj_free, method = "L-BFGS-B",
                   lower = lower[free], upper = upper[free],
                   control = list(parscale = ps, factr = 1e-9, maxit = 200)),
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
      theta_i[base] <- pmin(pmax(theta_i[base], lower[base]), upper[base])
      # Deviation params start at their origin (a = 0, b-family = 1); perturb them
      # additively over a broad symmetric range so the optimiser explores
      # interaction (a, b, ...) away from the reference model.
      ex <- intersect(spec$extra, free)        # only the FREE deviation params
      if (length(ex))
        theta_i[ex] <- theta0[ex] + stats::runif(length(ex), -10, 10)
    }
    res <- run_from(theta_i)
    if (!is.null(res) && is.finite(res$value) &&
        (is.null(best) || res$value < best$value)) best <- res
  }

  # If the deadline cut every start short before any optimiser returned cleanly,
  # fall back to the best feasible point evaluated so far (convergence = 99).
  if (is.null(best))
    best <- list(par = best_seen$par, value = best_seen$value, convergence = 99L)
  if (timed_out)
    warning(sprintf(
      "fit_model (%s/%s): time limit (%gs) reached; returning best-so-far%s.",
      reference, deviation, time_limit,
      if (identical(best$convergence, 99L)) " (convergence = 99)" else ""))

  # Clamp the optimum into the box (the Nelder-Mead fallback is unbounded), so
  # user-supplied and default constraints are always honoured in the result.
  par[free] <- pmin(pmax(best$par, lower[free]), upper[free])
  pred <- predict_with(par)
  obs <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  list(par = par, objective = objective_of(par), pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = best$convergence,
       reference = reference, deviation = deviation, response = response,
       conc_cols = conc_cols, n_chem = n_chem, fixed = fixed, kind = "mixture")
}

#' Forward-evaluate a mixture model at fixed parameters (Simulate)
#'
#' Mixture counterpart of [eval_single()]. Assembles the full parameter vector
#' from the frozen curve parameters plus caller-supplied interaction parameters
#' and evaluates the model without optimising, returning the enriched shape of
#' [fit_model()] (so the plotting layer consumes it unchanged) with an added
#' `simulated = TRUE` flag. Delegates to [fit_model()]'s all-fixed path.
#'
#' @param df Mixture data frame (`C1`/`C2`[/`C3`] + response columns).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param curve_params Named numeric of curve params (max, slope1, slope2,
#'   ec501, ec502).
#' @param interaction Named numeric of interaction params (e.g. `c(a = .., b = ..)`).
#'   Entries the chosen model does not use are dropped.
#' @return An enriched mixture fit list (as [fit_model()]) with `simulated = TRUE`.
#' @keywords internal
eval_mixture <- function(df, reference, deviation, response,
                         curve_params, interaction = numeric(0)) {
  n_chem <- length(intersect(c("C1", "C2", "C3"), names(df)))
  spec <- model_spec(reference, deviation, n_chem)
  use  <- interaction[intersect(names(interaction), spec$extra)]
  par  <- c(curve_params, use)
  fit  <- fit_model(df, reference, deviation, response,
                    start = par, fixed = spec$params)   # all fixed -> evaluate only
  fit$simulated <- TRUE
  fit
}

#' Joint refinement of a mixture fit ("Optimize all params")
#'
#' Re-fits EVERY parameter of an existing mixture fit at once, seeded from that
#' fit's parameters, optionally constrained by `lower`/`upper`. A parameter whose
#' lower and upper bounds are equal is pinned (held fixed) via
#' [split_fixed_bounds()]. Because the seed is the prior fit and the bisection
#' surface is non-smooth, multi-start (`n_starts`) is recommended. The returned
#' fit carries `joint = TRUE`.
#' @param fit An enriched mixture fit (from [fit_model()] / the staged analysis).
#' @param df The mixture data frame the fit was built from.
#' @param lower,upper Optional named numeric vectors of bounds over any model
#'   parameter (curve or interaction). Equal lower==upper pins that parameter.
#' @param n_starts,time_limit Forwarded to [fit_model()].
#' @return An enriched fit (as [fit_model()]) with `joint = TRUE`.
#' @keywords internal
refine_joint <- function(fit, df, lower = NULL, upper = NULL,
                         n_starts = 1, time_limit = 30) {
  spec  <- model_spec(fit$reference, fit$deviation, fit$n_chem)
  start <- fit$par[spec$params]
  sp    <- split_fixed_bounds(lower, upper, start)
  out <- fit_model(df, fit$reference, fit$deviation, fit$response,
                   start = sp$start, fixed = sp$fixed,
                   lower = sp$lower, upper = sp$upper,
                   n_starts = n_starts, time_limit = time_limit)
  out$joint <- TRUE
  out
}
