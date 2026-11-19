#' (De)potentiation: does a second stressor shift stressor 1's EC50?
#'
#' Fits stressor 1's three-parameter log-logistic curve ([ll3_predict()]) to
#' every series: the reference series (stressor 1 alone) and each series with
#' stressor 2 applied at some level. Each non-reference series is then tested
#' against the reference by a likelihood-ratio test with one degree of freedom:
#' run 1 fits both curves separately (`obj1`), run 2 fits them with one shared
#' EC50 while max and slope stay per curve (`obj2`). Continuous data use the
#' residual sum of squares and `X2 = n * log(obj2 / obj1)`; quantal data use the
#' binomial deviance and `X2 = obj2 - obj1`. `X2` is floored at 0 (run 2 is
#' nested in run 1, so a negative value is optimiser noise).
#'
#' Controls: a series without `C1 == 0` rows borrows the reference series'
#' controls, which then count in both curves of that comparison (as in the
#' Excel workbook this method comes from). A series with its own controls uses
#' those.
#'
#' A significant shift with `EC50_series / EC50_reference < 1` is reported as
#' `"more toxic"`, above 1 as `"less toxic"`; otherwise
#' `"no significant influence"`.
#' @param df Data frame with columns `Series` (curve label), `C1` (stressor-1
#'   concentration, >= 0) and `Res` (continuous) or `Exposed`/`Affected`
#'   (quantal). Labels are compared after `trimws(as.character())`.
#' @param reference Label of the reference series; defaults to the first series
#'   in `df`.
#' @param response "continuous" or "quantal".
#' @param alpha Significance level of the likelihood-ratio test.
#' @param stressor1 Display name of stressor 1, used in the summary text.
#' @return A list with `reference`, `response`, `alpha`, `series` (labels in
#'   file order); `curves` (run-1 fit per series: `par` = named
#'   `c(max, slope, ec50)`, `objective`, `n`, `convergence`,
#'   `controls_borrowed`, `ec50_extrapolated`, and `data` = the rows used,
#'   including borrowed controls, which keep the reference's `Series` label);
#'   `comparisons` (data frame, one row per non-reference series in file order:
#'   `series`, `n`, `obj1`, `obj2`, `chi`, `df`, `p`, `ec50_reference`,
#'   `ec50_series`, `ec50_shared`, `ratio`, `significant`, `direction`,
#'   `controls_borrowed`, `ec50_extrapolated`); `shared` (run-2 fit per
#'   comparison: `par` = named
#'   `c(max_reference, slope_reference, max_series, slope_series, ec50)`,
#'   `objective`, `convergence`); and `summary` (conclusion text).
#' @export
analyse_potentiation <- function(df, reference = NULL,
                                 response = c("continuous", "quantal"),
                                 alpha = 0.05, stressor1 = "Stressor 1") {
  response <- match.arg(response)
  df <- as.data.frame(df)
  cols <- c("Series", "C1",
            if (response == "continuous") "Res" else c("Exposed", "Affected"))
  missing <- setdiff(cols, names(df))
  if (length(missing))
    stop("analyse_potentiation: missing column(s) for ", response, " data: ",
         paste(missing, collapse = ", "))
  df <- df[cols]
  df$Series <- trimws(as.character(df$Series))

  series <- unique(df$Series)
  if (length(series) < 2)
    stop("analyse_potentiation needs at least two series (a reference and one ",
         "to compare); found ", length(series), ".")
  reference <- if (is.null(reference)) series[[1]] else
    trimws(as.character(reference))
  if (!reference %in% series)
    stop("analyse_potentiation: reference series '", reference,
         "' is not in the data (series: ", paste(series, collapse = ", "), ").")
  others <- setdiff(series, reference)                  # file order

  own <- lapply(stats::setNames(series, series), function(s)
    df[df$Series == s, , drop = FALSE])
  ref_controls <- own[[reference]][own[[reference]]$C1 == 0, , drop = FALSE]
  borrowed <- vapply(series, function(s)
    s != reference && !any(own[[s]]$C1 == 0) && nrow(ref_controls) > 0,
    logical(1))
  used <- lapply(stats::setNames(series, series), function(s) {
    d <- if (borrowed[[s]]) rbind(ref_controls, own[[s]]) else own[[s]]
    rownames(d) <- NULL
    d
  })
  too_small <- series[vapply(used, nrow, integer(1)) <= 3]
  if (length(too_small))
    stop("analyse_potentiation: each series needs more than 3 rows to fit a ",
         "curve; too few in: ", paste(too_small, collapse = ", "), ".")

  # Run 1: one curve per series; the reference fit is shared by every comparison.
  fits <- lapply(used, pot_fit_curve, response = response)
  # Run 2: shared EC50, per comparison.
  shared <- lapply(stats::setNames(others, others), function(k)
    pot_fit_shared(used[[reference]], used[[k]], fits[[reference]], fits[[k]],
                   response))
  # Run 2 is nested in run 1, so obj1 <= obj2 must hold. Re-seed each run-1
  # curve from the run-2 solution and keep the better fit: a descent from that
  # point can only lower its share of obj2, which guarantees the ordering.
  for (k in others) {
    sp <- shared[[k]]$par
    fits[[k]] <- pot_better(fits[[k]], pot_fit_curve(
      used[[k]], response, list(c(max = sp[["max_series"]],
                                  slope = sp[["slope_series"]],
                                  ec50 = sp[["ec50"]]))))
    fits[[reference]] <- pot_better(fits[[reference]], pot_fit_curve(
      used[[reference]], response, list(c(max = sp[["max_reference"]],
                                          slope = sp[["slope_reference"]],
                                          ec50 = sp[["ec50"]]))))
  }

  curves <- lapply(stats::setNames(series, series), function(s) {
    f <- fits[[s]]
    c(f, list(n = nrow(used[[s]]), controls_borrowed = borrowed[[s]],
              ec50_extrapolated = pot_extrapolated(f$par[["ec50"]], own[[s]]$C1),
              data = used[[s]]))
  })

  ref <- curves[[reference]]
  comparisons <- do.call(rbind, lapply(others, function(k) {
    cur <- curves[[k]]
    n <- ref$n + cur$n
    obj1 <- ref$objective + cur$objective
    obj2 <- shared[[k]]$objective
    chi <- if (isTRUE(obj2 > obj1))
      lr_test(obj2, obj1, 5, 6, n,
              if (response == "continuous") "continuous" else "binary")$chi
    else 0
    p <- stats::pchisq(chi, df = 1, lower.tail = FALSE)
    ratio <- cur$par[["ec50"]] / ref$par[["ec50"]]
    significant <- p < alpha
    direction <- if (significant && ratio < 1) "more toxic" else
      if (significant && ratio > 1) "less toxic" else "no significant influence"
    data.frame(series = k, n = n, obj1 = obj1, obj2 = obj2, chi = chi,
               df = 1L, p = p, ec50_reference = ref$par[["ec50"]],
               ec50_series = cur$par[["ec50"]],
               ec50_shared = shared[[k]]$par[["ec50"]], ratio = ratio,
               significant = significant, direction = direction,
               controls_borrowed = cur$controls_borrowed,
               ec50_extrapolated = cur$ec50_extrapolated,
               stringsAsFactors = FALSE)
  }))
  rownames(comparisons) <- NULL

  res <- list(reference = reference, response = response, alpha = alpha,
              series = series, curves = curves, comparisons = comparisons,
              shared = shared, summary = NA_character_)
  res$summary <- potentiation_summary(res, stressor1)
  res
}

# Parameter box shared by the separate and the shared-EC50 fits (fit_single()'s
# defaults; a quantal max is a probability, so it is capped at 1).
pot_bounds <- function(response) {
  list(lower = c(max = 1e-8, slope = 1e-3, ec50 = 1e-8),
       upper = c(max = if (response == "quantal") 1 else Inf, slope = 50,
                 ec50 = Inf))
}

# Objective of one curve on one series: residual SS (continuous) or binomial
# deviance with predictions clamped into (0, 1) (quantal).
pot_objective <- function(d, max, slope, ec50, response) {
  pred <- ll3_predict(d$C1, max, slope, ec50)
  if (response == "continuous") obj_ss(d$Res, pred)
  else obj_deviance(d$Exposed, d$Affected, pmin(pmax(pred, 1e-8), 1 - 1e-8))
}

# Exact gradient of pot_objective() in (max, slope, ec50). Steep curves (slope
# near its bound of 50, as in the example data) put the optimum in a narrow,
# curved valley where optim's finite-difference gradient is too coarse and
# L-BFGS-B stalls short of the minimum (level 3's shared fit stops ~2 SS units
# high); the exact gradient reaches it from the same seeds. With pred = max * w,
# w = 1 / (1 + (C/ec50)^slope): d pred / d slope = -max w (1 - w) log(C/ec50)
# and d pred / d ec50 = max w (1 - w) slope / ec50; C = 0 rows contribute
# nothing to either. A clamped quantal prediction has zero gradient.
pot_gradient <- function(d, max, slope, ec50, response) {
  w <- 1 / (1 + (d$C1 / ec50)^slope)
  pred <- max * w
  log_ratio <- ifelse(d$C1 > 0, log(d$C1 / ec50), 0)
  dpred <- cbind(max = w, slope = -max * w * (1 - w) * log_ratio,
                 ec50 = max * w * (1 - w) * slope / ec50)
  dobj <- if (response == "continuous") {
    -2 * (d$Res - pred)
  } else {
    p <- pmin(pmax(pred, 1e-8), 1 - 1e-8)
    ifelse(pred > 1e-8 & pred < 1 - 1e-8,
           -2 * (d$Affected / p - (d$Exposed - d$Affected) / (1 - p)), 0)
  }
  colSums(dobj * dpred)
}

# One L-BFGS-B descent, configured like fit_single() but with the exact
# gradient `gr`. A non-finite objective is turned into a large finite one so
# the optimiser can back off instead of failing; if it fails anyway the
# (clamped) start is returned with convergence = 99, so a re-seeded fit is
# never worse than its seed.
pot_optim <- function(start, fn, gr, lower, upper) {
  start <- pmin(pmax(start, lower), upper)
  safe_fn <- function(p) { v <- fn(p); if (is.finite(v)) v else 1e12 }
  res <- tryCatch(
    stats::optim(start, safe_fn, gr, method = "L-BFGS-B", lower = lower,
                 upper = upper,
                 control = list(parscale = pmax(abs(start), 1e-8),
                                factr = 1e-9, maxit = 1000)),
    error = function(e) NULL)
  if (is.null(res)) res <- list(par = start, value = safe_fn(start),
                                convergence = 99L)
  res$par <- stats::setNames(res$par, names(start))
  res
}

# Run-1 fit of one series: three free parameters, best of `starts` (default:
# fit_single()'s data-driven start).
pot_fit_curve <- function(d, response, starts = NULL) {
  b <- pot_bounds(response)
  if (is.null(starts)) {
    obs <- if (response == "continuous") d$Res else d$Affected / d$Exposed
    pos <- d$C1[d$C1 > 0]
    starts <- list(c(max = max(obs), slope = 1,
                     ec50 = if (length(pos)) stats::median(pos) else 1))
  }
  fn <- function(p) pot_objective(d, p[["max"]], p[["slope"]], p[["ec50"]],
                                  response)
  gr <- function(p) pot_gradient(d, p[["max"]], p[["slope"]], p[["ec50"]],
                                 response)
  best <- NULL
  for (st in starts) {
    res <- pot_optim(st[c("max", "slope", "ec50")], fn, gr, b$lower, b$upper)
    if (is.null(best) || res$value < best$value) best <- res
  }
  list(par = best$par, objective = best$value, convergence = best$convergence)
}

# Run-2 fit of one comparison: max and slope per curve, one shared EC50. Seeded
# from run 1 with the shared EC50 at the reference's, the series', and their
# geometric mean; the lowest objective wins.
pot_fit_shared <- function(d_ref, d_k, fit_ref, fit_k, response) {
  b <- pot_bounds(response)
  lower <- c(max_reference = b$lower[["max"]], slope_reference = b$lower[["slope"]],
             max_series = b$lower[["max"]], slope_series = b$lower[["slope"]],
             ec50 = b$lower[["ec50"]])
  upper <- c(max_reference = b$upper[["max"]], slope_reference = b$upper[["slope"]],
             max_series = b$upper[["max"]], slope_series = b$upper[["slope"]],
             ec50 = b$upper[["ec50"]])
  fn <- function(p)
    pot_objective(d_ref, p[["max_reference"]], p[["slope_reference"]],
                  p[["ec50"]], response) +
    pot_objective(d_k, p[["max_series"]], p[["slope_series"]], p[["ec50"]],
                  response)
  gr <- function(p) {
    g_ref <- pot_gradient(d_ref, p[["max_reference"]], p[["slope_reference"]],
                          p[["ec50"]], response)
    g_k <- pot_gradient(d_k, p[["max_series"]], p[["slope_series"]],
                        p[["ec50"]], response)
    c(g_ref[["max"]], g_ref[["slope"]], g_k[["max"]], g_k[["slope"]],
      g_ref[["ec50"]] + g_k[["ec50"]])
  }
  curve_par <- c(max_reference = fit_ref$par[["max"]],
                 slope_reference = fit_ref$par[["slope"]],
                 max_series = fit_k$par[["max"]],
                 slope_series = fit_k$par[["slope"]])
  e_ref <- fit_ref$par[["ec50"]]
  e_k <- fit_k$par[["ec50"]]
  best <- NULL
  for (e in c(e_ref, e_k, sqrt(e_ref * e_k))) {
    res <- pot_optim(c(curve_par, ec50 = e), fn, gr, lower, upper)
    if (is.null(best) || res$value < best$value) best <- res
  }
  list(par = best$par, objective = best$value, convergence = best$convergence)
}

# The lower-objective of two run-1 fits.
pot_better <- function(a, b) if (b$objective < a$objective) b else a

# Is the EC50 outside the series' own tested range (above its highest C1 or
# below its lowest positive C1)? A series with no positive dose cannot locate
# its EC50 at all, which counts as extrapolated.
pot_extrapolated <- function(ec50, conc) {
  pos <- conc[conc > 0]
  if (!length(pos)) return(TRUE)
  ec50 > max(pos) || ec50 < min(pos)
}
