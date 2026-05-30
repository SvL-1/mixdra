#' Heuristic / fitted seed for a single chemical's curve
#'
#' When the marginal series has enough distinct points (>= 4) a full
#' [fit_single()] is used; otherwise (tiny grids) it falls back to a cheap
#' heuristic (`max` = largest response, `slope` = 1, `ec50` = median positive
#' concentration). Either way it never errors, so seeding is robust on small
#' synthetic designs as well as full experimental data.
#' @keywords internal
seed_one <- function(conc, resp) {
  ok <- sum(!is.na(conc) & !is.na(resp))
  if (ok >= 4 && length(unique(conc[!is.na(conc)])) >= 4) {
    f <- tryCatch(fit_single(conc, resp), error = function(e) NULL)
    if (!is.null(f))
      return(c(max = f$par[["max"]], slope = f$par[["slope"]],
               ec50 = f$par[["ec50"]]))
  }
  pos <- conc[!is.na(conc) & conc > 0]
  c(max = max(resp, na.rm = TRUE),
    slope = 1,
    ec50 = if (length(pos)) stats::median(pos) else 1)
}

#' Seed starting values from per-chemical single-curve fits
#'
#' Fits/estimates one log-logistic curve per present chemical (chemical *i* from
#' the rows where every other concentration is 0) and names the seeds to match
#' the active registry (`slope1..n` / `ec50..n`).
#' @keywords internal
seed_from_singles <- function(df, response) {
  resp_col <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  n <- length(cols)
  pnames <- model_spec("CA", "reference", n)$params  # max, slope1..n, ec50..n

  seeds <- lapply(seq_along(cols), function(i) {
    others <- cols[-i]
    keep <- if (length(others) == 0) {
      rep(TRUE, nrow(df))
    } else {
      rowSums(df[others] == 0) == length(others)  # all other chemicals at 0
    }
    seed_one(df[[cols[i]]][keep], resp_col[keep])
  })

  maxv   <- mean(vapply(seeds, function(s) s[["max"]],   numeric(1)))
  slopes <- vapply(seeds, function(s) s[["slope"]], numeric(1))
  ec50s  <- vapply(seeds, function(s) s[["ec50"]],  numeric(1))
  stats::setNames(c(maxv, slopes, ec50s), pnames)
}

#' Stage-1 fit: estimate the curve parameters from single-compound data
#'
#' Fits the `reference` model to the control + single-compound rows only (rows
#' where at most one chemical is present). On those rows every mixture model
#' reduces exactly to the independent log-logistic curves, so this identifies the
#' shared `max` and per-chemical `slope`/`ec50` without any interaction
#' parameter. Falls back to the [seed_from_singles()] heuristic if the subset
#' cannot be fit (e.g. a chemical with too sparse a marginal series).
#' @keywords internal
fit_curve_from_singles <- function(df, reference, response,
                                   lower = NULL, upper = NULL,
                                   n_starts = 1, time_limit = 30) {
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  singles <- df[rowSums(df[cols] > 0) <= 1, , drop = FALSE]
  # Seed from the full frame (seed_from_singles already isolates each chemical's
  # marginal rows internally); fit_model then refines on the single-compound subset.
  seed <- seed_from_singles(df, response)
  fit <- tryCatch(
    fit_model(singles, reference, "reference", response, start = seed,
              lower = lower, upper = upper, n_starts = n_starts,
              time_limit = time_limit),
    error = function(e) NULL)
  if (is.null(fit) || !all(is.finite(fit$par[names(seed)]))) return(seed)
  fit$par[names(seed)]
}

#' Analyse a mixture: fit reference + deviations and compare
#'
#' Fitting is staged: the curve parameters are fixed from the single-compound
#' data, then only the interaction parameters (`a`, `b`, ...) are fitted to the
#' mixture data. Each fit's `df` therefore counts only its free interaction
#' parameters; likelihood-ratio tests use differences, so they are unaffected.
#'
#' @param df Mixture data frame (see [fit_model()]).
#' @param reference "CA" or "IA".
#' @param response "continuous" or "binary".
#' @param start Optional named vector of the curve parameters (`max`, `slope*`,
#'   `ec50*`) to hold fixed. When `NULL` (default) they are estimated from the
#'   single-compound data via [fit_curve_from_singles()].
#' @param alpha Significance threshold for model selection.
#' @param lower,upper Optional named numeric vectors of hard parameter bounds,
#'   forwarded to every [fit_model()] call (the four models share the base
#'   parameters). See [fit_model()] for the bound semantics.
#' @param n_starts Number of optimisation starts per model, forwarded to
#'   [fit_model()]. Defaults to 1 (single start) for speed during testing;
#'   raise it (e.g. 20) to multi-start and more reliably escape the local
#'   minima of the CA bisection surface at the cost of runtime.
#' @param time_limit Per-model wall-clock budget in seconds, forwarded to
#'   [fit_model()]. Applies to each of the four fits independently, so a full
#'   analysis can take up to `4 * time_limit`. `NULL` disables the limit.
#' @return A list: `fits` (named list of model fits), `comparison` (data frame of
#'   LR tests vs each model's parent), `chosen` (selected model name),
#'   `reference`, `response`.
#' @export
analyse_mixture <- function(df, reference, response = c("continuous", "binary"),
                            start = NULL, alpha = 0.05,
                            lower = NULL, upper = NULL, n_starts = 1,
                            time_limit = 30) {
  response <- match.arg(response)

  # Stage 1: estimate the curve parameters (max, slope*, ec50*) from the
  # single-compound data ONLY and hold them fixed. A user-supplied `start` is
  # taken as the fixed curve parameters directly. The curve params must not be
  # influenced by the mixture rows -- that is the whole point of staged fitting.
  base <- if (is.null(start)) {
    fit_curve_from_singles(df, reference, response, lower, upper,
                           n_starts, time_limit)
  } else {
    start
  }
  base_names <- names(base)

  # Stage 2: with the curve parameters fixed, fit only the interaction
  # parameters (a, b, ...) to the full mixture data, once per deviation.
  devs <- c("reference", "SA", "DR", "DL")
  fits <- lapply(devs, function(d)
    fit_model(df, reference, d, response, start = base, fixed = base_names,
              lower = lower, upper = upper, n_starts = n_starts,
              time_limit = time_limit))
  names(fits) <- devs

  n <- nrow(df)
  parent_of <- c(SA = "reference", DR = "SA", DL = "SA")
  comparison <- do.call(rbind, lapply(names(parent_of), function(m) {
    p <- parent_of[[m]]
    t <- lr_test(fits[[p]]$objective, fits[[m]]$objective,
                 fits[[p]]$df, fits[[m]]$df, n, response)
    data.frame(model = m, parent = p, chi = t$chi, df = t$df, p = t$p)
  }))

  list(fits = fits, comparison = comparison,
       chosen = select_parsimonious(fits, n, response, alpha),
       reference = reference, response = response)
}

#' Analyse a single chemical (dose-response curve only)
#' @param df Data frame with `C1` and `Res` (or `Exposed`/`Affected`).
#' @return The [fit_single()] result.
#' @export
analyse_single <- function(df) {
  resp <- if ("Res" %in% names(df)) df$Res else df$Affected / df$Exposed
  fit_single(df$C1, resp)
}
