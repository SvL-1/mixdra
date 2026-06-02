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

#' Model order for the joint selection chain
#'
#' Binary mixtures walk reference -> SA -> {DR, DL}; ternary mixtures support
#' only reference -> SA here (Advanced S/A is fitted separately).
#' @keywords internal
joint_chain_order <- function(n_chem) {
  if (n_chem == 2) c("reference", "SA", "DR", "DL") else c("reference", "SA")
}

#' Fit one mixture model jointly (curves + interaction), warm-started
#'
#' Unlike the staged [fit_curve_from_singles()] + fixed-curve fits used by
#' [analyse_mixture()], this leaves EVERY parameter free and optimises them
#' together (the Excel-faithful joint fit). It is seeded from `parent_fit` when
#' given (the parent in the nesting chain), else from `seed_curves`. The child's
#' new interaction parameter (the one its parent lacks) is seeded at its neutral
#' value `0` -- `a = 0` reduces SA to the reference, `b = 0` reduces DR/DL to SA
#' -- so the child starts at a point reproducing the parent's fit and can only
#' improve on it (keeping the nested LR comparison sound).
#' @param df Mixture data frame (see [fit_model()]).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param seed_curves Named numeric of curve params (max, slope*, ec50*) used as
#'   the start when there is no parent (the reference fit).
#' @param parent_fit The parent model's joint fit, or `NULL` for the reference.
#' @param n_starts,time_limit Forwarded to [fit_model()].
#' @return An enriched mixture fit (as [fit_model()]).
#' @keywords internal
joint_fit_one <- function(df, reference, deviation, response,
                          seed_curves, parent_fit = NULL,
                          n_starts = 1, time_limit = 30) {
  n_chem <- length(intersect(c("C1", "C2", "C3"), names(df)))
  spec   <- model_spec(reference, deviation, n_chem)
  start  <- if (is.null(parent_fit)) seed_curves else parent_fit$par
  # Neutral-seed any interaction parameter the child adds beyond its start
  # (a or b): 0 makes the new term vanish, reproducing the parent.
  new_extra <- setdiff(spec$extra, names(start))
  if (length(new_extra))
    start <- c(start, stats::setNames(rep(0, length(new_extra)), new_extra))
  fit_model(df, reference, deviation, response, start = start,
            n_starts = n_starts, time_limit = time_limit)
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
#'   forwarded to every [fit_model()] call (all deviations share the base
#'   parameters). See [fit_model()] for the bound semantics.
#' @param n_starts Number of optimisation starts per model, forwarded to
#'   [fit_model()]. Defaults to 1 (single start) for speed during testing;
#'   raise it (e.g. 20) to multi-start and more reliably escape the local
#'   minima of the CA bisection surface at the cost of runtime.
#' @param time_limit Per-model wall-clock budget in seconds, forwarded to
#'   [fit_model()]. Applies to each fit independently (four for binary:
#'   reference/SA/DR/DL; two for ternary: reference/SA), so a full analysis can
#'   take up to `n_fits * time_limit`. `NULL` disables the limit.
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
  # DR and DL are binary-only; ternary mixtures use reference / S/A only here
  # (Advanced S/A is fitted separately via analyse_ternary()).
  n_chem <- length(intersect(c("C1", "C2", "C3"), names(df)))
  devs <- if (n_chem == 2) c("reference", "SA", "DR", "DL") else c("reference", "SA")
  fits <- lapply(devs, function(d)
    fit_model(df, reference, d, response, start = base, fixed = base_names,
              lower = lower, upper = upper, n_starts = n_starts,
              time_limit = time_limit))
  names(fits) <- devs

  cmp <- compare_fits(fits, nrow(df), response, alpha)
  list(fits = fits, comparison = cmp$comparison, chosen = cmp$chosen,
       reference = reference, response = response)
}

#' Analyse a single chemical (dose-response curve only)
#' @param df Data frame with `C1` and `Res` (or `Exposed`/`Affected`).
#' @param lower,upper,start Optional named numeric vectors forwarded to
#'   [fit_single()] (`max`/`slope`/`ec50`).
#' @return The [fit_single()] result.
#' @export
analyse_single <- function(df, lower = NULL, upper = NULL, start = NULL) {
  resp <- if ("Res" %in% names(df)) df$Res else df$Affected / df$Exposed
  fit_single(df$C1, resp, lower = lower, upper = upper, start = start)
}
