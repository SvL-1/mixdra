#' Build a default mixture design grid anchored to the model's EC50s
#'
#' Produces a data frame of concentration columns (`C1`..`Cn`) covering: one
#' control row (all zero); a per-chemical single-compound geometric ladder
#' spanning `ec50 * mult` (others zero) so each marginal curve is identifiable
#' and the staged fit has single-compound rows; and fixed-ratio mixture rays,
#' each a geometric dose series in toxic units. Concentrations are scaled by each
#' chemical's EC50 (read from `par`).
#' @param par Named parameter vector containing the curve parameters; the number
#'   of `slope*` entries determines `n_chem` (2 or 3) and the `ec50*` entries set
#'   the concentration scale.
#' @param reference "CA" or "IA" (only used to resolve EC50 parameter names).
#' @param deviation Deviation name (only used to resolve parameter names).
#' @param ratios Optional list of mixture-fraction vectors (length `n_chem`).
#'   Defaults: binary `list(c(1,1), c(1,3), c(3,1))`; ternary
#'   `list(c(1,1,1), c(1,1,0), c(1,0,1), c(0,1,1))`.
#' @param n_per_ray Number of dose points per mixture ray.
#' @param mult Multipliers for the single-compound ladder (relative to EC50).
#' @return A data frame with columns `C1`..`Cn`.
#' @export
mixture_design <- function(par, reference = "CA", deviation = "reference",
                           ratios = NULL, n_per_ray = 7, mult = 2^(-3:3)) {
  n_chem <- sum(grepl("^slope[0-9]+$", names(par)))
  if (!n_chem %in% c(2, 3))
    stop("mixture_design: `par` must contain 2 or 3 slope* entries; got ", n_chem)
  spec <- model_spec(reference, deviation, n_chem)
  ec_names <- grep("^ec50", spec$params, value = TRUE)   # length n_chem, in order
  ec50 <- as.numeric(par[ec_names])
  if (anyNA(ec50))
    stop("mixture_design: `par` must supply every EC50 (",
         paste(ec_names, collapse = ", "), ")")
  cols <- paste0("C", seq_len(n_chem))

  zero_row <- function() stats::setNames(as.list(rep(0, n_chem)), cols)
  rows <- list(zero_row())                               # control

  for (i in seq_len(n_chem)) {                            # single-compound ladders
    for (cv in ec50[i] * mult) {
      r <- zero_row(); r[[cols[i]]] <- cv
      rows[[length(rows) + 1L]] <- r
    }
  }

  if (is.null(ratios))
    ratios <- if (n_chem == 2) list(c(1, 1), c(1, 3), c(3, 1))
              else list(c(1, 1, 1), c(1, 1, 0), c(1, 0, 1), c(0, 1, 1))

  dose <- 2^seq(-3, 3, length.out = n_per_ray)            # toxic-unit dose series
  for (w in ratios) {
    w <- w / sum(w)                                       # mixture fractions
    for (d in dose) {
      conc <- d * w * ec50                                # conc_i = d * w_i * ec50_i
      rows[[length(rows) + 1L]] <- stats::setNames(as.list(conc), cols)
    }
  }

  out <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(out) <- NULL
  out
}

#' Apply the response-appropriate noise layer to expected responses
#'
#' Continuous: additive Gaussian with SD = `cv * |mu|` (relative CV; `cv = 0`
#' returns `mu` exactly). Binary: when `group_size` is infinite, returns exact
#' (fractional) counts `Affected = 100 * mu` with `Exposed = 100`; otherwise
#' draws `Affected ~ Binomial(group_size, mu)`.
#' @param mu Numeric vector of expected responses (probabilities for binary).
#' @param response "continuous" or "binary".
#' @param cv Relative noise (continuous only).
#' @param group_size Binomial group size (binary only); `Inf` = exact.
#' @return A named list of response columns to bind onto the design frame.
#' @keywords internal
.apply_noise <- function(mu, response, cv, group_size) {
  if (response == "continuous") {
    res <- if (cv > 0) mu + stats::rnorm(length(mu), 0, cv * abs(mu)) else mu
    list(Res = res)
  } else if (is.infinite(group_size)) {
    list(Exposed = rep(100, length(mu)), Affected = 100 * mu)
  } else {
    list(Exposed = rep(group_size, length(mu)),
         Affected = stats::rbinom(length(mu), group_size, mu))
  }
}

#' Simulate a mixture dose-response dataset from known parameters
#'
#' Forward generator for verification by parameter recovery: builds (or accepts)
#' a design grid, predicts the noise-free expected response through the same
#' model the fitter uses ([mixture_predict()]), then applies an optional noise
#' layer. At `cv = 0` (continuous) or `group_size = Inf` (binary) the output is
#' deterministic and exact, so a downstream [analyse_mixture()] should recover
#' `par`.
#' @param par Named full parameter vector matching
#'   `model_spec(reference, deviation, n_chem)$params`.
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param design Optional concentration design (`C1`..`Cn`); defaults to
#'   [mixture_design()].
#' @param cv Relative Gaussian noise for continuous responses (0 = none).
#' @param group_size Binomial group size for binary responses; `Inf` (default)
#'   yields exact proportions.
#' @param reps Replicate each design row this many times before adding noise.
#'   With `cv = 0` (continuous) or `group_size = Inf` (binary) the replicates are
#'   exact duplicate rows.
#' @param seed Optional RNG seed for reproducible noisy draws. Restores the global
#'   RNG state on exit, so a seeded call does not disturb the caller's stream.
#' @return A data frame the engine consumes unchanged: `C1`..`Cn` plus `Res`
#'   (continuous) or `Exposed`/`Affected` (binary).
#' @export
simulate_mixture <- function(par, reference = "CA", deviation = "reference",
                             response = c("continuous", "binary"),
                             design = NULL, cv = 0, group_size = Inf,
                             reps = 1, seed = NULL) {
  response <- match.arg(response)
  n_chem <- sum(grepl("^slope[0-9]+$", names(par)))
  spec <- model_spec(reference, deviation, n_chem)
  extra <- setdiff(names(par), spec$params)
  if (length(extra))
    stop("simulate_mixture: `par` has unknown parameter(s) for ", reference, "/",
         deviation, ": ", paste(extra, collapse = ", "))
  if (response == "binary" && "max" %in% names(par) && par[["max"]] > 1)
    stop("simulate_mixture: binary `max` is a probability and must be <= 1")
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv),
            add = TRUE)
    set.seed(seed)
  }
  if (is.null(design)) design <- mixture_design(par, reference, deviation)
  if (reps > 1) design <- design[rep(seq_len(nrow(design)), reps), , drop = FALSE]
  rownames(design) <- NULL

  mu <- mixture_predict(design, par, reference, deviation)
  cbind(design, as.data.frame(.apply_noise(unname(mu), response, cv, group_size)))
}
