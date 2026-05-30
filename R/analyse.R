#' Seed starting values from per-chemical single-curve fits
#'
#' Fits one log-logistic curve per present chemical (chemical *i* seeded from the
#' rows where every other concentration is 0) and names the seeds to match the
#' active registry (`slope1..n` / `ec50..n`).
#' @keywords internal
seed_from_singles <- function(df, response) {
  resp_col <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  n <- length(cols)
  pnames <- model_spec("CA", "reference", n)$params  # max, slope1..n, ec50..n

  fits <- lapply(seq_along(cols), function(i) {
    others <- cols[-i]
    keep <- if (length(others) == 0) {
      rep(TRUE, nrow(df))
    } else {
      rowSums(df[others] == 0) == length(others)  # all other chemicals at 0
    }
    fit_single(df[[cols[i]]][keep], resp_col[keep])
  })

  maxv   <- mean(vapply(fits, function(f) f$par[["max"]],   numeric(1)))
  slopes <- vapply(fits, function(f) f$par[["slope"]], numeric(1))
  ec50s  <- vapply(fits, function(f) f$par[["ec50"]],  numeric(1))
  stats::setNames(c(maxv, slopes, ec50s), pnames)
}

#' Analyse a mixture: fit reference + deviations and compare
#'
#' @param df Mixture data frame (see [fit_model()]).
#' @param reference "CA" or "IA".
#' @param response "continuous" or "binary".
#' @param start Optional named starting vector; if NULL, seeded from single fits.
#' @param alpha Significance threshold for model selection.
#' @param n_starts Number of optimisation starts per model, forwarded to
#'   [fit_model()]. The default uses multi-start to reliably escape the local
#'   minima of the CA bisection surface.
#' @return A list: `fits` (named list of model fits), `comparison` (data frame of
#'   LR tests vs each model's parent), `chosen` (selected model name),
#'   `reference`, `response`.
#' @export
analyse_mixture <- function(df, reference, response = c("continuous", "binary"),
                            start = NULL, alpha = 0.05, n_starts = 20) {
  response <- match.arg(response)
  if (is.null(start)) start <- seed_from_singles(df, response)

  devs <- c("reference", "SA", "DR", "DL")
  fits <- lapply(devs, function(d)
    fit_model(df, reference, d, response, start = start, n_starts = n_starts))
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
