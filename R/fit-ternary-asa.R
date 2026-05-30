#' Classify mixture rows by how many chemicals are present
#' @param df Data frame with C1, C2 (and C3) concentration columns.
#' @return A factor: "control" (all 0), "single", "binary", or "ternary".
#' @keywords internal
classify_rows <- function(df) {
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  M <- as.matrix(df[cols])
  nz <- rowSums(M > 0)
  labs <- c("0" = "control", "1" = "single", "2" = "binary", "3" = "ternary")
  factor(unname(labs[as.character(nz)]),
         levels = c("control", "single", "binary", "ternary"))
}

#' Exact nominal mixture-ratio key for ternary rows
#'
#' A ratio is the fixed proportion C1:C2:C3 (a ray from the origin); a
#' dose-response series at one ratio holds the proportions constant. Users supply
#' NOMINAL concentrations, so grouping is exact — the only rounding (to `sig`
#' significant figures) neutralises floating-point representation across dose
#' levels, NOT measurement noise.
#' @param df Data frame of ternary rows (all of C1, C2, C3 > 0).
#' @param sig Significant figures for representation rounding (default 6).
#' @return Character vector of ratio keys, one per row.
#' @keywords internal
ternary_ratio_key <- function(df, sig = 6) {
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  M <- as.matrix(df[cols])
  p <- M / rowSums(M)
  apply(signif(p, sig), 1, paste, collapse = "_")
}

#' Staged Advanced S/A ternary fit (internal)
#'
#' Stage 1: base curve params from singles (+control). Stage 2: A1/A2/A3 from
#' binaries (base fixed, start 0). Stage 3a: overall A4 from all data (base +
#' A1/A2/A3 fixed). Stage 3b: individual A4 per ternary ratio.
#' @param df Mixture data frame (C1,C2,C3 + Res for continuous).
#' @param reference "CA" (IA not yet supported here).
#' @param response "continuous" (binary not yet validated for ASA).
#' @param lower,upper,n_starts,time_limit Forwarded to [fit_model()].
#' @return See [analyse_ternary()].
#' @keywords internal
fit_ternary_asa <- function(df, reference = "CA", response = "continuous",
                            lower = NULL, upper = NULL, n_starts = 1,
                            time_limit = 30) {
  if (reference != "CA")
    stop("fit_ternary_asa: only reference = 'CA' is implemented")
  cls <- classify_rows(df)
  base_params <- model_spec(reference, "reference", 3)$params

  # Stage 1 - base curve params from singles + control.
  singles <- df[cls %in% c("control", "single"), , drop = FALSE]
  seed <- seed_from_singles(df, response)
  f1 <- fit_model(singles, reference, "reference", response, start = seed,
                  lower = lower, upper = upper, n_starts = n_starts,
                  time_limit = time_limit)
  base_par <- f1$par[base_params]

  # Stage 2 - A1/A2/A3 from binaries; base + A4 held fixed; start all A at 0.
  binaries <- df[cls == "binary", , drop = FALSE]
  f2 <- fit_model(binaries, reference, "ASA", response,
                  start = c(base_par, A1 = 0, A2 = 0, A3 = 0, A4 = 0),
                  fixed = c(base_params, "A4"),
                  n_starts = n_starts, time_limit = time_limit)
  a123 <- f2$par[c("A1", "A2", "A3")]

  # Stage 3a - overall A4 from ALL residuals; everything else fixed.
  start3 <- c(base_par, a123, A4 = 0)
  f_overall <- fit_model(df, reference, "ASA", response, start = start3,
                         fixed = c(base_params, "A1", "A2", "A3"),
                         n_starts = n_starts, time_limit = time_limit)
  A4_overall <- f_overall$par[["A4"]]

  # Stage 3b - individual A4 per ternary ratio.
  tern_idx <- which(cls == "ternary")
  keys <- ternary_ratio_key(df[tern_idx, , drop = FALSE])
  groups <- split(tern_idx, keys)
  individual_fits <- lapply(names(groups), function(k) {
    sub <- df[groups[[k]], , drop = FALSE]
    fr <- fit_model(sub, reference, "ASA", response, start = start3,
                    fixed = c(base_params, "A1", "A2", "A3"),
                    n_starts = n_starts, time_limit = time_limit)
    prop <- as.numeric(strsplit(k, "_", fixed = TRUE)[[1]])
    list(ratio = k, A4 = fr$par[["A4"]], n = nrow(sub),
         C1 = prop[1], C2 = prop[2], C3 = prop[3], fit = fr)
  })
  names(individual_fits) <- names(groups)

  individual <- do.call(rbind, lapply(individual_fits, function(x)
    data.frame(ratio = x$ratio, C1 = x$C1, C2 = x$C2, C3 = x$C3,
               A4 = x$A4, n = x$n, stringsAsFactors = FALSE)))
  rownames(individual) <- NULL

  list(reference = reference, response = response,
       base = base_par, pairwise = a123, A4_overall = A4_overall,
       individual = individual,
       fits = list(base = f1, pairwise = f2, overall = f_overall,
                   individual = individual_fits))
}
