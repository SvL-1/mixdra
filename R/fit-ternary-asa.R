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
#' @param base Optional named curve params (`max, slope1-3, ec50_1-3`); when
#'   supplied the Stage-1 singles fit is skipped and these are held fixed.
#'   Defaults to `NULL` (fit from singles).
#' @param pairwise Optional named numeric vector (`A1, A2, A3`); when supplied
#'   the Stage-2 binaries fit is skipped and these are held fixed.
#'   Defaults to `NULL` (fit from binaries).
#' @param lower,upper,n_starts,time_limit Forwarded to [fit_model()].
#' @return See [analyse_ternary()].
#' @keywords internal
fit_ternary_asa <- function(df, reference = "CA", response = "continuous",
                            base = NULL, pairwise = NULL, lower = NULL,
                            upper = NULL, n_starts = 1, time_limit = 30) {
  if (reference != "CA")
    stop("fit_ternary_asa: only reference = 'CA' is implemented")
  cls <- classify_rows(df)
  base_params <- model_spec(reference, "reference", 3)$params

  # Stage 1 - base curve params: supplied (frozen, from the app's three single
  # curves) or fit from the singles + control. A supplied `base` is held fixed
  # downstream exactly as an auto-fit one would be.
  if (is.null(base)) {
    singles <- df[cls %in% c("control", "single"), , drop = FALSE]
    seed <- seed_from_singles(df, response)
    f1 <- fit_model(singles, reference, "reference", response, start = seed,
                    lower = lower, upper = upper, n_starts = n_starts,
                    time_limit = time_limit)
    base_par <- f1$par[base_params]
  } else {
    if (!is.numeric(base))
      stop("fit_ternary_asa: `base` must be a named numeric vector")
    miss <- setdiff(base_params, names(base))
    if (length(miss))
      stop("fit_ternary_asa: `base` is missing param(s): ",
           paste(miss, collapse = ", "))
    base_par <- base[base_params]
    f1 <- NULL
  }

  # Stage 2 - A1/A2/A3 from binaries; base + A4 held fixed; start all A at 0.
  # A supplied `pairwise` (the campaign's per-pair S/A values) skips this fit and
  # holds those values fixed downstream -- exact, because each binary row
  # activates exactly one A term (see the design doc, section 3).
  if (is.null(pairwise)) {
    binaries <- df[cls == "binary", , drop = FALSE]
    f2 <- fit_model(binaries, reference, "ASA", response,
                    start = c(base_par, A1 = 0, A2 = 0, A3 = 0, A4 = 0),
                    fixed = c(base_params, "A4"),
                    n_starts = n_starts, time_limit = time_limit)
    a123 <- f2$par[c("A1", "A2", "A3")]
  } else {
    if (!is.numeric(pairwise))
      stop("fit_ternary_asa: `pairwise` must be a named numeric vector")
    miss <- setdiff(c("A1", "A2", "A3"), names(pairwise))
    if (length(miss))
      stop("fit_ternary_asa: `pairwise` is missing param(s): ",
           paste(miss, collapse = ", "))
    a123 <- pairwise[c("A1", "A2", "A3")]
    f2 <- NULL
  }

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

  individual <- if (length(individual_fits) == 0) {
    data.frame(ratio = character(0), C1 = numeric(0), C2 = numeric(0),
               C3 = numeric(0), A4 = numeric(0), n = integer(0),
               stringsAsFactors = FALSE)
  } else {
    do.call(rbind, lapply(individual_fits, function(x)
      data.frame(ratio = x$ratio, C1 = x$C1, C2 = x$C2, C3 = x$C3,
                 A4 = x$A4, n = x$n, stringsAsFactors = FALSE)))
  }
  rownames(individual) <- NULL

  list(reference = reference, response = response,
       base = base_par, pairwise = a123, A4_overall = A4_overall,
       individual = individual,
       fits = list(base = f1, pairwise = f2, overall = f_overall,
                   individual = individual_fits))
}

#' Analyse a ternary mixture with the staged Advanced S/A workflow
#'
#' Fits base curve params (from singles), pairwise interactions A1/A2/A3 (from
#' binaries), an overall three-way A4 (from all data), and an individual A4 for
#' each ternary mixture ratio. Comparing per-ratio A4 against the overall A4
#' reveals interaction that "averages out" in the pooled fit (e.g. synergy at one
#' ratio cancelling antagonism at another).
#' @param df Ternary mixture data: C1, C2, C3 and `Res` (continuous).
#' @param reference "CA" (IA not yet supported).
#' @param response "continuous" (binary not yet validated).
#' @param base Optional named curve params (`max, slope1-3, ec50_1-3`); when
#'   supplied the Stage-1 singles fit is skipped and these are held fixed.
#'   Defaults to `NULL` (fit from singles).
#' @param pairwise Optional named numeric vector (`A1, A2, A3`); when supplied
#'   the Stage-2 binaries fit is skipped and these are held fixed.
#'   Defaults to `NULL` (fit from binaries).
#' @param lower,upper,n_starts,time_limit Forwarded to [fit_model()].
#' @return A list: `reference`, `response`, `base` (named curve params),
#'   `pairwise` (A1,A2,A3), `A4_overall`, `individual` (data frame: ratio,
#'   C1/C2/C3 proportions, A4, n), and `fits` (the underlying [fit_model()]
#'   results for base/pairwise/overall/individual).
#' @export
analyse_ternary <- function(df, reference = "CA",
                            response = c("continuous", "binary"),
                            base = NULL, pairwise = NULL, lower = NULL,
                            upper = NULL, n_starts = 1, time_limit = 30) {
  response <- match.arg(response)
  if (!all(c("C1", "C2", "C3") %in% names(df)))
    stop("analyse_ternary requires C1, C2 and C3 columns")
  fit_ternary_asa(df, reference = reference, response = response, base = base,
                  pairwise = pairwise, lower = lower, upper = upper,
                  n_starts = n_starts, time_limit = time_limit)
}

#' Per-ratio A4 effect-size readout
#'
#' For each ternary ratio, picks the observed ternary point whose reference (CA)
#' prediction is closest to 50% of `max` (the near-EC50 point Sam's workflow
#' targets) and reports the modelled response under three nested models: CA
#' (reference), CA+S/A (pairwise A1/A2/A3, A4 = 0), and CA+S/A+S/A (with the
#' ratio's individual A4). `a4_effect = pred_ASA - pred_SA` measures the size of
#' the effect introduced by that ratio's A4.
#' @param res An [analyse_ternary()] result.
#' @param df The data frame passed to [analyse_ternary()].
#' @return A data frame: ratio, C1, C2, C3 (the chosen point's concentrations),
#'   pred_CA, pred_SA, pred_ASA, a4_effect.
#' @export
ternary_effect_table <- function(res, df) {
  b <- res$base; p <- res$pairwise
  pred1 <- function(c1, c2, c3, A1, A2, A3, A4)
    ca_asa_tri(c1, c2, c3, b[["max"]], b[["slope1"]], b[["slope2"]],
               b[["slope3"]], b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]],
               A1, A2, A3, A4)

  cls <- classify_rows(df)
  tern_idx <- which(cls == "ternary")
  keys <- ternary_ratio_key(df[tern_idx, , drop = FALSE])

  out <- lapply(res$fits$individual, function(x) {
    rows <- tern_idx[keys == x$ratio]
    sub <- df[rows, , drop = FALSE]
    ca <- vapply(seq_len(nrow(sub)), function(i)
      pred1(sub$C1[i], sub$C2[i], sub$C3[i], 0, 0, 0, 0), numeric(1))
    j <- which.min(abs(ca - 0.5 * b[["max"]]))  # near-EC50 point
    c1 <- sub$C1[j]; c2 <- sub$C2[j]; c3 <- sub$C3[j]
    pCA  <- pred1(c1, c2, c3, 0, 0, 0, 0)
    pSA  <- pred1(c1, c2, c3, p[["A1"]], p[["A2"]], p[["A3"]], 0)
    pASA <- pred1(c1, c2, c3, p[["A1"]], p[["A2"]], p[["A3"]], x$A4)
    data.frame(ratio = x$ratio, C1 = c1, C2 = c2, C3 = c3,
               pred_CA = pCA, pred_SA = pSA, pred_ASA = pASA,
               a4_effect = pASA - pSA, stringsAsFactors = FALSE)
  })
  res_df <- if (length(out) == 0) {
    data.frame(ratio = character(0), C1 = numeric(0), C2 = numeric(0),
               C3 = numeric(0), pred_CA = numeric(0), pred_SA = numeric(0),
               pred_ASA = numeric(0), a4_effect = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    do.call(rbind, out)
  }
  rownames(res_df) <- NULL
  res_df
}
