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
