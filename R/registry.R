#' Look up the specification for a mixture model
#'
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", "DL", or "ASA".
#' @param n_chem Number of chemicals (2 for binary, 3 for ternary).
#' @return A list: `fn` (vectorised predictor), `params` (all free parameter
#'   names), `extra` (deviation parameters beyond the reference), `parent`
#'   (the deviation this one nests within, or NULL for the reference).
#' @keywords internal
model_spec <- function(reference, deviation, n_chem) {
  reference <- match.arg(reference, c("CA", "IA"))
  deviation <- match.arg(deviation, c("reference", "SA", "DR", "DL", "ASA"))
  if (deviation == "ASA" && n_chem != 3)
    stop("model_spec: deviation 'ASA' (Advanced S/A) is ternary-only (n_chem = 3)")
  if (deviation %in% c("DR", "DL") && n_chem != 2)
    stop("model_spec: dose-ratio (DR) and dose-level (DL) deviations are ",
         "binary-only (n_chem = 2); ternary mixtures use 'SA' / 'ASA'")
  if (!n_chem %in% c(2, 3))
    stop("model_spec: only n_chem = 2 or 3 are implemented")

  if (n_chem == 2) {
    base_params <- c("max", "slope1", "slope2", "ec501", "ec502")
    suffix <- "bi"
    extra <- switch(deviation,
                    reference = character(0),
                    SA = "a",
                    DR = c("a", "b"),
                    DL = c("a", "b"))
  } else {
    base_params <- c("max", "slope1", "slope2", "slope3",
                     "ec50_1", "ec50_2", "ec50_3")
    suffix <- "tri"
    # Ternary mixtures support only reference / S/A / Advanced S/A. DR and DL are
    # binary-only (rejected above), so they have no ternary entry here.
    extra <- switch(deviation,
                    reference = character(0),
                    SA = "a",
                    ASA = c("A1", "A2", "A3", "A4"))
  }

  parent <- switch(deviation,
                   reference = NULL,
                   SA = "reference",
                   DR = "SA",
                   DL = "SA",
                   ASA = "SA")
  fn <- if (deviation == "ASA") {
    get(paste0(tolower(reference), "_asa_", suffix, "_vec"), mode = "function")
  } else {
    make_adapter(reference, deviation, n_chem)
  }

  list(fn = fn,
       params = c(base_params, extra),
       extra = extra,
       parent = parent)
}
