#' Look up the specification for a mixture model
#'
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param n_chem Number of chemicals (2 for binary, 3 for ternary).
#' @return A list: `fn` (vectorised predictor), `params` (all free parameter
#'   names), `extra` (deviation parameters beyond the reference), `parent`
#'   (the deviation this one nests within, or NULL for the reference).
#' @keywords internal
model_spec <- function(reference, deviation, n_chem) {
  reference <- match.arg(reference, c("CA", "IA"))
  deviation <- match.arg(deviation, c("reference", "SA", "DR", "DL"))
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
    # Ternary deviation parameters differ from binary (see model_functions.R):
    # DR carries per-chemical b1/b2/b3 in addition to a; DL carries a/b.
    extra <- switch(deviation,
                    reference = character(0),
                    SA = "a",
                    DR = c("a", "b1", "b2", "b3"),
                    DL = c("a", "b"))
  }

  parent <- switch(deviation,
                   reference = NULL,
                   SA = "reference",
                   DR = "SA",
                   DL = "SA")
  dev_key <- switch(deviation,
                    reference = suffix,
                    SA = paste0("sa_", suffix),
                    DR = paste0("dr_", suffix),
                    DL = paste0("dl_", suffix))
  key <- paste0(tolower(reference), "_", dev_key, "_vec")

  list(fn = get(key, mode = "function"),
       params = c(base_params, extra),
       extra = extra,
       parent = parent)
}
