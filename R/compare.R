#' Likelihood-ratio test between a model and its nesting parent
#'
#' @param obj_parent,obj_child Objective values (SS for continuous, deviance for
#'   binary) of the simpler (parent) and more complex (child) models.
#' @param df_parent,df_child Number of free parameters in each model.
#' @param n Number of observations (used for the continuous statistic).
#' @param response "continuous" or "binary".
#' @return A list with `chi`, `df`, and `p` (upper-tail chi-squared p-value).
#' @export
lr_test <- function(obj_parent, obj_child, df_parent, df_child, n, response) {
  ddf <- df_child - df_parent
  chi <- if (response == "continuous") {
    n * log(obj_parent / obj_child)
  } else {
    obj_parent - obj_child
  }
  list(chi = chi, df = ddf, p = stats::pchisq(chi, df = ddf, lower.tail = FALSE))
}

#' Select the most parsimonious model from a set of nested fits
#'
#' Walks reference -> SA -> {DR, DL}; a more complex model is accepted only if
#' it significantly improves on its parent (LR test p < alpha).
#' @param fits Named list of fits (`reference`, `SA`, `DR`, `DL`), each with
#'   `objective` and `df`.
#' @param n Number of observations.
#' @param response "continuous" or "binary".
#' @param alpha Significance threshold.
#' @return The name of the selected model.
#' @export
select_parsimonious <- function(fits, n, response, alpha = 0.05) {
  improves <- function(child, parent) {
    t <- lr_test(fits[[parent]]$objective, fits[[child]]$objective,
                 fits[[parent]]$df, fits[[child]]$df, n, response)
    isTRUE(t$p < alpha)
  }
  best <- "reference"
  if (!is.null(fits$SA) && improves("SA", "reference")) {
    best <- "SA"
    dr_ok <- !is.null(fits$DR) && improves("DR", "SA")
    dl_ok <- !is.null(fits$DL) && improves("DL", "SA")
    if (dr_ok || dl_ok) {
      cand <- c(DR = if (dr_ok) fits$DR$objective else NA,
                DL = if (dl_ok) fits$DL$objective else NA)
      best <- names(which.min(cand))  # smaller objective = better fit
    }
  }
  best
}
