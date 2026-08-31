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
#' Walks reference -> SA -> DR and DL; a more complex model is accepted only if
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

#' Build the nested-model comparison frame and choose a model
#'
#' Regime-agnostic: works whether the fits are staged (curves fixed, only a/b
#' free) or fully joint (all parameters free). It only reads each fit's
#' `objective` and `df`, so the likelihood-ratio differences are identical in
#' form either way. The parent chain is `SA` vs `reference`, `DR`/`DL` vs `SA`;
#' models absent from `fits` are dropped (e.g. ternary reference + SA only).
#' @param fits Named list of fits (subset of `reference`, `SA`, `DR`, `DL`),
#'   each with `objective` and `df`. A child whose parent is absent from `fits`
#'   is silently dropped from the comparison rather than causing an error.
#' @param n Number of observations.
#' @param response "continuous" or "binary".
#' @param alpha Significance threshold for [select_parsimonious()].
#' @return A list: `comparison` (data frame of LR tests, or `NULL` if no
#'   non-reference model is present) and `chosen` (selected model name).
#' @export
compare_fits <- function(fits, n, response, alpha = 0.05) {
  parent_of <- c(SA = "reference", DR = "SA", DL = "SA")
  parent_of <- parent_of[names(parent_of) %in% names(fits) &
                          parent_of       %in% names(fits)]
  comparison <- if (length(parent_of)) {
    do.call(rbind, lapply(names(parent_of), function(m) {
      p <- parent_of[[m]]
      lrt <- lr_test(fits[[p]]$objective, fits[[m]]$objective,
                     fits[[p]]$df, fits[[m]]$df, n, response)
      data.frame(model = m, parent = p, chi = lrt$chi, df = lrt$df, p = lrt$p)
    }))
  } else {
    NULL
  }
  list(comparison = comparison,
       chosen = select_parsimonious(fits, n, response, alpha))
}
