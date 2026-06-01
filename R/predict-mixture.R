#' Evaluate a resolved model's predictor over a concentration list
#'
#' Low-level primitive shared by [fit_model()] and [mixture_predict()]: given the
#' vectorised model function, a named list of concentration columns (`c1`, `c2`,
#' ...), and a full named parameter vector, assemble the by-name call and return
#' the predictions. Keeping this in one place guarantees the fitter (inverse) and
#' the generator (forward) evaluate the identical model.
#' @param fn Vectorised predictor from `model_spec(...)$fn`.
#' @param conc Named list of concentration vectors (`c1`, `c2`, optionally `c3`).
#' @param par Named numeric vector of all model parameters, in `spec$params` order.
#' @return Numeric vector of predictions.
#' @keywords internal
.mixture_eval <- function(fn, conc, par) {
  args <- c(conc, as.list(par))
  names(args) <- c(names(conc), names(par))
  do.call(fn, args)
}

#' Noise-free expected response for a mixture model
#'
#' The forward evaluation shared by the fitter and the synthetic-data generator:
#' given a design frame, a full named parameter vector, and a model
#' (reference + deviation), return the model's expected response for every row.
#' Continuous models return the response value; binary models return a
#' probability in (0, 1). The `n_chem` is inferred from the `C1`/`C2`/`C3`
#' columns present.
#' @param df Data frame with `C1`, `C2` (and `C3` for ternary) columns.
#' @param par Named numeric vector of all parameters
#'   (`model_spec(reference, deviation, n_chem)$params`).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", "DL", or "ASA".
#' @return Numeric vector of expected responses, one per row of `df`.
#' @keywords internal
mixture_predict <- function(df, par, reference, deviation = "reference") {
  conc_cols <- intersect(c("C1", "C2", "C3"), names(df))
  spec <- model_spec(reference, deviation, length(conc_cols))
  missing <- setdiff(spec$params, names(par))
  if (length(missing))
    stop("mixture_predict: `par` is missing parameter(s): ",
         paste(missing, collapse = ", "))
  conc <- stats::setNames(as.list(df[conc_cols]), paste0("c", seq_along(conc_cols)))
  .mixture_eval(spec$fn, conc, par[spec$params])
}
