#' Approximate 95% confidence intervals for fitted parameters
#'
#' Uses the numerically-estimated Hessian of the objective at the optimum.
#' For continuous data the objective is rescaled to the log-likelihood
#' (sigma^2 = SS/n) so the Hessian yields standard errors; for binary data the
#' deviance Hessian is used directly.
#' @param fit A [fit_model()] result.
#' @param df,reference,deviation,response As passed to [fit_model()].
#' @param level Confidence level.
#' @return A data frame: `parameter`, `estimate`, `lower`, `upper`.
#' @export
param_ci <- function(fit, df, reference, deviation, response, level = 0.95) {
  n_chem <- sum(c("C1", "C2", "C3") %in% names(df))
  spec <- model_spec(reference, deviation, n_chem)
  free <- names(fit$par)
  conc <- df[intersect(c("C1", "C2", "C3"), names(df))]
  names(conc) <- paste0("c", seq_along(conc))

  nll <- function(theta) {
    args <- c(as.list(conc), as.list(stats::setNames(theta, free)))
    names(args) <- c(names(conc), free)
    pred <- do.call(spec$fn, args)
    if (response == "continuous") {
      n <- nrow(df); ss <- obj_ss(df$Res, pred)
      0.5 * n * log(ss / n)              # profile Gaussian negative log-likelihood
    } else {
      0.5 * obj_deviance(df$Exposed, df$Affected,
                         pmin(pmax(pred, 1e-8), 1 - 1e-8))
    }
  }
  H <- numDeriv::hessian(nll, fit$par)
  se <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(free)))
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(parameter = free, estimate = unname(fit$par),
             lower = unname(fit$par) - z * se,
             upper = unname(fit$par) + z * se)
}

#' Assemble the Table-2 style result block from an analysis
#' @param res An [analyse_mixture()] result.
#' @return A matrix: parameters + `objective`/`df` in rows, models in columns.
#' @export
result_table <- function(res) {
  models <- names(res$fits)
  all_params <- unique(unlist(lapply(res$fits, function(f) names(f$par))))
  rows <- c(all_params, "objective", "df")
  tab <- matrix(NA_real_, nrow = length(rows), ncol = length(models),
                dimnames = list(rows, models))
  for (m in models) {
    f <- res$fits[[m]]
    tab[names(f$par), m] <- f$par
    tab["objective", m] <- f$objective
    tab["df", m] <- f$df
  }
  tab
}
