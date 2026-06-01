#' Build a default mixture design grid anchored to the model's EC50s
#'
#' Produces a data frame of concentration columns (`C1`..`Cn`) covering: one
#' control row (all zero); a per-chemical single-compound geometric ladder
#' spanning `ec50 * mult` (others zero) so each marginal curve is identifiable
#' and the staged fit has single-compound rows; and fixed-ratio mixture rays,
#' each a geometric dose series in toxic units. Concentrations are scaled by each
#' chemical's EC50 (read from `par`).
#' @param par Named parameter vector containing the curve parameters; the number
#'   of `slope*` entries determines `n_chem` (2 or 3) and the `ec50*` entries set
#'   the concentration scale.
#' @param reference "CA" or "IA" (only used to resolve EC50 parameter names).
#' @param deviation Deviation name (only used to resolve parameter names).
#' @param ratios Optional list of mixture-fraction vectors (length `n_chem`).
#'   Defaults: binary `list(c(1,1), c(1,3), c(3,1))`; ternary
#'   `list(c(1,1,1), c(1,1,0), c(1,0,1), c(0,1,1))`.
#' @param n_per_ray Number of dose points per mixture ray.
#' @param mult Multipliers for the single-compound ladder (relative to EC50).
#' @return A data frame with columns `C1`..`Cn`.
#' @export
mixture_design <- function(par, reference = "CA", deviation = "reference",
                           ratios = NULL, n_per_ray = 7, mult = 2^(-3:3)) {
  n_chem <- sum(grepl("^slope[0-9]+$", names(par)))
  if (!n_chem %in% c(2, 3))
    stop("mixture_design: `par` must contain 2 or 3 slope* entries; got ", n_chem)
  spec <- model_spec(reference, deviation, n_chem)
  ec_names <- grep("^ec50", spec$params, value = TRUE)   # length n_chem, in order
  ec50 <- as.numeric(par[ec_names])
  if (anyNA(ec50))
    stop("mixture_design: `par` must supply every EC50 (",
         paste(ec_names, collapse = ", "), ")")
  cols <- paste0("C", seq_len(n_chem))

  zero_row <- function() stats::setNames(as.list(rep(0, n_chem)), cols)
  rows <- list(zero_row())                               # control

  for (i in seq_len(n_chem)) {                            # single-compound ladders
    for (cv in ec50[i] * mult) {
      r <- zero_row(); r[[cols[i]]] <- cv
      rows[[length(rows) + 1L]] <- r
    }
  }

  if (is.null(ratios))
    ratios <- if (n_chem == 2) list(c(1, 1), c(1, 3), c(3, 1))
              else list(c(1, 1, 1), c(1, 1, 0), c(1, 0, 1), c(0, 1, 1))

  dose <- 2^seq(-3, 3, length.out = n_per_ray)            # toxic-unit dose series
  for (w in ratios) {
    w <- w / sum(w)                                       # mixture fractions
    for (d in dose) {
      conc <- d * w * ec50                                # conc_i = d * w_i * ec50_i
      rows[[length(rows) + 1L]] <- stats::setNames(as.list(conc), cols)
    }
  }

  out <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(out) <- NULL
  out
}
