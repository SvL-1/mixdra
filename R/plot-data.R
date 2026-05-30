# Pure data builders for the plotting layer. No plotly here, so these are fully
# unit-testable. The renderers in R/plot.R consume their output.

#' Observed response vector from a mixture/single data frame
#'
#' Continuous data carry `Res`; binary (quantal) data carry `Affected`/`Exposed`
#' and the modelled response is the proportion `Affected / Exposed`.
#' @keywords internal
obs_response <- function(df) {
  if ("Res" %in% names(df)) df$Res else df$Affected / df$Exposed
}

#' Build the data for one chemical's (marginal) dose-response curve
#'
#' When every other chemical is at 0, both CA and IA reduce to the single
#' three-parameter log-logistic, so the marginal curve is just [ll3_predict()]
#' with that chemical's `slope`/`ec50` and the shared `max`.
#' @param fit An enriched fit (from [fit_model()] or [fit_single()]).
#' @param df The data frame the fit was built from.
#' @param chem Index of the chemical whose marginal to draw (mixtures only).
#' @return A list: `curve` (data frame conc/response), `observed` (data frame
#'   conc/response), `chem` (the concentration column name).
#' @keywords internal
dr_curve_data <- function(fit, df, chem = 1) {
  if (isTRUE(fit$kind == "single")) {
    conc_col <- "C1"
    mx <- fit$par[["max"]]; sl <- fit$par[["slope"]]; ec <- fit$par[["ec50"]]
    keep <- rep(TRUE, nrow(df))
  } else {
    conc_cols <- fit$conc_cols
    if (chem < 1 || chem > length(conc_cols))
      stop("`chem` must be between 1 and ", length(conc_cols), call. = FALSE)
    conc_col <- conc_cols[chem]
    others <- setdiff(conc_cols, conc_col)
    keep <- if (length(others) == 0) rep(TRUE, nrow(df))
            else rowSums(df[others] == 0) == length(others)
    mx <- fit$par[["max"]]
    sl <- fit$par[[paste0("slope", chem)]]
    ec <- fit$par[[paste0("ec50", chem)]]
  }
  conc <- df[[conc_col]][keep]
  resp <- obs_response(df)[keep]
  pos  <- conc[conc > 0]
  grid <- seq(min(pos), max(conc), length.out = 200)
  list(curve    = data.frame(conc = grid, response = ll3_predict(grid, mx, sl, ec)),
       observed = data.frame(conc = conc, response = resp),
       chem     = conc_col)
}

#' Observed vs predicted response for any fit
#'
#' Mixture fits carry their fitted values in `fit$pred`; single-chemical fits do
#' not, so predictions are recomputed from the log-logistic parameters.
#' @param fit An enriched fit from [fit_model()] or [fit_single()].
#' @param df The data frame the fit was built from.
#' @return A data frame with `observed` and `predicted` columns.
#' @keywords internal
obs_pred_data <- function(fit, df) {
  predicted <- if (!is.null(fit$pred)) fit$pred
               else ll3_predict(df$C1, fit$par[["max"]], fit$par[["slope"]],
                                fit$par[["ec50"]])
  data.frame(observed = obs_response(df), predicted = predicted)
}

#' Evaluate a fitted binary model over c1/c2 vectors
#'
#' Looks up the vectorised predictor for the fit's model
#' (`model_spec(reference, deviation, 2)$fn`) and supplies the fitted parameters
#' by name, so deviation parameters (`a`, `b`) are included automatically.
#' @keywords internal
predict_grid <- function(fit, c1, c2) {
  spec <- model_spec(fit$reference, fit$deviation, 2)
  args <- c(list(c1 = c1, c2 = c2), as.list(fit$par[spec$params]))
  do.call(spec$fn, args)
}

#' Build the response surface grid for a binary fit
#'
#' Mirrors Skylar's construction so the surface is not transposed:
#' `expand.grid(C1 = x_vals, C2 = y_vals)` (C1/x varies fastest), then a matrix
#' with `nrow = length(y_vals)`, `byrow = TRUE`, so `z[i, j]` is the response at
#' `x_vals[j], y_vals[i]` — exactly what `plotly::add_surface(x, y, z)` expects.
#' @param fit An enriched binary fit from [fit_model()].
#' @param df The data frame the fit was built from.
#' @param n Grid resolution per axis.
#' @return A list: `x_vals`, `y_vals`, `z` (matrix), `observed` (data frame
#'   x/y/z), `labels` (list of axis names).
#' @keywords internal
surface_grid_data <- function(fit, df, n = 100) {
  if (!isTRUE(fit$n_chem == 2))
    stop("surface/isobole require a binary (2-chemical) fit", call. = FALSE)
  cols <- fit$conc_cols
  x_vals <- seq(min(df[[cols[1]]]), max(df[[cols[1]]]), length.out = n)
  y_vals <- seq(min(df[[cols[2]]]), max(df[[cols[2]]]), length.out = n)
  grid <- expand.grid(C1 = x_vals, C2 = y_vals)   # C1 (x) varies fastest
  z <- predict_grid(fit, grid$C1, grid$C2)
  z_mat <- matrix(z, nrow = length(y_vals), ncol = length(x_vals), byrow = TRUE)
  list(x_vals = x_vals, y_vals = y_vals, z = z_mat,
       observed = data.frame(x = df[[cols[1]]], y = df[[cols[2]]],
                             z = obs_response(df)),
       labels = list(x = cols[1], y = cols[2]))
}

#' Build isobole (equal-response) contour paths for a binary fit
#'
#' Reuses [surface_grid_data()] and extracts contour lines with
#' [grDevices::contourLines()] at `levels * max` (the effect levels as fractions
#' of the control response). Optionally overlays a reference model's isoboles
#' (Skylar's deviation-vs-additivity comparison), using the reference fit's own
#' `max`. `contourLines` expects `z[i, j]` at `x[i], y[j]`, whereas the surface
#' grid is `z[y, x]`, so the matrix is transposed before extraction.
#' @param fit An enriched binary fit from [fit_model()].
#' @param df The data frame the fit was built from.
#' @param levels Effect levels as fractions of `max`, strictly in (0, 1).
#' @param reference_fit Optional enriched reference fit to overlay.
#' @param n Grid resolution per axis.
#' @return A data frame with columns `source` ("fit"/"reference"), `level`
#'   (absolute response), `group` (path id), `x`, `y`. Attribute `"labels"`
#'   holds the axis names.
#' @keywords internal
isobole_data <- function(fit, df, levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
                         reference_fit = NULL, n = 100) {
  if (any(levels <= 0 | levels >= 1))
    stop("`levels` must be fractions strictly between 0 and 1", call. = FALSE)
  g <- surface_grid_data(fit, df, n = n)
  contour_df <- function(z_mat, src, mx_src) {
    cl <- grDevices::contourLines(x = g$x_vals, y = g$y_vals, z = t(z_mat),
                                  levels = levels * mx_src)
    if (length(cl) == 0) return(NULL)
    do.call(rbind, lapply(seq_along(cl), function(k) {
      data.frame(source = src, level = cl[[k]]$level,
                 group = paste(src, k, sep = "_"),
                 x = cl[[k]]$x, y = cl[[k]]$y)
    }))
  }
  out <- contour_df(g$z, "fit", fit$par[["max"]])
  if (!is.null(reference_fit)) {
    gr <- surface_grid_data(reference_fit, df, n = n)
    out <- rbind(out, contour_df(gr$z, "reference", reference_fit$par[["max"]]))
  }
  attr(out, "labels") <- g$labels
  out
}
