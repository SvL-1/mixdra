# mixdra Plotting Layer (single + binary) Implementation Plan

**Goal:** Add four interactive plotly plot functions (dose-response curve, observed-vs-predicted, 3-D surface, 2-D isobole) to the `mixdra` package, built on the finished fitting engine.

**Architecture:** Each plot is split into a pure data-builder (numeric output, no plotly, rigorously unit-tested) in `R/plot-data.R` and a thin plotly renderer (smoke-tested) in `R/plot.R`. `fit_model()`/`fit_single()` are enriched with descriptor fields so plots take `(fit, df)` and are self-describing. The 3-D surface and 2-D isobole are ports of Skylar's validated binary script (`recieved/R scripts & excel/R scripts & excel/Example binary 3d curve & 2D isobles.R`), sourcing parameters from the fit object.

**Tech Stack:** R package; `plotly` (added to Suggests, guarded); base `grDevices::contourLines` for isobole paths; `testthat` (edition 3); `roxygen2` for docs.

**Design spec:** `docs/design/specs/2026-05-30-mixdra-binary-plots-design.md`

---

## Running R and tests here (read first)

- **In RStudio:** open `mixdra.Rproj`, `Ctrl+Shift+L` (load_all), then run the test command for each task in the Console.
- **From the CLI:** every non-interactive R call must prepend the user-library path or required packages are not found:
  ```r
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  ```
  e.g. `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="plot-data")'`
- **Test commands** below use `devtools::test(filter = "<name>")`, which runs `tests/testthat/test-<name>.R`. `devtools::load_all()` makes internal (non-exported) functions like `model_spec`, `seed_from_singles`, and the data builders directly callable in tests.
- Existing suite is 11 files / 49 tests, all green — the enrichment in Task 1 must keep it that way.

## File structure

| File | Responsibility |
|---|---|
| `R/fit.R` (modify) | `fit_model()` also records `reference`, `deviation`, `response`, `conc_cols`, `n_chem`, `kind="mixture"` |
| `R/single.R` (modify) | `fit_single()` tags its result `kind="single"` |
| `R/plot-data.R` (create) | Pure data builders + internal helpers (`obs_response`, `predict_grid`, `dr_curve_data`, `obs_pred_data`, `surface_grid_data`, `isobole_data`) |
| `R/plot.R` (create) | `require_plotly` guard + the four exported renderers |
| `tests/testthat/test-fit-enrich.R` (create) | Enrichment descriptor tests |
| `tests/testthat/test-plot-data.R` (create) | Rigorous builder tests (no plotly) |
| `tests/testthat/test-plot.R` (create) | Renderer smoke tests (`skip_if_not_installed("plotly")`) |
| `DESCRIPTION` (modify) | Add `plotly` to Suggests |
| `NAMESPACE` (regenerated) | Export the four renderers |

Note on scope refinement: the design listed `plot_surface(fit, df, axes, n)`. For a **binary** fit there is only one axis pairing (C1↔slope1/ec501, C2↔slope2/ec502); an `axes` swap would mislabel parameters. So the binary signature is `plot_surface(fit, df, n)` with no `axes` arg. Axis selection returns when ternary plotting is implemented (its own plan).

---

## Task 1: Enrich fit objects with model descriptors

**Files:**
- Modify: `R/fit.R` (the final `list(...)` return, around `R/fit.R:117-119`)
- Modify: `R/single.R` (the final `list(...)` return, `R/single.R:39`)
- Test: `tests/testthat/test-fit-enrich.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-fit-enrich.R`:

```r
test_that("fit_model records model descriptors on its result", {
  df <- read.csv(testthat::test_path("fixtures", "binary_mps_cpf_quantal.csv"))
  start <- seed_from_singles(df, "binary")
  set.seed(1)
  fit <- fit_model(df, "CA", "SA", "binary", start = start, n_starts = 1)
  expect_equal(fit$reference, "CA")
  expect_equal(fit$deviation, "SA")
  expect_equal(fit$response, "binary")
  expect_equal(fit$conc_cols, c("C1", "C2"))
  expect_equal(fit$n_chem, 2)
  expect_equal(fit$kind, "mixture")
})

test_that("fit_single tags its result kind", {
  f <- fit_single(c(0, 0.5, 1, 2, 4), c(10, 8, 5, 2, 1))
  expect_equal(f$kind, "single")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "fit-enrich")`
Expected: FAIL — `fit$reference` etc. are `NULL` (`expect_equal(NULL, "CA")` fails); `f$kind` is `NULL`.

- [ ] **Step 3: Implement the enrichment**

In `R/fit.R`, replace the final return list (currently):

```r
  list(par = par, objective = best$value, pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = best$convergence)
```

with:

```r
  list(par = par, objective = best$value, pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = best$convergence,
       reference = reference, deviation = deviation, response = response,
       conc_cols = conc_cols, n_chem = n_chem, kind = "mixture")
```

(`reference`, `deviation`, `response`, `conc_cols`, and `n_chem` are all already in scope inside `fit_model`.)

In `R/single.R`, replace the final return list (currently):

```r
  list(par = res$par, ssr = res$value, convergence = res$convergence)
```

with:

```r
  list(par = res$par, ssr = res$value, convergence = res$convergence,
       kind = "single")
```

- [ ] **Step 4: Run the new test and the full suite**

Run: `devtools::test(filter = "fit-enrich")`
Expected: PASS (2 tests).

Run: `devtools::test()`
Expected: PASS — all existing tests still green (the change is additive).

- [ ] **Step 5: Commit**

```bash
git add R/fit.R R/single.R tests/testthat/test-fit-enrich.R
git commit -m "feat: enrich fit objects with model descriptors for plotting"
```

---

## Task 2: Dose-response curve (`dr_curve_data` + `plot_dose_response`)

**Files:**
- Create: `R/plot-data.R`
- Create: `R/plot.R`
- Test: `tests/testthat/test-plot-data.R`, `tests/testthat/test-plot.R`

- [ ] **Step 1: Write the failing builder test**

Create `tests/testthat/test-plot-data.R`:

```r
test_that("dr_curve_data curve matches ll3_predict for the chosen chemical", {
  fit <- list(kind = "mixture", conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 2, slope2 = 1, ec501 = 0.5, ec502 = 0.3))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), C2 = c(0, 0, 0, 0),
                   Affected = c(10, 6, 3, 1), Exposed = rep(10, 4))
  d <- dr_curve_data(fit, df, chem = 1)
  expect_equal(d$curve$response, ll3_predict(d$curve$conc, 1, 2, 0.5))
  expect_equal(d$observed$response, df$Affected / df$Exposed)
  expect_equal(d$chem, "C1")
})

test_that("dr_curve_data handles a single-chemical fit", {
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0.5, 1, 2, 4), Res = c(8, 5, 2, 1))
  d <- dr_curve_data(fit, df)
  expect_equal(d$curve$response, ll3_predict(d$curve$conc, 10, 2, 1))
  expect_equal(d$observed$response, df$Res)
})
```

- [ ] **Step 2: Run builder test to verify it fails**

Run: `devtools::test(filter = "plot-data")`
Expected: FAIL — `could not find function "dr_curve_data"`.

- [ ] **Step 3: Implement the builder + shared helper**

Create `R/plot-data.R`:

```r
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
```

- [ ] **Step 4: Run builder test to verify it passes**

Run: `devtools::test(filter = "plot-data")`
Expected: PASS (2 tests).

- [ ] **Step 5: Write the failing renderer smoke test**

Create `tests/testthat/test-plot.R`:

```r
test_that("plot_dose_response builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "mixture", conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 2, slope2 = 1, ec501 = 0.5, ec502 = 0.3))
  df <- data.frame(C1 = c(0, 0.5, 1, 2), C2 = 0,
                   Affected = c(10, 6, 3, 1), Exposed = rep(10, 4))
  p <- plot_dose_response(fit, df, chem = 1)
  expect_s3_class(p, "plotly")
})
```

- [ ] **Step 6: Run renderer test to verify it fails**

Run: `devtools::test(filter = "plot")`
Expected: FAIL — `could not find function "plot_dose_response"`.

- [ ] **Step 7: Implement the guard + renderer**

Create `R/plot.R`:

```r
# Thin plotly renderers. Each calls a data builder from R/plot-data.R and
# assembles an interactive plotly object. `plotly` is a Suggested dependency, so
# every renderer first checks it is installed.

#' Error if plotly is not available
#' @keywords internal
require_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("The 'plotly' package is required for plotting. ",
         "Install it with install.packages(\"plotly\").", call. = FALSE)
}

#' Plot a chemical's dose-response curve
#'
#' Draws the fitted log-logistic marginal curve with the observed points
#' overlaid, as an interactive plotly object.
#' @param fit An enriched fit from [fit_model()] or [fit_single()].
#' @param df The data frame the fit was built from.
#' @param chem Index of the chemical whose marginal to draw (mixtures only).
#' @param log_x Use a log10 concentration axis (default `TRUE`). Control points
#'   at concentration 0 are not shown on a log axis.
#' @return A plotly object.
#' @export
plot_dose_response <- function(fit, df, chem = 1, log_x = TRUE) {
  require_plotly()
  d <- dr_curve_data(fit, df, chem)
  p <- plotly::plot_ly()
  p <- plotly::add_markers(p, x = d$observed$conc, y = d$observed$response,
                           name = "observed", marker = list(color = "black", size = 6))
  p <- plotly::add_lines(p, x = d$curve$conc, y = d$curve$response,
                         name = "fitted", line = list(color = "steelblue"))
  plotly::layout(p,
    xaxis = list(title = d$chem, type = if (log_x) "log" else "linear"),
    yaxis = list(title = "response"))
}
```

- [ ] **Step 8: Document and run the renderer test**

Run: `devtools::document()` (regenerates `NAMESPACE` + `man/`)
Run: `devtools::test(filter = "plot")`
Expected: PASS (1 test).

- [ ] **Step 9: Commit**

```bash
git add R/plot-data.R R/plot.R tests/testthat/test-plot-data.R tests/testthat/test-plot.R NAMESPACE man/
git commit -m "feat: add dose-response curve plot (data builder + plotly renderer)"
```

---

## Task 3: Observed-vs-predicted (`obs_pred_data` + `plot_obs_pred`)

**Files:**
- Modify: `R/plot-data.R` (append builder)
- Modify: `R/plot.R` (append renderer)
- Test: `tests/testthat/test-plot-data.R`, `tests/testthat/test-plot.R` (append)

- [ ] **Step 1: Write the failing builder tests**

Append to `tests/testthat/test-plot-data.R`:

```r
test_that("obs_pred_data returns fit$pred for a mixture fit", {
  fit <- list(kind = "mixture", pred = c(0.9, 0.5, 0.2))
  df <- data.frame(C1 = c(0, 1, 2), C2 = 0, Affected = c(9, 5, 2), Exposed = rep(10, 3))
  d <- obs_pred_data(fit, df)
  expect_equal(d$predicted, c(0.9, 0.5, 0.2))
  expect_equal(d$observed, c(0.9, 0.5, 0.2))
})

test_that("obs_pred_data computes predictions for a single fit", {
  fit <- list(kind = "single", par = c(max = 10, slope = 2, ec50 = 1))
  df <- data.frame(C1 = c(0.5, 1, 2), Res = c(8, 5, 2))
  d <- obs_pred_data(fit, df)
  expect_equal(d$predicted, ll3_predict(df$C1, 10, 2, 1))
  expect_equal(d$observed, df$Res)
})
```

- [ ] **Step 2: Run builder tests to verify they fail**

Run: `devtools::test(filter = "plot-data")`
Expected: FAIL — `could not find function "obs_pred_data"`.

- [ ] **Step 3: Implement the builder**

Append to `R/plot-data.R`:

```r
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
```

- [ ] **Step 4: Run builder tests to verify they pass**

Run: `devtools::test(filter = "plot-data")`
Expected: PASS (4 tests total).

- [ ] **Step 5: Write the failing renderer smoke test**

Append to `tests/testthat/test-plot.R`:

```r
test_that("plot_obs_pred builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(kind = "mixture", pred = c(0.9, 0.5, 0.2))
  df <- data.frame(C1 = c(0, 1, 2), C2 = 0, Affected = c(9, 5, 2), Exposed = rep(10, 3))
  p <- plot_obs_pred(fit, df)
  expect_s3_class(p, "plotly")
})
```

- [ ] **Step 6: Run renderer test to verify it fails**

Run: `devtools::test(filter = "plot")`
Expected: FAIL — `could not find function "plot_obs_pred"`.

- [ ] **Step 7: Implement the renderer**

Append to `R/plot.R`:

```r
#' Plot observed vs predicted response
#'
#' Scatter of observed against predicted values with a 1:1 reference line, as a
#' plotly object. Works for any fit.
#' @param fit An enriched fit from [fit_model()] or [fit_single()].
#' @param df The data frame the fit was built from.
#' @return A plotly object.
#' @export
plot_obs_pred <- function(fit, df) {
  require_plotly()
  d <- obs_pred_data(fit, df)
  lim <- range(c(d$observed, d$predicted), na.rm = TRUE)
  p <- plotly::plot_ly()
  p <- plotly::add_markers(p, x = d$observed, y = d$predicted, name = "points",
                           marker = list(color = "black", size = 6))
  p <- plotly::add_lines(p, x = lim, y = lim, name = "1:1",
                         line = list(color = "grey", dash = "dash"))
  plotly::layout(p, xaxis = list(title = "observed"),
                 yaxis = list(title = "predicted"))
}
```

- [ ] **Step 8: Document and run the renderer test**

Run: `devtools::document()`
Run: `devtools::test(filter = "plot")`
Expected: PASS (2 tests total).

- [ ] **Step 9: Commit**

```bash
git add R/plot-data.R R/plot.R tests/testthat/test-plot-data.R tests/testthat/test-plot.R NAMESPACE man/
git commit -m "feat: add observed-vs-predicted plot"
```

---

## Task 4: 3-D surface (`predict_grid` + `surface_grid_data` + `plot_surface`)

**Files:**
- Modify: `R/plot-data.R` (append `predict_grid`, `surface_grid_data`)
- Modify: `R/plot.R` (append renderer)
- Test: `tests/testthat/test-plot-data.R`, `tests/testthat/test-plot.R` (append)

- [ ] **Step 1: Write the failing builder tests**

Append to `tests/testthat/test-plot-data.R`:

```r
test_that("surface_grid_data has correct dims, orientation, and control corner", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 0.9, slope1 = 1.5, slope2 = 1.2, ec501 = 0.1, ec502 = 0.2))
  df <- data.frame(C1 = c(0, 0.1, 0.3), C2 = c(0, 0.2, 0.4),
                   Affected = c(9, 5, 2), Exposed = rep(10, 3))
  g <- surface_grid_data(fit, df, n = 5)
  expect_equal(dim(g$z), c(5, 5))                 # nrow = length(y), ncol = length(x)
  i <- 3; j <- 4                                  # orientation: z[i,j] == f(x[j], y[i])
  expect_equal(g$z[i, j], predict_grid(fit, g$x_vals[j], g$y_vals[i]))
  expect_equal(g$z[1, 1], 0.9)                    # control corner (0,0) == max
})

test_that("surface_grid_data errors on a non-binary fit", {
  fit <- list(n_chem = 1, conc_cols = "C1")
  df <- data.frame(C1 = 1:3, Res = 1:3)
  expect_error(surface_grid_data(fit, df), "binary")
})

test_that("surface_grid_data does not error when some cells are NaN", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 0.5, slope1 = 5, slope2 = 5, ec501 = 0.01, ec502 = 0.01))
  df <- data.frame(C1 = c(0, 1), C2 = c(0, 1), Affected = c(5, 0), Exposed = c(10, 10))
  expect_silent(g <- surface_grid_data(fit, df, n = 10))
  expect_equal(dim(g$z), c(10, 10))
})
```

- [ ] **Step 2: Run builder tests to verify they fail**

Run: `devtools::test(filter = "plot-data")`
Expected: FAIL — `could not find function "surface_grid_data"`.

- [ ] **Step 3: Implement the builders**

Append to `R/plot-data.R`:

```r
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
```

- [ ] **Step 4: Run builder tests to verify they pass**

Run: `devtools::test(filter = "plot-data")`
Expected: PASS (7 tests total).

- [ ] **Step 5: Write the failing renderer smoke test**

Append to `tests/testthat/test-plot.R`:

```r
test_that("plot_surface builds a plotly object", {
  skip_if_not_installed("plotly")
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 0.9, slope1 = 1.5, slope2 = 1.2, ec501 = 0.1, ec502 = 0.2))
  df <- data.frame(C1 = c(0, 0.1, 0.3), C2 = c(0, 0.2, 0.4),
                   Affected = c(9, 5, 2), Exposed = rep(10, 3))
  p <- plot_surface(fit, df, n = 10)
  expect_s3_class(p, "plotly")
})
```

- [ ] **Step 6: Run renderer test to verify it fails**

Run: `devtools::test(filter = "plot")`
Expected: FAIL — `could not find function "plot_surface"`.

- [ ] **Step 7: Implement the renderer**

Append to `R/plot.R`:

```r
#' Plot the fitted 3-D response surface of a binary mixture
#'
#' Observed points (`scatter3d`) overlaid on the fitted response surface, as an
#' interactive plotly object. Binary fits only.
#' @param fit An enriched binary fit from [fit_model()].
#' @param df The data frame the fit was built from.
#' @param n Grid resolution per axis (default 100).
#' @return A plotly object.
#' @export
plot_surface <- function(fit, df, n = 100) {
  require_plotly()
  g <- surface_grid_data(fit, df, n = n)
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, x = g$observed$x, y = g$observed$y, z = g$observed$z,
                         type = "scatter3d", mode = "markers",
                         marker = list(size = 3, color = "blue"), name = "observed")
  p <- plotly::add_surface(p, x = g$x_vals, y = g$y_vals, z = g$z,
                           opacity = 0.8, showscale = FALSE)
  plotly::layout(p, scene = list(
    xaxis = list(title = g$labels$x),
    yaxis = list(title = g$labels$y),
    zaxis = list(title = "response")))
}
```

- [ ] **Step 8: Document and run the renderer test**

Run: `devtools::document()`
Run: `devtools::test(filter = "plot")`
Expected: PASS (3 tests total).

- [ ] **Step 9: Commit**

```bash
git add R/plot-data.R R/plot.R tests/testthat/test-plot-data.R tests/testthat/test-plot.R NAMESPACE man/
git commit -m "feat: add 3-D response surface plot"
```

---

## Task 5: 2-D isobole (`isobole_data` + `plot_isobole`)

**Files:**
- Modify: `R/plot-data.R` (append builder)
- Modify: `R/plot.R` (append renderer)
- Test: `tests/testthat/test-plot-data.R`, `tests/testthat/test-plot.R` (append)

- [ ] **Step 1: Write the failing builder tests**

Append to `tests/testthat/test-plot-data.R`:

```r
test_that("isobole_data returns contour paths at the requested effect levels", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 0.5, 1), C2 = c(0, 0.5, 1),
                   Affected = c(10, 5, 1), Exposed = rep(10, 3))
  d <- isobole_data(fit, df, levels = c(0.25, 0.5, 0.75), n = 40)
  expect_true(all(c("source", "level", "group", "x", "y") %in% names(d)))
  expect_setequal(unique(d$source), "fit")
  expect_setequal(sort(unique(d$level)), sort(c(0.25, 0.5, 0.75) * 1))
})

test_that("isobole_data adds a reference set when reference_fit is given", {
  fit <- list(reference = "CA", deviation = "SA", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5, a = 2))
  ref <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 0.5, 1), C2 = c(0, 0.5, 1),
                   Affected = c(10, 5, 1), Exposed = rep(10, 3))
  d <- isobole_data(fit, df, levels = 0.5, reference_fit = ref, n = 40)
  expect_setequal(unique(d$source), c("fit", "reference"))
})

test_that("isobole_data rejects out-of-range levels", {
  fit <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 1), C2 = c(0, 1), Affected = c(10, 1), Exposed = c(10, 10))
  expect_error(isobole_data(fit, df, levels = c(0.5, 1.5)), "between 0 and 1")
})
```

- [ ] **Step 2: Run builder tests to verify they fail**

Run: `devtools::test(filter = "plot-data")`
Expected: FAIL — `could not find function "isobole_data"`.

- [ ] **Step 3: Implement the builder**

Append to `R/plot-data.R`:

```r
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
```

- [ ] **Step 4: Run builder tests to verify they pass**

Run: `devtools::test(filter = "plot-data")`
Expected: PASS (10 tests total).

- [ ] **Step 5: Write the failing renderer smoke test**

Append to `tests/testthat/test-plot.R`:

```r
test_that("plot_isobole builds a plotly object, including a reference overlay", {
  skip_if_not_installed("plotly")
  fit <- list(reference = "CA", deviation = "SA", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5, a = 2))
  ref <- list(reference = "CA", deviation = "reference", n_chem = 2,
              conc_cols = c("C1", "C2"),
              par = c(max = 1, slope1 = 3, slope2 = 3, ec501 = 0.5, ec502 = 0.5))
  df <- data.frame(C1 = c(0, 0.5, 1), C2 = c(0, 0.5, 1),
                   Affected = c(10, 5, 1), Exposed = rep(10, 3))
  p <- plot_isobole(fit, df, levels = c(0.25, 0.5, 0.75), reference_fit = ref, n = 30)
  expect_s3_class(p, "plotly")
})
```

- [ ] **Step 6: Run renderer test to verify it fails**

Run: `devtools::test(filter = "plot")`
Expected: FAIL — `could not find function "plot_isobole"`.

- [ ] **Step 7: Implement the renderer**

Append to `R/plot.R`:

```r
#' Plot 2-D isoboles (equal-response contours) of a binary mixture
#'
#' Draws contour lines of equal response in the (C1, C2) plane at the requested
#' effect levels (solid, black). If `reference_fit` is supplied, its isoboles are
#' overlaid (dashed, red) so departures from additivity are visible. Binary fits
#' only.
#' @param fit An enriched binary fit from [fit_model()].
#' @param df The data frame the fit was built from.
#' @param levels Effect levels as fractions of `max` (default
#'   `c(0.1, 0.25, 0.5, 0.75, 0.9)`).
#' @param reference_fit Optional enriched reference fit to overlay (dashed).
#' @param n Grid resolution per axis.
#' @return A plotly object.
#' @export
plot_isobole <- function(fit, df, levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
                         reference_fit = NULL, n = 100) {
  require_plotly()
  d <- isobole_data(fit, df, levels = levels, reference_fit = reference_fit, n = n)
  labs <- attr(d, "labels")
  p <- plotly::plot_ly()
  for (grp in unique(d$group)) {
    seg <- d[d$group == grp, ]
    is_ref <- seg$source[1] == "reference"
    p <- plotly::add_lines(p, x = seg$x, y = seg$y, showlegend = FALSE,
                           name = paste0(seg$source[1], " ", signif(seg$level[1], 3)),
                           line = list(color = if (is_ref) "red" else "black",
                                       dash  = if (is_ref) "dash" else "solid"))
  }
  plotly::layout(p, xaxis = list(title = labs$x), yaxis = list(title = labs$y))
}
```

- [ ] **Step 8: Document and run the renderer test**

Run: `devtools::document()`
Run: `devtools::test(filter = "plot")`
Expected: PASS (4 tests total).

- [ ] **Step 9: Commit**

```bash
git add R/plot-data.R R/plot.R tests/testthat/test-plot-data.R tests/testthat/test-plot.R NAMESPACE man/
git commit -m "feat: add 2-D isobole plot with reference overlay"
```

---

## Task 6: Declare the plotly dependency and verify the whole package

**Files:**
- Modify: `DESCRIPTION`
- Test: full suite + `R CMD check`

- [ ] **Step 1: Add plotly to Suggests**

In `DESCRIPTION`, the `Suggests:` field currently reads:

```
Suggests:
    readxl,
    testthat (>= 3.0.0)
```

Change it to:

```
Suggests:
    plotly,
    readxl,
    testthat (>= 3.0.0)
```

- [ ] **Step 2: Regenerate docs and confirm exports**

Run: `devtools::document()`
Then confirm the four renderers are exported. Run:
`Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); cat(readLines("NAMESPACE"), sep="\n")'`
Expected: `NAMESPACE` contains `export(plot_dose_response)`, `export(plot_obs_pred)`, `export(plot_surface)`, `export(plot_isobole)`.

- [ ] **Step 3: Run the full test suite**

Run: `devtools::test()`
Expected: PASS — the original 49 tests plus the new enrichment (2), builder (10), and renderer (4) tests, 0 failures / 0 warnings.

- [ ] **Step 4: Run R CMD check**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); rcmdcheck::rcmdcheck(args = "--no-manual")'`
Expected: 0 errors, 0 warnings. (A NOTE about plotly being a Suggested package used conditionally is acceptable; the renderers guard with `requireNamespace`.)

- [ ] **Step 5: Commit**

```bash
git add DESCRIPTION NAMESPACE man/
git commit -m "build: declare plotly Suggests dependency for the plotting layer"
```

---

## Self-review notes (for the implementer)

- **Orientation is the subtle bit.** `surface_grid_data` returns `z[i, j]` = response at `(x_vals[j], y_vals[i])` (Task 4 pins this). `isobole_data` transposes (`t(z_mat)`) because `contourLines` uses the opposite convention. If the surface or isobole looks rotated/mirrored, this is where to look.
- **Builders never need plotly**; only renderers do. Builder tests must pass even if plotly is not installed. Renderer tests `skip_if_not_installed("plotly")`.
- **Synthetic fits in tests** are hand-built lists (not real `fit_model` output) to keep builder tests fast and deterministic — they exercise the builder math, not the optimiser. Only Task 1 calls the real fitter (with `n_starts = 1` for speed).
- **Out of scope** (later plans): ternary isoplane + ΣTU/z-value plots, the Shiny app, Excel/CSV I/O, README/vignette, static image export.
