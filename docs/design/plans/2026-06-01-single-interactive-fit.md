# Interactive Single Chemical Tab Implementation Plan

**Goal:** Turn the Single Chemical tab into an explorable fitting tool — visible model equation, per-parameter explanation + lower/upper bounds + editable values, and two actions: **Autofit parameters** (optimise) and **Simulate** (draw the curve for the typed values), with a live SSR readout.

**Architecture:** Three small additive engine/IO changes (`fit_single` gains optional bounds/start; new `eval_single` forward evaluator; `analyse_single` forwards the args; `collect_bounds` takes a parameter-set argument), then a rewrite of `R/app-single.R`'s UI + server around a single `current_fit` reactive that both buttons write to. Pure plotting layer is untouched (a simulated fit is a normal single-fit-shaped list).

**Tech Stack:** R package; `shiny` + `bslib` + `plotly` + `DT` (Suggests); `testthat` (edition 3); `roxygen2`.

**Design spec:** `docs/design/specs/2026-06-01-single-interactive-fit-design.md` (supersedes the `2026-06-01-single-model-help-design.md` accordion, commit `ace9b8a`).

---

## Running R and tests here (read first)

- **CLI:** every non-interactive R call MUST prepend the user-library path:
  ```r
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  ```
  e.g. `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="single")'`
- **DO NOT run the full suite** (`devtools::test()`) — it is too slow. Run only the affected file via `devtools::test(filter="<name>")`.
- **Rtools is not installed** — `load_all`/`test`/`document` work; `R CMD check`/install do not. No browser/`shinytest2` tests.
- Commits **must not** add a `Co-Authored-By` trailer.

## File structure

| File | Responsibility |
|---|---|
| `R/single.R` (modify) | `fit_single()` gains optional `lower`/`upper`/`start`; new internal `eval_single()` forward evaluator. |
| `R/analyse.R` (modify) | `analyse_single()` forwards `lower`/`upper`/`start` to `fit_single()`. |
| `R/app-io.R` (modify) | `collect_bounds()` gains a `params` argument (default = binary set). |
| `R/app-single.R` (rewrite) | Remove `model_help_single()` + DT table; add `single_model_equation()` + `param_row()`; rewrite `single_ui()` and `single_server()`. |
| `tests/testthat/test-single.R` (modify) | Tests for `fit_single` bounds/start + clamp, and `eval_single`. |
| `tests/testthat/test-app-io.R` (modify) | Test for `collect_bounds` with the single param set. |
| `tests/testthat/test-app-modules.R` (modify) | Rewrite single_server tests for Autofit/Simulate; replace model-help UI test. |
| `man/*` (generated) | `eval_single.Rd`, `param_row.Rd`, `single_model_equation.Rd` added; `fit_single.Rd`/`analyse_single.Rd`/`collect_bounds.Rd` updated; `model_help_single.Rd` removed. |

**Verified engine facts:**
- `fit_single` currently is `fit_single(conc, resp)` with hardcoded `lower=c(max=1e-8,slope=1e-3,ec50=1e-8)`, `upper=c(max=Inf,slope=50,ec50=Inf)`, computed `start`, returning `list(par, ssr, convergence, kind="single")`.
- Plot builders (`dr_curve_data`, `obs_pred_data` in `R/plot-data.R`) need only `fit$kind=="single"` and `fit$par[["max"/"slope"/"ec50"]]`. `obs_response(df)` (internal, `R/plot-data.R`) gives the observed vector for continuous (`Res`) or quantal (`Affected/Exposed`).
- `ll3_predict(conc, max, slope, ec50)` is exported.

---

## Task 1: `fit_single` bounds/start + `eval_single`

**Files:**
- Modify: `R/single.R`
- Test: `tests/testthat/test-single.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-single.R`:

```r
test_that("fit_single honours a partial upper bound", {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  resp <- ll3_predict(conc, 100, 2, 0.5)   # true slope 2
  fit <- fit_single(conc, resp, upper = c(slope = 1.5))
  expect_lte(fit$par[["slope"]], 1.5 + 1e-6)
})

test_that("fit_single clamps an out-of-bounds start instead of erroring", {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  resp <- ll3_predict(conc, 100, 2, 0.5)
  expect_error(
    fit_single(conc, resp, lower = c(ec50 = 1), start = c(ec50 = 0.001)),
    NA)
})

test_that("eval_single returns the fit shape without optimising", {
  conc <- c(0, 0.5, 1, 2, 4)
  resp <- ll3_predict(conc, 10, 2, 1)
  ev <- eval_single(conc, resp, 10, 2, 1)
  expect_equal(ev$kind, "single")
  expect_equal(unname(ev$par[c("max", "slope", "ec50")]), c(10, 2, 1))
  expect_equal(ev$ssr, 0)
  expect_true(is.na(ev$convergence))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="single")'`
Expected: FAIL — `unused argument (upper = ...)` and `could not find function "eval_single"`.

- [ ] **Step 3: Replace `fit_single` and add `eval_single`**

In `R/single.R`, replace the entire `fit_single` function (currently `fit_single <- function(conc, resp) { ... }`) with:

```r
#' Fit a three-parameter log-logistic curve to single-chemical data
#'
#' Minimises the residual sum of squares. `lower`, `upper`, and `start` are
#' optional **named, partial** overrides of the built-in defaults
#' (`lower = max 1e-8 / slope 1e-3 / ec50 1e-8`, `upper = max Inf / slope 50 /
#' ec50 Inf`, and a data-driven start). The effective start is clamped into
#' `[lower, upper]` so the optimiser never receives an out-of-bounds seed.
#' @param conc Numeric vector of concentrations.
#' @param resp Numeric vector of responses (same length as `conc`).
#' @param lower,upper,start Optional named numeric vectors (`max`/`slope`/`ec50`)
#'   overriding the corresponding defaults.
#' @return A list with `par` (named: max, slope, ec50), `ssr`, `convergence`, and
#'   `kind = "single"`.
#' @export
fit_single <- function(conc, resp, lower = NULL, upper = NULL, start = NULL) {
  stopifnot(length(conc) == length(resp), length(conc) > 3)
  pos <- conc[conc > 0]
  def_start <- c(max = max(resp, na.rm = TRUE), slope = 1, ec50 = stats::median(pos))
  def_lower <- c(max = 1e-8, slope = 1e-3, ec50 = 1e-8)
  def_upper <- c(max = Inf,  slope = 50,   ec50 = Inf)

  lo <- def_lower; if (!is.null(lower)) lo[names(lower)] <- lower
  hi <- def_upper; if (!is.null(upper)) hi[names(upper)] <- upper
  st <- def_start; if (!is.null(start)) st[names(start)] <- start
  st <- pmin(pmax(st, lo), hi)   # keep the seed inside the box

  obj <- function(p) {
    pred <- ll3_predict(conc, p[["max"]], p[["slope"]], p[["ec50"]])
    sum((resp - pred)^2)
  }
  # parscale normalises the three very differently-scaled parameters (max ~ 1e2,
  # slope ~ 1, ec50 ~ 1e-1) so L-BFGS-B's relative tolerance bites uniformly.
  res <- stats::optim(st, obj, method = "L-BFGS-B",
                      lower = lo, upper = hi,
                      control = list(parscale = pmax(abs(st), 1e-8),
                                     factr = 1e-9, maxit = 1000))
  list(par = res$par, ssr = res$value, convergence = res$convergence,
       kind = "single")
}

#' Evaluate the log-logistic curve at given parameters (no optimisation)
#'
#' Forward evaluation used by the Shiny app's "Simulate" action: computes the
#' predicted response and residual sum of squares for caller-supplied parameters,
#' returning the same shape as [fit_single()] so the plotting layer consumes it
#' unchanged.
#' @param conc,resp Concentration and observed-response vectors (same length).
#' @param max,slope,ec50 Curve parameters to evaluate.
#' @return A list with `par` (named max/slope/ec50), `ssr`, `convergence = NA`,
#'   and `kind = "single"`.
#' @keywords internal
eval_single <- function(conc, resp, max, slope, ec50) {
  pred <- ll3_predict(conc, max, slope, ec50)
  list(par = c(max = max, slope = slope, ec50 = ec50),
       ssr = sum((resp - pred)^2),
       convergence = NA_integer_, kind = "single")
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="single")'`
Expected: PASS — the three new tests plus the pre-existing `test-single.R` tests (`ll3_predict`, `fit_single recovers known parameters`) stay green (no-arg calls hit the unchanged defaults).

- [ ] **Step 5: Commit**

```bash
git add R/single.R tests/testthat/test-single.R
git commit -m "feat: fit_single optional bounds/start + eval_single forward evaluator"
```

---

## Task 2: `analyse_single` forwards bounds/start

**Files:**
- Modify: `R/analyse.R:150-153`
- Test: `tests/testthat/test-single.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-single.R`:

```r
test_that("analyse_single forwards bounds to fit_single", {
  d <- data.frame(C1 = c(0, 0.1, 0.3, 1, 3, 10),
                  Res = ll3_predict(c(0, 0.1, 0.3, 1, 3, 10), 100, 2, 0.5))
  fit <- analyse_single(d, upper = c(slope = 1.5))
  expect_lte(fit$par[["slope"]], 1.5 + 1e-6)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="single")'`
Expected: FAIL — `unused argument (upper = c(slope = 1.5))`.

- [ ] **Step 3: Replace `analyse_single`**

In `R/analyse.R`, replace the function (currently lines 150–153) with:

```r
analyse_single <- function(df, lower = NULL, upper = NULL, start = NULL) {
  resp <- if ("Res" %in% names(df)) df$Res else df$Affected / df$Exposed
  fit_single(df$C1, resp, lower = lower, upper = upper, start = start)
}
```

Also update its roxygen block immediately above (currently `#' @return The [fit_single()] result.`) to document the new params — replace the param/return lines with:

```r
#' @param df Data frame with `C1` and `Res` (or `Exposed`/`Affected`).
#' @param lower,upper,start Optional named numeric vectors forwarded to
#'   [fit_single()] (`max`/`slope`/`ec50`).
#' @return The [fit_single()] result.
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="single")'`
Expected: PASS. Also run `devtools::test(filter="nchem")` (it calls `analyse_single(d)` with no extra args) — Expected: PASS (defaults unchanged).

- [ ] **Step 5: Commit**

```bash
git add R/analyse.R tests/testthat/test-single.R
git commit -m "feat: analyse_single forwards optional bounds/start"
```

---

## Task 3: `collect_bounds` parameter-set argument

**Files:**
- Modify: `R/app-io.R:113-132`
- Test: `tests/testthat/test-app-io.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-app-io.R`:

```r
test_that("collect_bounds reads a custom parameter set (single tab)", {
  vals <- list(lo_max = NA, hi_max = NA, lo_slope = NA, hi_slope = 1.5,
               lo_ec50 = 0.01, hi_ec50 = NA)
  b <- collect_bounds(vals, c("max", "slope", "ec50"))
  expect_equal(b$upper, c(slope = 1.5))
  expect_equal(b$lower, c(ec50 = 0.01))
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="app-io")'`
Expected: FAIL — `collect_bounds` reads the hardcoded binary params, so `b$upper`/`b$lower` are `NULL` (slope/ec50 not in the binary set), failing the `expect_equal`.

- [ ] **Step 3: Add the `params` argument**

In `R/app-io.R`, replace the `collect_bounds` function (lines 113–132, including its roxygen block) with:

```r
#' Assemble lower/upper bound vectors from Advanced-panel inputs
#'
#' Reads `lo_<param>` / `hi_<param>` values for `params`; blank/NA entries are
#' dropped. To FIX a parameter, set its lower and upper to the same value. The
#' default `params` is the binary base set; the single tab passes its own.
#' @param values Named list (e.g. a Shiny `input`) holding `lo_*`/`hi_*` numbers.
#' @param params Character vector of parameter names to read.
#' @return A list with `lower` and `upper` named numeric vectors (or NULL).
#' @keywords internal
collect_bounds <- function(values,
                           params = c("max", "slope1", "slope2", "ec501", "ec502")) {
  pick <- function(prefix) {
    v <- vapply(params, function(p) {
      x <- values[[paste0(prefix, p)]]
      if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else as.numeric(x)
    }, numeric(1))
    v <- v[!is.na(v)]
    if (length(v)) v else NULL
  }
  list(lower = pick("lo_"), upper = pick("hi_"))
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="app-io")'`
Expected: PASS — the new test plus the two pre-existing `collect_bounds` tests (which call `collect_bounds(vals)` with the default binary set) stay green.

- [ ] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat: collect_bounds accepts a custom parameter set"
```

---

## Task 4: Rewrite `single_ui` / `single_server` (equation, bounds, Autofit/Simulate)

**Files:**
- Modify (rewrite): `R/app-single.R`
- Test: `tests/testthat/test-app-modules.R`

This task depends on Tasks 1–3. The `param_row` helper mirrors `bound_row()` in `R/app-binary.R`.

- [ ] **Step 1: Write/replace the failing tests**

In `tests/testthat/test-app-modules.R`, **delete** the two blocks added previously:
- `test_that("model_help_single explains the log-logistic model and its parameters", { ... })`
- `test_that("single_ui embeds the model help panel", { ... })`

Also **delete** the existing block `test_that("single_server validates, fits, and exposes a single fit", { ... })` (it drives the removed `input$fit` / `fit_r()`). Keep `test_that("single_server surfaces validation errors and withholds a fit", ...)` and `test_that("single_ui builds a Shiny UI fragment", ...)` unchanged.

Then append these tests:

```r
test_that("single_server autofit fits and exposes the fit via current_fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(
      data.frame(Conc = c(0, 0.1, 0.3, 1, 3, 10),
                 Res  = ll3_predict(c(0, 0.1, 0.3, 1, 3, 10), 100, 2, 0.5)),
      path, row.names = FALSE)
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "single.csv"),
                      val_max = NA, val_slope = NA, val_ec50 = NA,
                      lo_max = NA, hi_max = NA, lo_slope = NA, hi_slope = NA,
                      lo_ec50 = NA, hi_ec50 = NA)
    expect_length(errs(), 0)
    session$setInputs(autofit = 1)
    expect_equal(current_fit()$kind, "single")
    expect_equal(unname(round(current_fit()$par[["ec50"]], 2)), 0.5)
  })
})

test_that("single_server simulate evaluates the typed parameter values", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(
      data.frame(Conc = c(0, 0.1, 0.3, 1, 3, 10),
                 Res  = ll3_predict(c(0, 0.1, 0.3, 1, 3, 10), 100, 2, 0.5)),
      path, row.names = FALSE)
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "single.csv"),
                      val_max = 100, val_slope = 2, val_ec50 = 0.5,
                      simulate = 1)
    expect_equal(unname(current_fit()$par[c("max", "slope", "ec50")]), c(100, 2, 0.5))
    expect_equal(current_fit()$ssr, 0)
  })
})

test_that("single_ui shows the model equation and Autofit/Simulate buttons", {
  html <- as.character(single_ui("single"))
  expect_match(html, "Y = max / (1 + (C / EC50)", fixed = TRUE)
  expect_match(html, "Autofit parameters", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="app-modules")'`
Expected: FAIL — `could not find function "model_help_single"` is gone, but `current_fit` / `input$autofit` not wired and the equation/buttons absent.

- [ ] **Step 3: Replace the whole of `R/app-single.R`**

Overwrite `R/app-single.R` with exactly this content:

```r
# Single Chemical stage: upload one chemical's dose-response data, then either
# Autofit a three-parameter log-logistic curve via analyse_single() or Simulate
# the curve for caller-entered parameter values via eval_single(). Shows the
# curve, observed-vs-predicted, an editable parameter grid with bounds, and a
# live SSR readout.

#' Concentration-axis label from the shared meta store
#'
#' With `chem_field` (e.g. "chem1") uses that chemical's name; without it uses a
#' generic "Concentration". Appends the unit when present. Used by the single
#' stage (generic) and the binary stage (per chemical).
#' @keywords internal
axis_label <- function(meta, chem_field = NULL) {
  nm <- if (!is.null(chem_field)) meta[[chem_field]] else NULL
  unit <- meta$unit
  base <- if (!is.null(nm) && nzchar(nm)) nm else "Concentration"
  if (!is.null(unit) && nzchar(unit)) paste0(base, " (", unit, ")") else base
}

#' Visible model-equation header for the Single Chemical tab
#' @return A `shiny` tag.
#' @keywords internal
single_model_equation <- function() {
  shiny::tags$p(
    shiny::tags$b("Model: "),
    shiny::tags$code("Y = max / (1 + (C / EC50)", shiny::tags$sup("slope"), ")")
  )
}

#' One parameter row: label, plain-English meaning, and lower/upper/value inputs
#' @param ns Module namespace function.
#' @param param Parameter key (`max`/`slope`/`ec50`); drives input ids.
#' @param label Display label.
#' @param meaning One-line explanation.
#' @param hi_default Default for the upper-bound input (NA = blank).
#' @keywords internal
param_row <- function(ns, param, label, meaning, hi_default = NA) {
  shiny::fluidRow(
    shiny::column(2, shiny::tags$b(label)),
    shiny::column(4, shiny::tags$small(meaning)),
    shiny::column(2, shiny::numericInput(ns(paste0("lo_", param)), NULL, value = NA)),
    shiny::column(2, shiny::numericInput(ns(paste0("hi_", param)), NULL, value = hi_default)),
    shiny::column(2, shiny::numericInput(ns(paste0("val_", param)), NULL, value = NA))
  )
}

#' Single Chemical stage UI
#' @param id Module id.
#' @keywords internal
single_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::uiOutput(ns("errors"))
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Dose-response curve"),
                  plotly::plotlyOutput(ns("dr"))),
      bslib::card(bslib::card_header("Observed vs predicted"),
                  plotly::plotlyOutput(ns("op")))
    ),
    bslib::card(
      bslib::card_header("Parameters"),
      single_model_equation(),
      shiny::fluidRow(
        shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
        shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Value")))
      ),
      param_row(ns, "max", "max", "Response at C = 0 (control / upper plateau)."),
      param_row(ns, "slope", "slope", "Steepness of the decline (> 0 = decreasing).",
                hi_default = 50),
      param_row(ns, "ec50", "EC50", "Concentration that halves the response."),
      shiny::uiOutput(ns("diagnostics")),
      shiny::div(
        shiny::actionButton(ns("autofit"), "Autofit parameters", class = "btn-primary"),
        shiny::actionButton(ns("simulate"), "Simulate")
      )
    )
  )
}

#' Single Chemical stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
single_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() paste0("single_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("single", input$response), file, row.names = FALSE)
    )

    parsed <- shiny::reactive({
      shiny::req(input$file)
      read_upload(input$file$datapath)
    })

    errs <- shiny::reactive({
      shiny::req(input$file)
      validate_upload(parsed(), "single", input$response)
    })

    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    fit_df <- shiny::reactive({
      shiny::req(input$file, length(errs()) == 0)
      to_engine_df(parsed(), "single")
    })

    current_fit <- shiny::reactiveVal(NULL)

    # The Value column read as a named numeric (blank -> NA).
    current_values <- function() {
      raw <- list(max = input$val_max, slope = input$val_slope, ec50 = input$val_ec50)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

    shiny::observeEvent(input$autofit, {
      shiny::req(length(errs()) == 0)
      b <- collect_bounds(shiny::reactiveValuesToList(input), c("max", "slope", "ec50"))
      if (!is.null(b$lower) && !is.null(b$upper)) {
        common <- intersect(names(b$lower), names(b$upper))
        if (length(common) && any(b$lower[common] > b$upper[common])) {
          shiny::showNotification("Lower bound exceeds upper bound.", type = "error")
          return()
        }
      }
      vals <- current_values()
      start <- vals[!is.na(vals)]
      if (!length(start)) start <- NULL
      fit <- tryCatch(
        analyse_single(fit_df(), lower = b$lower, upper = b$upper, start = start),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      shiny::updateNumericInput(session, "val_max",   value = round(fit$par[["max"]], 4))
      shiny::updateNumericInput(session, "val_slope", value = round(fit$par[["slope"]], 4))
      shiny::updateNumericInput(session, "val_ec50",  value = round(fit$par[["ec50"]], 4))
      current_fit(fit)
    })

    shiny::observeEvent(input$simulate, {
      shiny::req(input$file, length(errs()) == 0)
      vals <- current_values()
      if (any(is.na(vals))) {
        shiny::showNotification("Enter max, slope and EC50 to simulate.", type = "warning")
        return()
      }
      resp <- obs_response(fit_df())
      current_fit(eval_single(fit_df()$C1, resp,
                              vals[["max"]], vals[["slope"]], vals[["ec50"]]))
    })

    output$dr <- plotly::renderPlotly({
      shiny::req(current_fit())
      p <- plot_dose_response(current_fit(), fit_df())
      plotly::layout(p, xaxis = list(title = axis_label(meta)),
                     yaxis = list(title = if (!is.null(meta$endpoint) && nzchar(meta$endpoint))
                                            meta$endpoint else "Response"))
    })
    output$op <- plotly::renderPlotly({
      shiny::req(current_fit())
      plot_obs_pred(current_fit(), fit_df())
    })
    output$diagnostics <- shiny::renderUI({
      shiny::req(current_fit())
      shiny::tags$p(
        shiny::tags$b("SSR: "), round(current_fit()$ssr, 2),
        "   |   ", shiny::tags$b("n: "), nrow(fit_df())
      )
    })
  })
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="app-modules")'`
Expected: PASS — the rewritten single_server tests, the new single_ui test, and the untouched intro/binary tests all green.

- [ ] **Step 5: Regenerate docs and remove the orphaned man page**

Run: `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::document()'`
Then ensure the removed helper's man page is gone (roxygen may or may not delete it):

```bash
git rm --ignore-unmatch man/model_help_single.Rd
```

Expected: `man/eval_single.Rd`, `man/param_row.Rd`, `man/single_model_equation.Rd` created; `man/fit_single.Rd`, `man/analyse_single.Rd`, `man/collect_bounds.Rd` updated; `NAMESPACE` unchanged (no new exports — `fit_single`/`analyse_single` were already exported).

- [ ] **Step 6: Commit**

```bash
git add R/app-single.R tests/testthat/test-app-modules.R man/ NAMESPACE
git commit -m "feat: interactive Single tab — equation, bounds, Autofit/Simulate"
```

---

## Manual verification (optional, after the plan)

In RStudio after `devtools::load_all()`:

```r
run_app()
```

On the Single Chemical tab: upload `single_cpf_continuous.csv` (Continuous), press
**Autofit parameters** and confirm the curve renders, the Value fields fill with
fitted numbers, and the SSR line shows. Then edit a Value (e.g. halve EC50), press
**Simulate**, and confirm the curve moves and SSR changes without re-optimising.
Set an Upper bound on `slope` (e.g. 1.5) and Autofit to confirm the bound is honoured.

---

## Self-review notes (for the implementer)

- **`obs_response()` is an internal package function** (`R/plot-data.R`); call it
  unqualified from `single_server` — `load_all()` makes it visible.
- **`fit_df` is now a plain `reactive()`** (not `eventReactive`), because both
  buttons and all outputs depend on it; outputs still gate on `current_fit()` so
  nothing draws before an action.
- **Bounds/start are named, partial.** `collect_bounds` drops blanks; an all-blank
  Advanced section yields `NULL` bounds and the engine defaults apply.
- **Autofit overwrites the Value fields** with the fitted numbers (intended), so
  the user can tweak and Simulate next.
- **Do not touch the Binary tab** or `bound_row()` in `R/app-binary.R`.
- **Never run the full test suite** — only `devtools::test(filter=...)` per file.
```
