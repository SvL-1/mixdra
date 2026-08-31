# Binary tab — staged, checkpointed UX: Implementation Plan

**Goal:** Make the Binary tab show its real 3-stage pipeline — fit two single curves → freeze → fit interaction models → diagnostics — by extracting the Single-tab fitter into a reusable Shiny sub-module and mounting it twice.

**Architecture:** Factor the interactive fitter out of `single_server` into a new `curve_fit_panel` sub-module (`curve_fit_ui` / `curve_fit_server`) that does not own its data source. The Single tab becomes a thin wrapper that mounts one panel; the Binary tab mounts two (one per chemical's marginal series), then gates a `frozen` checkpoint that runs `analyse_mixture(start = frozen_curve_params)` and reveals the comparison + diagnostics sections. Two pure helpers (`marginal_df`, `assemble_curve_params`) carry the data-splitting and curve-parameter assembly so the Shiny code stays thin.

**Tech Stack:** R, shiny, bslib, plotly, DT, testthat (`shiny::testServer`), devtools.

**Spec:** `docs/design/specs/2026-06-01-binary-staged-ux-design.md`

---

## File structure

| File | Change | Responsibility |
|---|---|---|
| `R/app-io.R` | Modify (add 2 helpers) | Pure data helpers: `marginal_df()` (split a binary frame into one chemical's single series), `assemble_curve_params()` (combine two single fits into the frozen curve-parameter vector). No Shiny. |
| `R/app-curve-fit.R` | **Create** | The reusable interactive fitter sub-module: `curve_fit_ui()` / `curve_fit_server()`, plus the panel-owned UI helpers `param_row()` and `single_model_equation()` (moved here from `app-single.R`). Owns the params grid, Autofit, Simulate, the DR + obs-vs-pred plots and the SSR/n diagnostics. Returns a reactive of the current fit. |
| `R/app-single.R` | Modify | Becomes a thin wrapper: upload/template/response/validation controls + `fit_df`, mounting one `curve_fit_panel`. `axis_label()` stays here. |
| `R/app-binary.R` | Modify (restructure) | Three progressive sections; mounts two panels in Stage 1, holds the `frozen` checkpoint state, runs `analyse_mixture(start=…)`, reveals Stage 2 (comparison + model picker) and Stage 3 (surface, isobole, obs-vs-pred, results, CIs). `bound_row()` and its Advanced curve-bound inputs are removed. |
| `tests/testthat/test-app-io.R` | Modify (add tests) | Unit tests for `marginal_df` + `assemble_curve_params`. |
| `tests/testthat/test-curve-fit.R` | **Create** | `testServer` tests for `curve_fit_server` (Autofit, Simulate, returned-fit reactive) + UI smoke. |
| `tests/testthat/test-app-modules.R` | Modify | Relocate the single Autofit/Simulate behaviour tests to `test-curve-fit.R`; trim `single_server` test to wiring; rewrite the `binary_server` test for the staged flow. |

All new functions are internal (`@keywords internal`, not exported) — NAMESPACE is unchanged, so no `devtools::document()` run is needed.

**Test runner note:** Per project convention, run only the filtered file under work, never the full suite. Each command below uses `devtools::test(filter='<name>')`, which matches `tests/testthat/test-<name>.R`.

---

## Task 0: Branch setup

**Files:** none (git only)

- [ ] **Step 1: Create the feature branch off `mixdra-engine`**

Run:
```bash
git checkout mixdra-engine
git checkout -b binary-staged-ux
git status
```
Expected: on branch `binary-staged-ux`, working tree clean apart from the untracked `single_cpf_continuous.csv` already present.

---

## Task 1: Pure helper — `marginal_df()`

Splits a binary engine frame into a single chemical's marginal series: keep the rows where the *other* chemical is at 0 (this includes the shared control), drop the other concentration column, and rename this chemical's column to `C1` so the curve-fit panel can consume it like single-chemical data.

**Files:**
- Modify: `R/app-io.R` (add function after `to_engine_df`, before `collect_bounds`)
- Test: `tests/testthat/test-app-io.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-app-io.R`:
```r
test_that("marginal_df extracts a chemical's single series, renaming its conc to C1", {
  df <- data.frame(C1 = c(0, 1, 2, 0, 0, 3),
                   C2 = c(0, 0, 0, 1, 2, 4),
                   Res = c(100, 60, 40, 70, 50, 10))

  m1 <- marginal_df(df, 1)               # rows where C2 == 0
  expect_equal(m1$C1, c(0, 1, 2))
  expect_equal(m1$Res, c(100, 60, 40))
  expect_false("C2" %in% names(m1))

  m2 <- marginal_df(df, 2)               # rows where C1 == 0; C2 renamed to C1
  expect_equal(m2$C1, c(0, 1, 2))
  expect_equal(m2$Res, c(100, 70, 50))
  expect_false("C2" %in% names(m2))
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::test(filter='app-io')"`
Expected: FAIL — `could not find function "marginal_df"`.

- [ ] **Step 3: Write the minimal implementation**

In `R/app-io.R`, add after the `to_engine_df` function (immediately before `collect_bounds`'s roxygen block):
```r
#' One chemical's single-compound series from a binary frame
#'
#' Keeps the rows where the *other* chemical's concentration is 0 (so the shared
#' control row is included), drops the other concentration column, and renames
#' this chemical's concentration column to `C1`. The result has the shape a
#' single-chemical fitter expects (`C1` + response columns).
#' @param df Binary engine data frame (`C1`, `C2`, response columns).
#' @param chem 1 or 2 — which chemical's marginal series to extract.
#' @return A data frame with `C1` and the response columns.
#' @keywords internal
marginal_df <- function(df, chem) {
  this  <- paste0("C", chem)
  other <- paste0("C", if (chem == 1) 2 else 1)
  out <- df[df[[other]] == 0, , drop = FALSE]
  out[[other]] <- NULL
  names(out)[names(out) == this] <- "C1"
  rownames(out) <- NULL
  out
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e "devtools::test(filter='app-io')"`
Expected: PASS — all `test-app-io` tests green.

- [ ] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat: marginal_df helper to split a binary frame into one chemical's single series"
```

---

## Task 2: Pure helper — `assemble_curve_params()`

Combines two single-chemical fits into the frozen curve-parameter vector the engine expects: a shared `max` (the average of the two per-chemical `max`, matching `seed_from_singles()`), plus per-chemical `slope1/slope2` and `ec501/ec502`. Names match `model_spec("CA","reference",2)$params`.

**Files:**
- Modify: `R/app-io.R` (add after `marginal_df`)
- Test: `tests/testthat/test-app-io.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-app-io.R`:
```r
test_that("assemble_curve_params averages max and keeps per-chemical slope/ec50", {
  f1 <- list(par = c(max = 700, slope = 2, ec50 = 1))
  f2 <- list(par = c(max = 600, slope = 1, ec50 = 5))
  p <- assemble_curve_params(f1, f2)
  expect_equal(names(p), c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_equal(unname(p[["max"]]),    650)   # mean(700, 600)
  expect_equal(unname(p[["slope1"]]), 2)
  expect_equal(unname(p[["slope2"]]), 1)
  expect_equal(unname(p[["ec501"]]),  1)
  expect_equal(unname(p[["ec502"]]),  5)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::test(filter='app-io')"`
Expected: FAIL — `could not find function "assemble_curve_params"`.

- [ ] **Step 3: Write the minimal implementation**

In `R/app-io.R`, add immediately after `marginal_df`:
```r
#' Frozen curve-parameter vector from two single-chemical fits
#'
#' Builds the named vector `analyse_mixture(start = …)` holds fixed: a shared
#' `max` (the average of the two per-chemical fits, matching the engine's
#' [seed_from_singles()] behaviour) plus per-chemical `slope1/slope2` and
#' `ec501/ec502`. Names match the binary registry's base parameters.
#' @param fit1,fit2 Single-fit results (each a list with `par = c(max, slope, ec50)`).
#' @return A named numeric vector: `max`, `slope1`, `slope2`, `ec501`, `ec502`.
#' @keywords internal
assemble_curve_params <- function(fit1, fit2) {
  c(max    = mean(c(fit1$par[["max"]], fit2$par[["max"]])),
    slope1 = fit1$par[["slope"]],
    slope2 = fit2$par[["slope"]],
    ec501  = fit1$par[["ec50"]],
    ec502  = fit2$par[["ec50"]])
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e "devtools::test(filter='app-io')"`
Expected: PASS — all `test-app-io` tests green.

- [ ] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat: assemble_curve_params helper to build the frozen curve-parameter vector"
```

---

## Task 3: Extract the reusable `curve_fit_panel` sub-module

Create `R/app-curve-fit.R` with `curve_fit_ui()` / `curve_fit_server()`. The server takes its data as an **injected reactive** (`fit_df`), the shared `meta` store, and an optional `chem_field` for the axis label. It owns the params grid, Autofit, Simulate, the two plots and the diagnostics, and **returns the current-fit reactive**. The panel does NOT know about upload validation — the parent only ever feeds it a valid `fit_df` (its `req()` keeps it dormant otherwise).

The UI helpers `param_row()` and `single_model_equation()` move here from `R/app-single.R` (they belong to the panel). `axis_label()` stays in `R/app-single.R`.

**Files:**
- Create: `R/app-curve-fit.R`
- Modify: `R/app-single.R` (delete `param_row` + `single_model_equation` definitions only — they relocate to the new file; leave `axis_label` and everything else for Task 4)
- Test: `tests/testthat/test-curve-fit.R` (create)

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-curve-fit.R`:
```r
skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

# A clean continuous single-chemical frame with a known curve (ec50 = 0.5).
.curve_df <- function() {
  conc <- c(0, 0.1, 0.3, 1, 3, 10)
  data.frame(C1 = conc, Res = ll3_predict(conc, 100, 2, 0.5))
}

test_that("curve_fit_server Autofit fits and exposes the fit via the returned reactive", {
  meta <- shiny::reactiveValues()
  df <- .curve_df()
  shiny::testServer(curve_fit_server,
                    args = list(fit_df = shiny::reactive(df), meta = meta), {
    session$setInputs(val_max = NA, val_slope = NA, val_ec50 = NA,
                      lo_max = NA, hi_max = NA, lo_slope = NA, hi_slope = NA,
                      lo_ec50 = NA, hi_ec50 = NA, autofit = 1)
    expect_equal(current_fit()$kind, "single")
    expect_equal(unname(round(current_fit()$par[["ec50"]], 2)), 0.5)
  })
})

test_that("curve_fit_server Simulate evaluates the typed parameter values", {
  meta <- shiny::reactiveValues()
  df <- .curve_df()
  shiny::testServer(curve_fit_server,
                    args = list(fit_df = shiny::reactive(df), meta = meta), {
    session$setInputs(val_max = 100, val_slope = 2, val_ec50 = 0.5, simulate = 1)
    expect_equal(unname(current_fit()$par[c("max", "slope", "ec50")]), c(100, 2, 0.5))
    expect_equal(current_fit()$ssr, 0)
  })
})

test_that("curve_fit_ui shows the model equation and Autofit/Simulate buttons", {
  html <- as.character(curve_fit_ui("curve"))
  expect_match(html, "Y = max / (1 + (C / EC50)", fixed = TRUE)
  expect_match(html, "Autofit parameters", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::test(filter='curve-fit')"`
Expected: FAIL — `could not find function "curve_fit_server"` / `"curve_fit_ui"`.

- [ ] **Step 3: Create the sub-module**

Create `R/app-curve-fit.R`:
```r
# Reusable interactive curve-fit panel. Given an injected reactive `fit_df`
# (a single chemical's `C1` + response columns), it fits a 3-parameter
# log-logistic curve via analyse_single() (Autofit) or evaluates caller-entered
# values via eval_single() (Simulate), and shows the curve, observed-vs-predicted,
# an editable parameter grid with bounds, and a live SSR/n readout. It returns a
# reactive of the current fit so a parent (Single or Binary tab) can read it.
# The panel owns no data source and knows nothing about upload validation.

#' Visible model-equation header for the curve-fit panel
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

#' Curve-fit panel UI (plots + parameter grid + Autofit/Simulate)
#' @param id Module id.
#' @keywords internal
curve_fit_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
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

#' Curve-fit panel server
#' @param id Module id.
#' @param fit_df A reactive returning the fit data frame (`C1` + response columns).
#' @param meta Shared reactiveValues for experiment metadata (axis labels).
#' @param chem_field Optional meta field for the x-axis label (e.g. "chem1"); NULL
#'   uses a generic "Concentration" label.
#' @return A reactive returning the current fit (a [fit_single()]/[eval_single()]
#'   result), or NULL before any fit.
#' @keywords internal
curve_fit_server <- function(id, fit_df, meta, chem_field = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    current_fit <- shiny::reactiveVal(NULL)

    # The Value column read as a named numeric (blank -> NA).
    current_values <- function() {
      raw <- list(max = input$val_max, slope = input$val_slope, ec50 = input$val_ec50)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

    shiny::observeEvent(input$autofit, {
      df <- fit_df()
      shiny::req(df)
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
        analyse_single(df, lower = b$lower, upper = b$upper, start = start),
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
      df <- fit_df()
      shiny::req(df)
      vals <- current_values()
      if (any(is.na(vals))) {
        shiny::showNotification("Enter max, slope and EC50 to simulate.", type = "warning")
        return()
      }
      resp <- obs_response(df)
      current_fit(eval_single(df$C1, resp,
                              vals[["max"]], vals[["slope"]], vals[["ec50"]]))
    })

    output$dr <- plotly::renderPlotly({
      shiny::req(current_fit())
      p <- plot_dose_response(current_fit(), fit_df())
      plotly::layout(p, xaxis = list(title = axis_label(meta, chem_field)),
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

    current_fit
  })
}
```

- [ ] **Step 4: Remove the relocated helpers from `R/app-single.R`**

In `R/app-single.R`, delete the `single_model_equation()` definition (its roxygen block + body, lines ~20–28) and the `param_row()` definition (its roxygen block + body, lines ~30–45). They now live in `R/app-curve-fit.R`. Leave `axis_label()` and everything from `single_ui` onward untouched in this task.

- [ ] **Step 5: Run the curve-fit tests to verify they pass**

Run: `Rscript -e "devtools::test(filter='curve-fit')"`
Expected: PASS — 3 tests green (Autofit recovers ec50 ≈ 0.5; Simulate gives ssr 0; UI smoke matches).

- [ ] **Step 6: Confirm no breakage from the helper move**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: PASS — `single_ui`/`binary_ui` still build (they reference the relocated `param_row`/`single_model_equation`, now found in the new file). The single Autofit/Simulate tests still pass here too (they target `single_server`, unchanged in this task).

- [ ] **Step 7: Commit**

```bash
git add R/app-curve-fit.R R/app-single.R tests/testthat/test-curve-fit.R
git commit -m "refactor: extract reusable curve_fit_panel sub-module from single tab"
```

---

## Task 4: Rewire the Single tab to mount the panel

Make `single_server` a thin wrapper that keeps upload/template/response/validation + `fit_df`, then mounts one `curve_fit_panel`. `single_ui` embeds `curve_fit_ui` in the main area. Behaviour is unchanged for the Single tab; the Autofit/Simulate behaviour is now tested in `test-curve-fit.R`, so the duplicate tests are removed from `test-app-modules.R`.

**Files:**
- Modify: `R/app-single.R` (rewrite `single_ui` body's main area + `single_server` body)
- Modify: `tests/testthat/test-app-modules.R` (remove the two relocated single behaviour tests; keep the validation + UI tests)

- [ ] **Step 1: Update the single behaviour tests (remove the relocated ones)**

In `tests/testthat/test-app-modules.R`, **delete** these two test blocks (their behaviour now lives in `test-curve-fit.R`):
- `test_that("single_server autofit fits and exposes the fit via current_fit", { … })`
- `test_that("single_server simulate evaluates the typed parameter values", { … })`

**Keep** unchanged:
- `test_that("single_server surfaces validation errors and withholds a fit", { … })`
- `test_that("single_ui builds a Shiny UI fragment", { … })`
- `test_that("single_ui shows the model equation and Autofit/Simulate buttons", { … })` — this remains valid because `single_ui` embeds `curve_fit_ui`, whose HTML contains the equation and both buttons.

- [ ] **Step 2: Run the trimmed module tests (green baseline before the refactor)**

`single_server`'s rewrite is a behaviour-preserving refactor (the Autofit/Simulate behaviour is already covered by `test-curve-fit.R`), so the kept tests must stay green throughout. Establish the baseline first. Run:
`Rscript -e "devtools::test(filter='app-modules')"`
Expected: PASS for the kept single tests (the old `single_server` still works). This confirms the file is consistent before we change the source. Proceed to the rewrite.

- [ ] **Step 3: Rewrite `single_ui`'s main area to embed the panel**

In `R/app-single.R`, replace the entire `single_ui` function with:
```r
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
    curve_fit_ui(ns("curve"))
  )
}
```

- [ ] **Step 4: Rewrite `single_server` to mount the panel**

In `R/app-single.R`, replace the entire `single_server` function with:
```r
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

    curve_fit_server("curve", fit_df = fit_df, meta = meta)
  })
}
```

- [ ] **Step 5: Run the module tests to verify they pass**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: PASS — validation test still surfaces the "distinct" error; `single_ui` smoke tests still match (equation + buttons present via the embedded panel).

- [ ] **Step 6: Re-run the curve-fit tests (regression)**

Run: `Rscript -e "devtools::test(filter='curve-fit')"`
Expected: PASS — unaffected.

- [ ] **Step 7: Commit**

```bash
git add R/app-single.R tests/testthat/test-app-modules.R
git commit -m "refactor: single tab mounts the curve_fit_panel (behaviour unchanged)"
```

---

## Task 5: Restructure the Binary tab UI into three stages

Rebuild `binary_ui` as: sidebar (response, reference, template, file, thorough, Advanced = `n_starts`/`alpha`/`time_limit` only), then Stage 1 (two embedded panels + Freeze button), then Stage 2 and Stage 3 wrapped in `conditionalPanel`s gated on a `frozen` output flag. The standalone `dr1`/`dr2` outputs are removed (those curves now live in the Stage 1 panels). The `bound_row()` helper and the curve-parameter bound inputs are deleted.

**Files:**
- Modify: `R/app-binary.R` (replace `bound_row` + `binary_ui`)
- Test: `tests/testthat/test-app-modules.R` (the existing `binary_ui` smoke test stays; no new UI assertions needed here — Stage gating is covered by the server test in Task 6)

- [ ] **Step 1: Delete the `bound_row` helper**

In `R/app-binary.R`, delete the `bound_row()` roxygen block + function (lines ~5–13). It is no longer used.

- [ ] **Step 2: Replace `binary_ui`**

In `R/app-binary.R`, replace the entire `binary_ui` function with:
```r
#' Binary Mixture stage UI
#' @param id Module id.
#' @keywords internal
binary_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
      shiny::uiOutput(ns("thorough_note")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced",
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1)
        )
      ),
      shiny::uiOutput(ns("errors"))
    ),

    # Stage 1 -- two single-chemical curve panels + the freeze checkpoint.
    bslib::card(
      bslib::card_header("Stage 1 · Single curves"),
      shiny::p("Review and adjust each chemical's dose-response curve, then freeze ",
               "to fit the interaction. The mixture model uses one shared ",
               shiny::tags$code("max"), " (the average of the two fits)."),
      bslib::layout_columns(
        shiny::div(shiny::h5(shiny::textOutput(ns("chem1_title"))),
                   curve_fit_ui(ns("chem1"))),
        shiny::div(shiny::h5(shiny::textOutput(ns("chem2_title"))),
                   curve_fit_ui(ns("chem2")))
      ),
      shiny::actionButton(ns("freeze"), "Freeze curves → fit interactions",
                          class = "btn-primary"),
      shiny::uiOutput(ns("freeze_note"))
    ),

    # Stage 2 -- interaction-model comparison + model picker (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Interaction models"),
        shiny::selectInput(ns("model"), "Model to display", choices = NULL),
        shiny::helpText("The most parsimonious model is pre-selected; ",
                        "override above to inspect another."),
        DT::DTOutput(ns("comparison"))
      )
    ),

    # Stage 3 -- diagnostics for the displayed model (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 3 · Diagnostics"),
        bslib::layout_columns(
          bslib::card(bslib::card_header("3-D response surface"),
                      plotly::plotlyOutput(ns("surface"))),
          bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                      plotly::plotlyOutput(ns("isobole")))
        ),
        bslib::card(bslib::card_header("Observed vs predicted"),
                    plotly::plotlyOutput(ns("op"))),
        bslib::layout_columns(
          bslib::card(bslib::card_header("Results table"), DT::DTOutput(ns("results"))),
          bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                      DT::DTOutput(ns("cis")))
        )
      )
    )
  )
}
```

- [ ] **Step 3: Verify the UI still builds (server rewrite comes in Task 6)**

The `binary_server` from before still references removed inputs/outputs (`input$fit`, `dr1`, `dr2`, bound inputs), so the server tests will not pass yet — that is expected and handled in Task 6. This step only checks the UI fragment builds. Run:
`Rscript -e "devtools::load_all(); invisible(as.character(binary_ui('binary'))); cat('UI OK\n')"`
Expected: prints `UI OK` with no error.

- [ ] **Step 4: Commit**

```bash
git add R/app-binary.R
git commit -m "feat: stage the binary tab UI (Stage 1 panels, frozen-gated Stages 2-3)"
```

---

## Task 6: Restructure the Binary tab server (staged flow + checkpoint)

Rewrite `binary_server` to: mount two `curve_fit_panel`s on the per-chemical marginal series; expose a `frozen` output flag; assemble `curve_params` from the two fits; gate `analyse_mixture(start = curve_params)` behind the Freeze button; reveal the comparison + diagnostics; let the picker override the displayed model; and invalidate the freeze when a curve or the header model changes.

**Files:**
- Modify: `R/app-binary.R` (replace `binary_server`)
- Test: `tests/testthat/test-app-modules.R` (replace the `binary_server` test)

- [ ] **Step 1: Replace the `binary_server` test with the staged-flow test**

In `tests/testthat/test-app-modules.R`, replace the existing `test_that("binary_server fits the four models and exposes the chosen model", { … })` block with:
```r
test_that("binary_server stages the fit: split, freeze gate, average max, override, invalidate", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      file = list(datapath = csv, name = "binary.csv"))
    expect_length(errs(), 0)

    # marginal split: each panel's series drops the other chemical's column
    expect_false("C2" %in% names(m1()))
    expect_false("C2" %in% names(m2()))

    # nothing frozen yet
    expect_false(isTRUE(frozen()))

    # supply both single curves deterministically via Simulate, with distinct max
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    # frozen curve parameters: shared max is the average; per-chemical slope/ec50 kept
    expect_equal(unname(curve_params()[["max"]]), 650)
    expect_equal(unname(curve_params()[["slope1"]]), 2)
    expect_equal(unname(curve_params()[["ec502"]]), 5)

    # freeze -> fits the four models and reveals the result
    session$setInputs(freeze = 1)
    expect_true(isTRUE(frozen()))
    expect_setequal(names(res_r()$fits), c("reference", "SA", "DR", "DL"))
    expect_equal(shown_fit()$deviation, res_r()$chosen)

    # picker overrides the displayed model
    session$setInputs(model = "DR")
    expect_equal(shown_fit()$deviation, "DR")

    # editing a single curve after freezing invalidates the freeze
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_false(isTRUE(frozen()))
  })
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: FAIL — the current `binary_server` has no `m1`/`m2`/`frozen`/`curve_params`/`freeze` wiring (errors like `could not find function`/`object 'm1' not found`, or `input$fit`-based reactive never produces `frozen`).

- [ ] **Step 3: Replace `binary_server`**

In `R/app-binary.R`, replace the entire `binary_server` function with:
```r
#' Binary Mixture stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
binary_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() paste0("binary_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("binary", input$response), file, row.names = FALSE)
    )

    output$thorough_note <- shiny::renderUI({
      if (isTRUE(input$thorough))
        shiny::div(class = "text-warning",
                   shiny::tags$small("Multi-start fitting may take several minutes."))
    })

    output$chem1_title <- shiny::renderText({
      nm <- meta$chem1; if (!is.null(nm) && nzchar(nm)) nm else "Chemical 1"
    })
    output$chem2_title <- shiny::renderText({
      nm <- meta$chem2; if (!is.null(nm) && nzchar(nm)) nm else "Chemical 2"
    })

    parsed <- shiny::reactive({
      shiny::req(input$file); read_upload(input$file$datapath)
    })
    errs <- shiny::reactive({
      shiny::req(input$file); validate_upload(parsed(), "binary", input$response)
    })
    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    engine_df <- shiny::reactive({
      shiny::req(input$file, length(errs()) == 0)
      to_engine_df(parsed(), "binary")
    })
    m1 <- shiny::reactive(marginal_df(engine_df(), 1))
    m2 <- shiny::reactive(marginal_df(engine_df(), 2))

    # Stage 1: two embedded single-chemical fitters, one per marginal series.
    fit1 <- curve_fit_server("chem1", fit_df = m1, meta = meta, chem_field = "chem1")
    fit2 <- curve_fit_server("chem2", fit_df = m2, meta = meta, chem_field = "chem2")

    # Checkpoint state. `frozen` gates Stages 2-3 (exposed to the UI as an output).
    frozen <- shiny::reactiveVal(FALSE)
    output$frozen <- shiny::reactive(isTRUE(frozen()))
    shiny::outputOptions(output, "frozen", suspendWhenHidden = FALSE)

    # Frozen curve-parameter vector (shared max = average of the two fits).
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2())
      assemble_curve_params(fit1(), fit2())
    })

    output$freeze_note <- shiny::renderUI({
      if (is.null(fit1()) || is.null(fit2()))
        shiny::div(class = "text-muted",
                   shiny::tags$small(
                     "Fit both single curves (Autofit or Simulate) before freezing."))
    })

    # Freeze requires both curves; marking frozen triggers the interaction fit.
    shiny::observeEvent(input$freeze, {
      if (is.null(fit1()) || is.null(fit2())) {
        shiny::showNotification("Fit both single curves before freezing.", type = "warning")
        return()
      }
      frozen(TRUE)
    })

    # A changed curve or header model invalidates the freeze (result can't go stale).
    shiny::observeEvent(
      list(fit1(), fit2(), input$reference, input$response),
      { if (isTRUE(frozen())) frozen(FALSE) },
      ignoreInit = TRUE)

    # Stage 2/3: fit reference + deviations with the frozen curve params fixed.
    res_r <- shiny::eventReactive(input$freeze, {
      shiny::req(fit1(), fit2())
      df <- engine_df()
      engine_response <- if (input$response == "quantal") "binary" else "continuous"
      n_starts <- if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts
      start <- curve_params()
      shiny::withProgress(message = "Fitting interaction models...", value = 0.5, {
        tryCatch(
          analyse_mixture(df, reference = input$reference, response = engine_response,
                          start = start, alpha = input$alpha,
                          n_starts = n_starts, time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
            NULL
          })
      })
    })

    # Refresh the model picker after each fit, defaulting to the chosen model.
    shiny::observeEvent(res_r(), {
      shiny::req(res_r())
      shiny::updateSelectInput(session, "model",
                               choices = names(res_r()$fits), selected = res_r()$chosen)
    })

    # Fall back to the chosen model until the picker's input has populated.
    shown_fit <- shiny::reactive({
      shiny::req(res_r())
      m <- input$model
      if (is.null(m) || !m %in% names(res_r()$fits)) m <- res_r()$chosen
      res_r()$fits[[m]]
    })

    output$surface <- plotly::renderPlotly({
      shiny::req(frozen()); plot_surface(shown_fit(), engine_df())
    })
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen())
      plot_isobole(shown_fit(), engine_df(), reference_fit = res_r()$fits$reference)
    })
    output$op <- plotly::renderPlotly({
      shiny::req(frozen()); plot_obs_pred(shown_fit(), engine_df())
    })

    output$results <- DT::renderDT({
      shiny::req(frozen(), res_r())
      tab <- round(result_table(res_r()), 4)
      DT::datatable(as.data.frame(tab), options = list(dom = "t"))
    })
    output$comparison <- DT::renderDT({
      shiny::req(frozen(), res_r())
      DT::datatable(res_r()$comparison, rownames = FALSE, options = list(dom = "t"))
    })
    output$cis <- DT::renderDT({
      shiny::req(frozen(), shown_fit())
      f <- shown_fit()
      DT::datatable(param_ci(f, engine_df(), f$reference, f$deviation, f$response),
                    rownames = FALSE, options = list(dom = "t"))
    })

    res_r   # return for testability
  })
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: PASS — split correct; `frozen` false before freeze; `curve_params()[["max"]]` == 650; freeze yields four fits; picker override gives `DR`; re-simulating chem1 flips `frozen` back to false.

- [ ] **Step 5: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat: staged binary server with frozen checkpoint and interaction fit"
```

---

## Task 7: Integration verification

Confirm the three touched test files are green together and the app still launches end-to-end, then commit any final touch-ups.

**Files:** none (verification only)

- [ ] **Step 1: Run all touched test files**

Run:
```bash
Rscript -e "devtools::test(filter='app-io')"
Rscript -e "devtools::test(filter='curve-fit')"
Rscript -e "devtools::test(filter='app-modules')"
```
Expected: all PASS, 0 failures.

- [ ] **Step 2: Confirm the package loads cleanly**

Run: `Rscript -e "devtools::load_all(); cat('load OK\n')"`
Expected: prints `load OK` with no errors or warnings about undefined functions (`marginal_df`, `assemble_curve_params`, `curve_fit_ui`, `curve_fit_server` all resolve).

- [ ] **Step 3: Manual app smoke test**

The Stage gating uses `conditionalPanel`, which only evaluates in a live browser session (not under `testServer`). Verify it manually once. Launch the app (the run helper is in `R/app-run.R`):
```bash
Rscript -e "devtools::load_all(); run_app()"
```
> If the launcher is named differently, find it with `Grep` for `shinyApp` in `R/app-run.R`.

In the **Binary Mixture** tab: download the binary template, re-upload it, Autofit both Stage 1 curves, click **Freeze curves → fit interactions**, and confirm Stage 2 (comparison + picker) and Stage 3 (surface, isobole, obs-vs-pred, results, CIs) appear. Then re-Autofit one curve and confirm Stages 2–3 collapse. Close the app.

- [ ] **Step 4: Final commit (if anything changed during smoke test)**

Only if you adjusted source during the smoke test:
```bash
git add -A
git commit -m "fix: binary staged-ux smoke-test adjustments"
```
Otherwise skip.

---

## Self-review notes (coverage against the spec)

- **§2 decisions** — Inspect+checkpoint flow, progressive one-page layout (Stage 1 always visible; Stages 2–3 in `conditionalPanel`), Stage 1 fully interactive via the embedded panel (Task 3–6), 3-D surface end-only (Stage 3, Task 5–6), shared `max` = average (Task 2 + Task 6 `curve_params`), reference/response up-front in header (Task 5). ✓
- **§3 three-stage flow** — marginal split (`marginal_df`, Task 1; `m1`/`m2`, Task 6), Freeze checkpoint (Task 6), comparison + override picker (Task 5–6), Stage 3 diagnostics (Task 5–6). ✓
- **§4 reusable panel** — `curve_fit_ui`/`curve_fit_server` returning the fit reactive; Single tab thin wrapper; Binary mounts two (Tasks 3, 4, 6). ✓
- **§5 data flow** — `analyse_mixture(start = curve_params)` with the frozen vector (Task 6). ✓
- **§6 shared max** — averaged in `assemble_curve_params`; UI note in Stage 1 header (Task 2, Task 5). ✓
- **§7 checkpoint/invalidation** — `frozen` reactiveVal gates Stages 2–3; invalidation observer on `fit1/fit2/reference/response` (Task 6). ✓
- **§8 errors/edge cases** — validation reused (`errs`), sparse series surfaced by the panel's Autofit `tryCatch`, mixture-fit failure notification, Advanced reduced to `n_starts`/`alpha`/`time_limit` with curve bound rows dropped (Tasks 5, 6). ✓
- **§9 parked** — live baseline surface intentionally not implemented. ✓
- **§10 testing** — panel `testServer` tests (Task 3), binary orchestration `testServer` test (Task 6), single-tab regression (Task 4). ✓
