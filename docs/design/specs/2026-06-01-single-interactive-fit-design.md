# Design: Interactive Single Chemical tab (Autofit / Simulate + bounds)

**Date:** 2026-06-01
**Status:** Approved (brainstorm)
**Supersedes:** `2026-06-01-single-model-help-design.md` (the collapsed "About this
model" accordion, commit `ace9b8a`). That panel is removed by this redesign.
**Scope:** Single Chemical tab of the mixdra Shiny app, plus the additive engine
changes its features require. Binary tab untouched.

## Problem

The Single Chemical tab is a black box: upload, press **Fit**, read a static
parameter table (`max, slope, ec50, SSR, n`) with no indication of what model
those parameters belong to and no way to interrogate the fit. A researcher
cannot pin a known control `max`, explore how a parameter reshapes the curve, or
sanity-check the optimiser against a hand-chosen set of values.

## Goal

Turn the tab into an explorable fitting tool:

1. Show the fitted model equation in plain sight.
2. Fold the per-parameter explanation into the parameter display.
3. Let the user set lower/upper bounds per parameter (as the Binary tab already
   allows) and edit the parameter values directly.
4. Offer two explicit actions: **Autofit parameters** (optimise) and **Simulate**
   (draw the curve for the current values without optimising), with SSR as live
   feedback.

## The model (ground truth)

`analyse_single(df)` → `fit_single(conc, resp)` fits the three-parameter
log-logistic curve `Y = max / (1 + (C / ec50)^slope)` by least squares.
**Quantal data use the same equation**, fit on the proportion
`Affected / Exposed` (so `max` is then the control proportion, ≈ 1). One equation
header is therefore correct for both response types.

`ll3_predict(conc, max, slope, ec50)` is an exported forward model. The plot
builders (`dr_curve_data`, `obs_pred_data`) require only `fit$kind == "single"`
and `fit$par[["max"/"slope"/"ec50"]]`; single fits have no `fit$pred` (predictions
are recomputed from `par`). So a "simulated" fit is just a fit-shaped list.

## Design

### 1. Engine — `R/single.R`

- Extend `fit_single(conc, resp, lower = NULL, upper = NULL, start = NULL)`:
  - Defaults unchanged: `lower = c(max = 1e-8, slope = 1e-3, ec50 = 1e-8)`,
    `upper = c(max = Inf, slope = 50, ec50 = Inf)`, computed `start`.
  - Supplied `lower`/`upper`/`start` are **named and partial**; each named entry
    overrides the corresponding default. The effective `start` is clamped into
    `[lower, upper]` so L-BFGS-B never receives an out-of-bounds seed.
  - Return shape unchanged: `list(par, ssr, convergence, kind = "single")`.
  - Backward compatible: existing callers (`seed_one`, `fit_curve_from_singles`
    via `seed_one`, `simulate_*` tests) pass nothing and get identical behaviour.
- Add `eval_single(conc, resp, max, slope, ec50)`:
  - No optimisation. Computes `pred <- ll3_predict(...)`,
    `ssr <- sum((resp - pred)^2)`.
  - Returns the **same shape** as `fit_single`:
    `list(par = c(max, slope, ec50), ssr, convergence = NA_integer_, kind = "single")`.
  - `@keywords internal` (only the app consumes it; not exported — YAGNI).

### 2. Engine — `R/analyse.R`

- `analyse_single(df, lower = NULL, upper = NULL, start = NULL)` forwards the
  three optional args to `fit_single`. Quantal column handling unchanged.

### 3. `R/app-io.R`

- Generalise `collect_bounds(values, params = c("max","slope1","slope2","ec501","ec502"))`
  to take the parameter set as an argument. Default = the binary set, so the
  Binary tab's call site is unaffected. The Single tab calls
  `collect_bounds(values, c("max","slope","ec50"))`.

### 4. UI — `R/app-single.R`

- **Remove** `model_help_single()` and its accordion, and the read-only DT
  Parameters table (`output$params`).
- New internal helper `param_row(ns, param, label, meaning, hi_default = NA)`:
  a `fluidRow` with columns **Parameter (label) · Meaning (text) · Lower
  (`numericInput` `lo_<param>`) · Upper (`hi_<param>`) · Value (`val_<param>`)**.
  Mirrors the Binary tab's `bound_row()`.
- Parameters card content, in order:
  - Model header (HTML): `Model: Y = max / (1 + (C / EC50)`<sup>`slope`</sup>`)`.
  - A header row of column titles, then `param_row()` for `max`, `slope`, `ec50`
    (with `slope`'s upper defaulting to 50 to echo the engine default).
  - A live diagnostics line: `SSR` and `n` (`uiOutput`/`textOutput`).
  - Buttons: `actionButton("autofit", "Autofit parameters", class = "btn-primary")`
    and `actionButton("simulate", "Simulate")`.
- The sidebar keeps response type, template download, file upload, and the error
  area; the lone "Fit" button is removed (replaced by the two card buttons).

### 5. Server — `single_server`

- The three `val_*` inputs are the single source of truth. A
  `current_fit <- shiny::reactiveVal(NULL)` drives every output.
- Keep `parsed`, `errs`, `output$errors`, `fit_df` (renames `Conc`→`C1`),
  `output$template`.
- `shiny::observeEvent(input$autofit, …)`:
  - `req(length(errs()) == 0)`.
  - Validate `lower ≤ upper` for each parameter where both are supplied; on a
    violation, `showNotification(..., type = "error")` and return.
  - `b <- collect_bounds(reactiveValuesToList(input), c("max","slope","ec50"))`.
  - `start <- ` the supplied `val_*` values (drop NAs; `NULL` if none).
  - `fit <- tryCatch(analyse_single(fit_df(), lower = b$lower, upper = b$upper,
    start = start), error = notify+NULL)`.
  - On success: `updateNumericInput` each `val_*` to the fitted value, then
    `current_fit(fit)`.
- `shiny::observeEvent(input$simulate, …)`:
  - `req(input$file)` and all three `val_*` present (else notify and return).
  - `resp <- obs_response(fit_df())`.
  - `current_fit(eval_single(fit_df()$C1, resp, val_max, val_slope, val_ec50))`.
- Outputs:
  - `output$dr` / `output$op`: `req(current_fit())`, same `plot_dose_response` /
    `plot_obs_pred` calls and axis labels as today.
  - `output$diagnostics`: `req(current_fit())`, render `SSR` (rounded) and
    `n = nrow(fit_df())`.

### 6. Interaction summary

Autofit optimises and **writes the fitted numbers back into the Value fields**
(user-approved), so the user can then tweak a value and press Simulate to see the
hand-adjusted curve; the SSR line updates either way.

## Testing (filtered runs only — never the full suite)

- `test-single.R`: `fit_single` honours partial `lower`/`upper`/`start` and clamps
  an out-of-bounds start; `eval_single` returns the fit shape with the supplied
  `par` and the correct `ssr` (and `kind == "single"`).
- `test-app-io.R`: `collect_bounds(values, c("max","slope","ec50"))` picks single
  params; the default (binary) param set still works.
- `test-app-modules.R`: rewrite the two `single_server` tests — Autofit
  (`session$setInputs(autofit = 1)`) populates `current_fit()` and recovers the
  known EC50; Simulate with typed `val_*` yields `current_fit()$par` equal to the
  typed values and a finite `ssr`. Replace the `model_help_single` / "About this
  model" tests with a `single_ui` test asserting the model equation and both
  button labels render.

Run with `devtools::test(filter = "single")`, `"app-io"`, `"app-modules"`
individually.

## Out of scope

- Binary tab — unchanged.
- No new exported API beyond the extended `fit_single` / `analyse_single`
  signatures; `eval_single` stays internal.
- No quantal-specific bound enforcement (the user can set `max` upper ≈ 1
  manually).
- No persistence of entered values between sessions.
