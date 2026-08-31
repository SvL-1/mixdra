# Binary tab — interaction-fit workspace (pick model · Autofit/Simulate · Find best): design

**Date:** 2026-06-02
**Status:** ✅ Approved design (pre-implementation). Implemented on branch
`binary-explain-fit` (continuing the legibility work already committed there:
`interaction_help()`, relocated controls, inline explanation).
Builds on the staged Binary tab (`R/app-binary.R`), the curve-fit panel
(`R/app-curve-fit.R`), and the engine (`R/fit.R` `fit_model()`,
`R/analyse.R` `analyse_mixture()`, `R/compare.R` `select_parsimonious()`,
`R/registry.R` `model_spec()`).

---

## 1. Purpose

The current Freeze step is opaque: pressing one button silently freezes the
curves *and* runs a four-model fit-and-select. The result — the interaction
parameters `a`/`b` and which deviation model was chosen — is never something the
user drives or inspects.

This redesign turns Stage 2 into an **interaction workspace that mirrors the
single-curve panel**: the user picks a model, then **Autofits** that model's
`a`/`b`, or **Simulates** it with hand-entered `a`/`b`. The automatic
four-model comparison becomes an explicit, optional **"Find best model"**
action. Nothing computes until the user acts.

This directly answers the two complaints: the step is renamed and re-shaped
around its actual output (`a`/`b`), and it now has the Autofit/Simulate rhythm
the curve panel already uses.

## 2. Decisions (locked)

| Decision | Choice |
|---|---|
| Scope | **Binary tab only.** Single tab unchanged; ternary deferred. |
| Stage 2 shape | A **per-model workspace**: model picker → editable `a`/`b` grid → **Autofit** / **Simulate** buttons, mirroring the curve panel. |
| Model selection | Automatic four-model comparison is demoted to an explicit **"Find best model"** button (today's `analyse_mixture` + `select_parsimonious`). It is the only consumer of `alpha`. |
| Freeze | The Freeze button **only locks the curves and reveals Stage 2** — it no longer fits. Renamed from "Freeze curves → fit interactions" to **"Freeze curves"**. |
| Initial state | Stage 2 **starts empty** (model picked, no fit). The user presses Autofit / Simulate / Find best. |
| `a`/`b` grid | **Same column layout as the curve panel** (Parameter · Meaning · Lower · Upper · Value), but the Lower/Upper cells are **inert "—"** placeholders — `a`/`b` are unconstrained and `fit_model()` rejects bounds on non-curve params. Only **Value** is editable. |
| Reuse | The inline `interaction_help()` block and the relocated `n_starts`/`thorough`/`time_limit`/`alpha` controls already on this branch are kept and repurposed as the Stage 2 fit options. |
| Engine | One small new internal `eval_mixture()` (mirrors `eval_single()`), delegating to `fit_model()`'s existing all-fixed evaluation path. Autofit reuses `fit_model()`; Find best reuses `analyse_mixture()`. |

## 3. Layout

```
Stage 1 · Single curves
  [ chem 1 curve panel ]      [ chem 2 curve panel ]
  [ Freeze curves ]                    ← locks both curves + reveals Stage 2

Stage 2 · Interaction model           (revealed once frozen)
  Model:  [ reference ▾ ]              ← reference / S/A (SA) / dose-ratio (DR) / dose-level (DL)
  <interaction_help block for the selected model + reference>

  Parameters
   Parameter | Meaning                         | Lower | Upper | Value
     a       | overall strength & direction    |  —    |  —    | [   ]
     b       | shift with the mixture ratio    |  —    |  —    | [   ]   ← row shown for DR/DL only
  [ Autofit (a, b) ]   [ Simulate ]
  <fit objective + n readout>

  Fit options:  n_starts [1]   ☐ Thorough fit   time_limit (s) [30]

  ── Find best model ──
  α [0.05]   [ Find best model ]
  [ comparison table: model / parent / chi / df / p — appears after Find best, chosen flagged ]

Stage 3 · Diagnostics                 (revealed once frozen)
  3-D surface | 2-D isobole | observed vs predicted | results table | CIs
  (all reflect the current displayed model + its current a/b)
```

Per selected model the editable params are: **reference** → none (grid hidden,
Autofit just evaluates, Simulate disabled); **S/A** → `a` only; **DR** → `a`,
`b`; **DL** → `a`, `b`.

## 4. Behaviour

Let `curve_params` be the frozen curve vector (`max`, `slope1`, `slope2`,
`ec501`, `ec502`; shared `max` = average of the two single fits — unchanged from
today). `reference` = `input$reference` (CA/IA). `dev` = `input$model` (the
selected deviation). `engine_response` = `"binary"` for quantal, else
`"continuous"`.

- **Autofit (a, b)** — fits *only `dev`'s* interaction params with the curves
  fixed:
  ```
  fit_model(df, reference, dev, engine_response,
            start = curve_params, fixed = names(curve_params),
            n_starts = <thorough ? max(n_starts,20) : n_starts>,
            time_limit = time_limit)
  ```
  On success, write the fitted `a`(/`b`) into the Value inputs and set the
  current fit. (For `dev == "reference"` there are no interaction params, so
  `fit_model` hits its all-fixed branch and simply evaluates.)
- **Simulate** — evaluates `dev` with the **entered** `a`(/`b`), no refit:
  ```
  eval_mixture(df, reference, dev, engine_response,
               curve_params, interaction = c(a = <val_a>, b = <val_b>))
  ```
  Redraws Stage 3 + the objective readout. Requires the needed values to be
  present (`a` for S/A; `a`,`b` for DR/DL) — otherwise a warning notification,
  no-op. Disabled/no-op for `reference` (nothing to simulate).
- **Find best model** — the current automatic path, now explicit:
  ```
  analyse_mixture(df, reference, engine_response, start = curve_params,
                  alpha = alpha, n_starts = <…>, time_limit = time_limit)
  ```
  On success: populate the comparison table, switch the **Model** picker to
  `res$chosen`, and fill that model's `a`/`b` into the Value inputs (so the user
  lands on the chosen model, ready to inspect/Simulate).

The **current displayed fit** (feeding Stage 3 + the objective readout) is the
most recent of Autofit / Simulate / Find-best, for the currently selected
model. Switching the **Model** picker clears the grid to that model's params and
drops the displayed fit back to "none" until the user Autofits/Simulates it (or
it was the Find-best winner).

## 5. Engine: `eval_mixture()`

New internal helper in **`R/fit.R`**, alongside `fit_model()` (it is the
mixture counterpart of `eval_single()`).

```r
#' Forward-evaluate a mixture model at fixed parameters (Simulate)
#'
#' Mirror of [eval_single()] for binary/ternary mixtures: assembles the full
#' parameter vector from the frozen curve params plus caller-supplied
#' interaction params and evaluates the model without optimising, returning the
#' same enriched shape as [fit_model()] so the plotting layer consumes it
#' unchanged. Delegates to [fit_model()]'s all-fixed evaluation path.
#'
#' @param df Mixture data frame (`C1`/`C2`[/`C3`] + response columns).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param curve_params Named curve vector (max, slope1, slope2, ec501, ec502).
#' @param interaction Named numeric of interaction params (e.g. c(a = .., b = ..));
#'   ignored entries that the model does not use are dropped.
#' @return An enriched mixture fit list (as [fit_model()]) with `simulated = TRUE`.
#' @keywords internal
eval_mixture <- function(df, reference, deviation, response,
                         curve_params, interaction = numeric(0)) {
  spec <- model_spec(reference, deviation, length(intersect(c("C1","C2","C3"), names(df))))
  par  <- c(curve_params, interaction[intersect(names(interaction), spec$extra)])
  fit  <- fit_model(df, reference, deviation, response,
                    start = par, fixed = spec$params)   # all fixed -> evaluate only
  fit$simulated <- TRUE
  fit
}
```

`fit_model`'s all-fixed branch (when `free` is empty) already returns
`list(par, objective, pred, residuals, df = 0, n, convergence = 0, reference,
deviation, response, conc_cols, n_chem, kind = "mixture")`, so the existing
`plot_surface`/`plot_isobole`/`plot_obs_pred`/`result_table` consume it
unchanged. `simulated = TRUE` lets the CI panel show a note instead of computing
CIs on hand-entered values.

## 6. Architecture

- **`R/app-binary.R` — UI (`binary_ui`):**
  - Stage 1: rename the button to `"Freeze curves"`; remove the intro sentence
    and the Fit-options group added in the prior commit from Stage 1 (they move
    to Stage 2).
  - Stage 2 (in the existing `conditionalPanel(condition = "output.frozen")`):
    - `selectInput(ns("model"), "Model", choices = c("reference","SA","DR","DL"))`
      with readable labels (`"No interaction (reference)"`, `"Similar action
      (S/A)"`, `"Dose-ratio (DR)"`, `"Dose-level (DL)"`).
    - `uiOutput(ns("interaction_help"))` (kept).
    - An `a`/`b` grid built by a new helper
      `interaction_param_row(ns, param, label, meaning)` — same five-column
      layout as `param_row()` but with static `"—"` in the Lower/Upper columns
      and only a `val_<param>` numeric input. The `b` row is wrapped in
      `conditionalPanel` keyed on the model (`input.model == 'DR' ||
      input.model == 'DL'`); the whole grid is hidden for `reference`.
    - `actionButton(ns("autofit"), "Autofit (a, b)", class="btn-primary")` and
      `actionButton(ns("simulate"), "Simulate")`.
    - `uiOutput(ns("objective"))` (fit objective + n; label "SSR" for
      continuous, "deviance" for quantal).
    - **Fit options** group: `n_starts`, `thorough` (+ `thorough_note`),
      `time_limit` (ids unchanged).
    - **Find best** group: `alpha` + `actionButton(ns("find_best"), "Find best
      model")`, then the `DT::DTOutput(ns("comparison"))`.
  - Stage 3: unchanged outputs (`surface`, `isobole`, `op`, `results`, `cis`).
- **`R/app-binary.R` — server (`binary_server`):**
  - `frozen` gate unchanged; **remove** the `res_r <- eventReactive(frozen(),
    …)` auto-fit. Freeze only locks curves.
  - `current_fit <- reactiveVal(NULL)` — the displayed fit (Autofit/Simulate/
    Find-best result for the current model). Reset to `NULL` when `input$model`,
    `input$reference`, `input$response`, or a curve changes.
  - `last_compare <- reactiveVal(NULL)` — the most recent `analyse_mixture`
    result (for the comparison table); reset on the same invalidations.
  - `observeEvent(input$autofit, …)` → guard `req(frozen(), curve_params())`;
    call `fit_model(...)`; on success update `val_a`/`val_b` and
    `current_fit(fit)`.
  - `observeEvent(input$simulate, …)` → guard frozen + required values present;
    call `eval_mixture(...)`; set `current_fit`.
  - `observeEvent(input$find_best, …)` → guard frozen; call `analyse_mixture`;
    `last_compare(res)`; `updateSelectInput(model, selected = res$chosen)`;
    fill `val_a`/`val_b` from `res$fits[[res$chosen]]`; `current_fit(res$fits[[res$chosen]])`.
  - Reuse `interaction_help` renderUI, now driven by `interaction_help(input$reference, input$model)`.
  - Stage 3 plots/tables gate on `req(frozen(), current_fit())` and use
    `current_fit()`. `comparison` gates on `req(last_compare())`. `cis` shows the
    note when `isTRUE(current_fit()$simulated)`.
- **`R/app-curve-fit.R`:** add `interaction_param_row()` next to `param_row()`
  (it is a sibling presentational helper for the same grid style).
- **`R/fit.R`:** add `eval_mixture()`.

## 7. Data flow

```
Freeze curves ─► frozen(TRUE) ─► curve_params()  (shared-max curve vector)

input$model ─┬─► interaction_help(reference, model)  → Stage 2 explanation
             └─► grid visibility (a; b for DR/DL; none for reference)

Autofit  ─► fit_model(dev, fixed curves) ─► fill val_a/val_b ─► current_fit
Simulate ─► eval_mixture(dev, entered a/b) ───────────────────► current_fit (simulated)
Find best─► analyse_mixture ─► last_compare (table) + pick chosen
            └─► set model picker + val_a/val_b ───────────────► current_fit

current_fit ─► Stage 3 surface/isobole/obs-pred/results/CIs + objective readout
last_compare ─► comparison table
```

## 8. Error handling / edge cases

- **Not frozen:** Stage 2/3 hidden (`conditionalPanel` on `output.frozen`); the
  action observers also `req(frozen())`.
- **Reference model:** grid hidden; Autofit evaluates (objective only, no params
  written); Simulate is a no-op with an explanatory notification.
- **Simulate with blank values:** warning notification ("Enter a (and b) to
  simulate."), no-op — mirrors the curve panel's Simulate guard.
- **Fit failure / time-out:** `tryCatch` → error notification, `current_fit`
  unchanged (mirrors curve panel + today's `res_r`).
- **CIs on a simulated fit:** the `cis` panel shows "CIs apply to fitted
  parameters; this is a simulated set." instead of a table.
- **Switching model after a fit:** `current_fit`/grid reset so Stage 3 never
  shows a model that doesn't match the picker; the user re-Autofits/Simulates.
- **`alpha` only affects Find best**; Autofit/Simulate ignore it.

## 9. Testing

- **`eval_mixture()` unit tests** (new `tests/testthat/test-eval-mixture.R`):
  on a noiseless CA-reference surface,
  `eval_mixture(…, "reference")` gives `objective ≈ 0`; with a known S/A `a`,
  `eval_mixture(…, "SA", interaction = c(a = 5))` reproduces `ca_sa_bi_vec`
  predictions (residuals ≈ 0) and returns `simulated = TRUE`, `kind = "mixture"`,
  the full named `par`, and a `pred` of length `nrow(df)`.
- **`interaction_param_row()`** smoke test: builds a row whose HTML contains the
  `val_<param>` input id and the `"—"` placeholders (no `lo_`/`hi_` inputs).
- **`interaction_help()`**: existing tests stay green (now driven by
  `input$model` directly).
- **`testServer` (binary)** — update the staged-flow test to the new flow:
  freeze (curves only) → `frozen()` TRUE, `current_fit()` NULL; select `model =
  "SA"`, `autofit` → `current_fit()$deviation == "SA"` and `val_a` populated;
  `simulate` with a hand-entered `a` → `current_fit()$simulated` TRUE and the
  objective changes; `find_best` → `last_compare()` has four rows, the picker is
  set to `res$chosen`, and `current_fit()$deviation == chosen`.
- **UI smoke (`binary_ui`)**: Freeze button reads "Freeze curves"; Stage 2 HTML
  contains the model `selectInput`, `binary-val_a`, the Autofit/Simulate/Find
  best buttons, and the relocated `binary-n_starts`/`binary-alpha` ids.

## 10. Out of scope / tracked follow-ups

- **Single tab "About this model" panel** (`model_help_single()`): still a
  separate tracked follow-up; not touched here.
- **Ternary tab**: this workspace is binary-only; the ternary app remains
  deferred. `eval_mixture()` is written generically (n-chem via `model_spec`) so
  a future ternary workspace can reuse it, but no ternary UI is built now.
- **Bounds on `a`/`b`**: deliberately not offered (engine leaves them
  unconstrained); the inert Lower/Upper cells are cosmetic only.
- The **scroll fix** (`page_navbar(fillable = FALSE)`) remains a separate
  working-tree change, unrelated to this spec.
