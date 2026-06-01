# Binary tab — staged, checkpointed UX: design

**Date:** 2026-06-01
**Status:** ✅ Approved design (pre-implementation). To be implemented on a new
branch `binary-staged-ux` off `mixdra-engine`.
Builds on the interactive Single Chemical tab (`R/app-single.R`, commit `1dc7ba0`)
and the staged engine (`docs/superpowers/specs/2026-05-30-staged-fitting-design.md`,
`analyse_mixture()` in `R/analyse.R`).

---

## 1. Purpose

The current Binary tab (`R/app-binary.R`) presents the whole analysis as one
atomic **Fit** button that dumps every output at once (model picker, DR curves,
3-D surface, isobole, obs-vs-pred, results, comparison, CIs). But the fit is
genuinely a **3-stage pipeline**, and that structure is invisible — so the user
can't tell what happened, can't sanity-check the foundation before trusting the
result, and can't see that a binary fit is, at its base, just **two
single-chemical models**.

This redesign makes the pipeline visible and adds one consequential checkpoint,
without turning the tab into a click-heavy wizard.

## 2. Decisions (locked)

| Decision | Choice |
|---|---|
| Overall goal | **Inspect + checkpoint** (not read-only narrative, not full manual wizard) |
| Layout | **Progressive sections on one scrollable page** — stages unlock/expand in order, earlier ones stay visible above |
| Stage 1 interaction | **Fully interactive** — reuse the Single Chemical fitter (params, bounds, Autofit, Simulate), instantiated **twice** |
| 3-D surface | **End only** (Stage 3), as today. The "live baseline surface" idea is **parked** (see §9) |
| Shared `max` | **Average the two per-chemical fits** at freeze (matches current `analyse_mixture` behaviour); UI notes that `max` is shared |
| Reference (CA/IA) + response type | **Up-front controls** in the tab header, as today |

## 3. The three-stage flow

A binary upload contains a single-chemical series for each chemical plus mixture
rows. The tab splits it and walks three sections:

1. **Stage 1 · Single curves.** Split the upload into each chemical's **marginal
   series** (rows where the *other* chemical's concentration is 0) and the
   **mixture rows** (both > 0). Feed each marginal series to an embedded copy of
   the interactive single-chemical fitter. The user reviews/adjusts each curve
   exactly as on the Single tab. Checkpoint button: **"Freeze curves → fit
   interactions"**.
2. **Stage 2 · Interaction models.** On freeze, hold the curve parameters fixed
   and fit only the interaction parameters (`a`, `b`) for all four deviations
   (reference / S/A / DR / DL). Show the comparison table + LR tests with the
   most-parsimonious model flagged. The user **accepts the chosen model or
   overrides** it via a picker.
3. **Stage 3 · Diagnostics.** For the displayed model: 3-D response surface, 2-D
   isobole, observed-vs-predicted, results table, confidence intervals. These are
   today's outputs, relocated to the final section.

Stages 2–3 are hidden/disabled until the Stage 1 checkpoint passes.

## 4. Architecture — extract a reusable curve-fit panel

The core structural change is to **factor the Single Chemical fitter into an
inner, reusable Shiny sub-module** that does not own its data source:

- **New sub-module `curve_fit_panel` (`curve_fit_ui` / `curve_fit_server`).**
  Inputs (injected by the parent): a reactive **fit data frame** (`C1`, response),
  the **response type**, the **meta** store, and an **axis label**. It owns the
  params grid + meaning text + bounds, Autofit (`analyse_single`), Simulate
  (`eval_single`), the DR + obs-vs-pred plots, and the SSR/n diagnostics. It
  **returns a reactive of the current fit** (`par = max/slope/ec50`, `ssr`), so a
  parent can read the fitted curve.
- **Single tab (`single_server`) becomes a thin wrapper:** it keeps the upload,
  template, response-type, and validation controls, derives the fit data frame,
  and mounts one `curve_fit_panel`. No behaviour change for the Single tab.
- **Binary tab (`binary_server`) is restructured** into the three sections,
  mounting **two** `curve_fit_panel`s in Stage 1 (one per marginal series).

This keeps each unit single-purpose and independently testable, and is the reuse
the user explicitly asked for ("reuse the interface we use for the single
chemical … params, autofit, simulate").

## 5. Data flow

```
upload (binary CSV)
  └─ split:
       chem-1 marginal = rows where C2 == 0      → curve_fit_panel #1 → fit1 (max,slope,ec50)
       chem-2 marginal = rows where C1 == 0      → curve_fit_panel #2 → fit2 (max,slope,ec50)
       mixture rows = C1>0 & C2>0  (kept for the engine; full df is passed)
  └─ [Freeze] assemble curve-parameter vector:
       max    = mean(fit1$max, fit2$max)         # §6, decision (a)
       slope1 = fit1$slope ; ec501 = fit1$ec50
       slope2 = fit2$slope ; ec502 = fit2$ec50
  └─ analyse_mixture(full_df, reference, response, start = curve_params,
                     alpha, lower, upper, n_starts, time_limit)
       → list(fits, comparison, chosen, …)        # Stage 2 inputs
  └─ displayed model = override picker value, else $chosen
       → plot_surface / plot_isobole / plot_obs_pred / result_table / param_ci   # Stage 3
```

`analyse_mixture(start = …)` already fixes the supplied curve parameters and fits
only the interaction parameters — so Stage 1 in the UI simply substitutes the
user's reviewed curves for the engine's automatic `fit_curve_from_singles()`.

## 6. Shared `max` reconciliation

The engine uses one shared `max` across both chemicals; two independent single
fits each estimate their own. **Decision (a):** at freeze, set the frozen
`max = mean(fit1$max, fit2$max)` (what `seed_from_singles()` already does). The
two Stage 1 panels still display their own fitted `max`; a short note states that
the mixture model uses a shared `max` (the average). No engine change.

## 7. Checkpoint state & invalidation

- A `frozen` reactive value gates Stages 2–3 (hidden/disabled until set).
- **Freeze** captures the current curve parameters and triggers the
  `analyse_mixture` run (with the existing "Fitting models…" progress).
- Editing either single curve (re-Autofit / Simulate / changed bounds) after
  freezing **invalidates** Stages 2–3: they collapse and require re-freezing, so
  the displayed result can never be stale relative to the shown curves.
- Reference-model / response-type changes in the header also invalidate the
  freeze (they change the fit entirely).

## 8. Error handling & edge cases

- **Validation** reuses `validate_upload(…, "binary", response)`; errors show in
  the header and block freezing.
- **Sparse marginal series** (a chemical with too few distinct points to fit):
  the embedded panel surfaces the single-fit error in place; freezing is blocked
  until both curves have a finite fit. (`analyse_single` already falls back
  gracefully; the panel shows the SSR/n diagnostics so the user sees the problem.)
- **Fit failure** in `analyse_mixture` shows the existing notification and leaves
  Stages 2–3 unpopulated.
- **Advanced bounds** (`n_starts`, `alpha`, `time_limit`, parameter bounds) stay
  in the header "Advanced" accordion (as today). They apply to the interaction
  fit only, since the curve parameters are already frozen — so the bound rows for
  `max`/`slope*`/`ec50*` are dropped from this accordion (those are now set in
  Stage 1), leaving `n_starts`/`alpha`/`time_limit`.

## 9. Out of scope / parked

- **Live baseline surface** (the additive-baseline 3-D surface shown after Stage 1
  and morphing to the fitted surface): deferred. Revisit once the staged flow is
  in use. The data to build it (`mix_response` with `a = 0`) already exists.
- **Ternary tab**: unchanged here; the same `curve_fit_panel` could later back a
  3-panel Stage 1, but that is a separate effort.

## 10. Testing

- `curve_fit_panel` via `shiny::testServer`: injected data → Autofit populates the
  value inputs and the returned fit; Simulate uses entered values; returned-fit
  reactive matches `analyse_single`.
- Binary orchestration via `testServer`: marginal split is correct; freeze gates
  Stage 2 (no `analyse_mixture` call before freeze); frozen `max` is the average;
  override picker switches the displayed model; editing a curve invalidates the
  freeze.
- Single tab regression: behaviour unchanged after extraction (existing
  `test-app-single` / module tests stay green).
