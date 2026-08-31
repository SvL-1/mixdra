# Binary tab: per-model joint refine (post-selection parameter polish)

**Date:** 2026-06-03
**Status:** Design proposed; pending user review → implementation plan
**Builds on:** `2026-06-02-binary-model-selection-loop-ux-design.md` (the 3-stage
staged compare-all loop) and the `excel-joint-fitting` / `mixdra-staged-fitting`
memories.

## Context

The binary tab selects the interaction model with the **staged** method (curves
fixed from the single compounds, only `a`/`b` fitted per model), and Sam signed
off that selection **must** stay staged — a joint curve+interaction fit can
re-absorb a real interaction into the curves and mask it (the quantal
`reference→DR` flip). The live "Fit & compare all models" loop (Stage 2) and its
results table embody that staged verdict.

Separately, the **Excel workbook fits each model jointly** (all params at once,
one SSR), and the team wants those Excel-comparable parameter values back for the
model they care about. That joint fit existed previously as the "Optimize all
params (joint)" card (`refine_joint()`), removed from the UI in the 3-stage
redesign. The engine function `refine_joint()` still exists.

The key realization: a joint fit is safe **as a post-selection polish of one
model's reported parameters**, but unsafe **as the selection mechanism**. This
feature re-adds the former while protecting the latter.

## Goal

Let the user refine the parameters of a **selected** interaction model — either by
**optimizing jointly** (all params at once, seeded from the staged fit, Excel-style)
or by **editing by hand** — and see the result in the diagnostics, **without ever
changing the staged "which model wins" verdict.**

## Non-goals (the bright line)

- **No joint model selection.** The Stage 2 comparison table's `p`-value and
  highlighted-winner columns are computed **once** by the staged "Compare all"
  run and stay **frozen**. A joint refine never recomputes them. (This is Sam's
  constraint; it is the whole reason for the separation below.)
- **No write-back of joint curves to Stage 1.** Stage 1 remains the staged seed
  / single source of truth, so the staged verdict is always reproducible. A
  refine's moved curves live only in that model's refined fit.
- No new optimizer; reuse `refine_joint()` (`R/fit.R`).

## Design

### Two separate stores (the crux)

- `fits_store` (unchanged) — the **staged** per-model fits. Drives the Stage 2
  comparison table and, with `last_compare`, the frozen verdict. Written only by
  "Compare all" (and the manual staged Autofit/Simulate, as today).
- `refined_fits` (new) — a `reactiveVal(list())` of **joint-refined** fits, keyed
  by model name. Written only by the joint-optimize action. **Never read by the
  comparison table.** This is what keeps the verdict pure.

The model currently in focus is `input$model` (the existing picker). The
diagnostics and CIs show the **refined** fit for that model when one exists, else
the staged fit:

```r
display_fit <- reactive({
  m <- input$model
  refined_fits()[[m]] %||% fits_store()[[m]]   # refined wins when present
})
```

### UI — a "Tune / refine selected model" panel under the Stage 2 table

The user's instinct is that the table (with all params) is where tuning belongs,
so the controls sit **directly under the comparison table** and act on the
selected model. This folds the old Stage-3 "Explore by hand" accordion up next to
the table (it no longer lives separately in Stage 3).

```
Stage 2 · Compare interaction models
┌ Model comparison (STAGED verdict — decides the winner; read-only) ───┐
│ Model      max  sl1 sl2 ec501 ec502  a     b    SSR/Dev  df  p    ★  │
│ reference  786  21  1.5 .094  .33    —     —    1240     0   —       │
│ SA         786  21  1.5 .094  .33    1.55  —    1189     1  .004     │
│ DR ◀sel    786  21  1.5 .094  .33    1.62  .21  1180     2  .15      │
│ DL         786  21  1.5 .094  .33    1.49 -.08  1185     2  .31   ★  │
└──────────────────────────────────────────────────────────────────────┘
 Tune / refine selected model: DR
   a [1.62]  b [0.21]    [ Autofit a/b (staged) ]  [ Simulate ]
                         [ Optimize all params (joint) ]
   Joint refine: SSR 1180 → 1156  ✓   (curves now free; Excel-style)
   Refined (joint) params:  max 791  sl1 19.8 …  a 1.71  b 0.18
```

- **Autofit a/b (staged)** and **Simulate** — as today (curves fixed); write the
  staged `fits_store` entry for the model (so they *can* update the table row —
  they are staged, consistent with the verdict).
- **Optimize all params (joint)** — calls `refine_joint(fits_store()[[m]],
  engine_df(), n_starts, time_limit)`, seeded from the model's staged fit; stores
  the result in `refined_fits()[[m]]`; shows a **before → after** objective
  readout and the refined parameter values. **Does not touch `fits_store` or
  `last_compare`** → the table verdict is untouched.
- A small badge on the Stage 3 diagnostics (“showing joint-refined fit”) when
  `refined_fits()[[m]]` is in use, so it's never ambiguous which fit is on screen.

### Stage 3 — diagnostics for the (possibly refined) selected model

Stage 3 stays the diagnostics area (surface / isobole / obs-pred / CIs) but now
renders `display_fit()` — i.e. the joint-refined fit when present, else staged.
CIs: a joint-refined fit has all params free, so `param_ci()` yields real
intervals for every parameter (an improvement over the staged reference fit's
df=0 NaN CIs, for that model). The manual a/b grid moves up into the tune panel.

### Invalidation

- **Reference / response change:** clear `refined_fits` (and `fits_store`, as
  today) — the whole basis changed.
- **A Stage-1 curve edit:** clear `refined_fits` (the staged seed moved, so any
  joint refine is stale and must be re-run). The existing re-eval keeps the
  staged `fits_store` live for the table; refined fits are dropped, not
  re-evaluated.
- **Re-running "Compare all":** clear `refined_fits` (fresh staged baseline).

### What stays staged vs joint (summary)

| Thing | Method | Who writes it |
|---|---|---|
| Comparison table `p` / winner ★ | **staged**, frozen | `compare_all` only |
| Comparison table params/SSR rows | staged | `compare_all` + manual Autofit/Simulate |
| Selected model's diagnostics/CIs | refined if present, else staged | `optimize_all` (refined) |
| Refined params + before→after readout | **joint** | `optimize_all` only |

## Engine

No new engine function. `refine_joint(fit, df, lower=NULL, upper=NULL, n_starts,
time_limit)` already re-fits every parameter seeded from `fit`, returns a fit with
`joint = TRUE`. We call it with no bounds (all free). `param_ci()` already handles
an all-free joint fit.

## Testing / verification

1. **testServer happy path:** after `compare_all` + drain, `last_compare()` and
   the staged `fits_store` are set; select a model, fire `optimize_all`;
   `refined_fits()[[m]]$joint` is TRUE, its objective ≤ the staged objective, and
   **`last_compare()` and `fits_store()` are byte-for-byte unchanged** (verdict
   frozen). `display_fit()` returns the refined fit.
2. **Verdict isolation:** the `chosen` model and the comparison `p`-values are
   identical before and after a joint refine of any model.
3. **Invalidation:** a Stage-1 curve edit clears `refined_fits`; re-running
   `compare_all` clears `refined_fits`.
4. **Manual path:** Autofit/Simulate still update the staged `fits_store` entry
   and the table row.
5. **CIs:** a joint-refined model yields finite CIs for all params.

## Risks / notes

- Joint refine is slower and needs multistart on hard data — the existing
  `n_starts` / Thorough / `time_limit` controls apply. Before→after readout makes
  a bad local minimum visible (seeded from staged, it can only match or improve).
- Refining the **reference** model jointly = fitting curves to the mixture data
  with no interaction (all 5 curve params free) — valid, just the no-interaction
  joint curve fit.
- The staged table row and the refined params for the same model will differ
  (that's expected and the point); the badge + separate "Refined (joint) params"
  readout keep it unambiguous.

## Related memory

- `mixdra-staged-fitting` — staged selection is Sam's required method (the
  constraint this design protects).
- `excel-joint-fitting` — the joint/Excel method + that `refine_joint()` already
  exists; this re-surfaces it as a post-selection polish.
