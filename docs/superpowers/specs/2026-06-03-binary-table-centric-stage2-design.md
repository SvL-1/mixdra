# Binary tab: table-centric Stage 2 (auto-fill + row-select + one Optimize)

**Date:** 2026-06-03
**Status:** Design proposed; pending user review → implementation plan
**Supersedes the UI of:** `2026-06-02-binary-model-selection-loop-ux-design.md` and
`2026-06-03-binary-joint-refine-polish-design.md`. **Reuses their server core**
(`fits_store`, `refined_fits`, `display_fit`, the `optimize_all` observer, the
live stepper, `compare_fits`/`selection_chain_order`). The staged-selection
constraint from `mixdra-staged-fitting` (Sam) still governs.

## Context

Today's Stage 2 has a "Fit & compare all models" button + a results table, and a
separate "Tune / refine selected model" panel below it (model picker + a/b grid +
Autofit/Simulate + Optimize-all). The user finds the multiple surfaces
complicated and wants the **table itself** to be the single hub: it should fill
in automatically, you click a row to see that model, and one button optimises the
selected row.

## Goal

Collapse the Stage 2 table + picker + tune panel into **one live table with a
single action bar**:

1. The comparison table **auto-populates** as soon as both single curves are fit
   (and re-runs when a curve changes) — no "Compare all" button.
2. **Click a row** to make that model the one shown in the Stage 3 diagnostics
   (replaces the model picker).
3. **One "Optimize all params (joint)" button** acts on the selected row; its
   joint result appears as an extra column + drives the diagnostics.

## The staged-selection constraint (unchanged, still the bright line)

The winner (★) and the `p`-values are computed by `compare_fits` from the
**staged** `fits_store` only. The joint column is **display-only** and is NEVER
fed into `compare_fits`/`last_compare`. Model *selection* stays staged (Sam);
joint stays a per-model reporting polish.

## Design

### Stage 2 — one live table + an action bar

```
Stage 2 · Interaction models          (fills in automatically; α [0.05])
┌ Model      a     b     SSR    SSR(joint)  df   p      ★ ┐
│ reference  —     —     1240   —           0    —        │
│ SA         1.55  —     1189   —           1   .004      │
│▶DR  (sel)  1.62  0.21  1180   1156 ✓      2   .15       │
│ DL         1.49 −0.08  1185   —           2   .31    ★  │
└──────────────────────────────────────────────────────────┘
  [ Optimize all params (joint) ]   acts on the selected row
  Joint refine SSR (vs staged): 1180 → 1156  ✓
  Shared curves (from Stage 1): max 786 · sl1 21.1 · sl2 1.54 · ec50₁ .094 · ec50₂ .334
```

- **Auto-fill.** When both single curves exist (`frozen()`), the staged loop runs
  automatically — reference first (free: curves evaluated, no interaction), then
  SA/DR/DL — filling rows live (the existing one-model-per-tick stepper, just
  triggered automatically instead of by a button). It re-runs when `curve_params()`
  changes (a chemical re-fit). Curve changes are discrete button events from the
  Stage-1 panels (not slider drags), so re-running is cheap and not janky.
- **One row per model.** Columns: `Model · a · b · SSR/Dev · SSR(joint) · df · p`,
  winner row highlighted. The five **curve** params are identical across the
  staged rows (shared, fixed from the singles), so they're shown ONCE as a caption
  under the table rather than repeated per row. `SSR(joint)` is blank until that
  row is optimised; a `✓` marks an optimised row.
- **Row-click selection** (DT `selection = "single"`). The selected row is the
  "displayed model": it drives the Stage 3 diagnostics and is what the Optimize
  button acts on. After the loop completes, the **winner row is auto-selected**.
- **Action bar** under the table: a single **"Optimize all params (joint)"**
  button (acts on the selected row → `refine_joint` seeded from that model's
  staged fit → fills its `SSR(joint)` cell + drives the diagnostics), plus the
  before→after readout. `alpha` lives by the table (it sets the ★).

### Stage 3 — diagnostics (essentially unchanged)

Surface / isobole / obs-pred / CIs for the **selected** model, rendering
`display_fit()` (joint when that row is optimised, else staged) with the existing
"joint-refined" badge.

### Deliberate removals (please confirm in review)

To actually simplify, this design **removes**:
- the **"Fit & compare all models" button** (auto-fill replaces it);
- the **model picker** `selectInput` (row-click replaces it);
- the **manual a/b grid + "Autofit (a, b)" + "Simulate"** controls. Autofit is now
  redundant (every model is auto-fit); the hand-enter-a/b "Simulate" explore path
  is dropped. The user's notion of "simulate a row" is satisfied by **clicking the
  row** (it shows that model in the diagnostics). *If you want to keep hand-entry
  of hypothetical a/b values, say so and we'll retain a minimal version.*

### What the server reuses vs changes

**Reused unchanged:** `fits_store`, `refined_fits`, `refine_pre/post`,
`display_fit()` logic, the `optimize_all` → `refine_joint` flow, the
`refine_readout`/`refined_badge` outputs, the stepper, `compare_fits`,
`selection_chain_order`, all invalidation (reference/response, curve-change,
manual-supersede no longer applies since manual is gone).

**Changes:**
- **Auto-trigger** the stepper: replace `observeEvent(input$compare_all,…)` with
  an `observeEvent(curve_params(), …)` (guarded on `frozen()`) that clears the
  stores and arms the loop — so the table fills on freeze and refills on a curve
  change.
- **Selected model from the table:** a `selected_model <- reactive(...)` derived
  from `input$results_rows_selected` mapped to the row order
  (`intersect(c("reference","SA","DR","DL"), names(fits_store()))`); default to
  `last_compare()$chosen` once the loop finishes (set via DT proxy selection).
  Everywhere that read `input$model` now reads `selected_model()` —
  `current_fit`, `display_fit`, `optimize_all`, `interaction_help`.
- **Joint column:** `output$results` gains an `SSR(joint)`/`Dev(joint)` column
  populated from `refined_fits()` (display-only; the `p`/★ still come from staged
  `last_compare`). A `✓`/flag marks optimised rows.
- **Remove** the `observeEvent(input$model,…)`, `observeEvent(input$autofit,…)`,
  `observeEvent(input$simulate,…)`, the `read_ab()` helper, and the `val_a`/`val_b`
  `updateNumericInput` calls.

## Testing / verification

1. **Auto-fill:** in testServer, after both curves are simulated (`frozen()`),
   advancing the timer fills `fits_store` with all four staged fits and sets
   `last_compare()` — with NO `compare_all` input fired.
2. **Re-fill on curve change:** changing a chemical curve re-runs the loop
   (fits_store repopulated, last_compare refreshed at the new curves).
3. **Row-select drives display:** setting `input$results_rows_selected <- k`
   makes `selected_model()` the k-th model and `current_fit()`/`display_fit()`
   follow it.
4. **Verdict isolation (unchanged guarantee):** optimising the selected row sets
   `refined_fits()[[m]]`, fills the joint column, and leaves
   `last_compare()$comparison`/`$chosen` and the staged `fits_store` entry
   unchanged.
5. **UI:** the table output has `selection="single"`; the Optimize button +
   readout + badge are present; the compare_all button, model picker, a/b grid,
   Autofit and Simulate are GONE.

## Risks / notes

- **DT row-selection ↔ model mapping** must track the fixed row order; default-
  selecting the winner uses `DT::dataTableProxy` + `selectRows`. This is the main
  new plumbing.
- **Auto-fill compute:** 4 staged fits per curve change. Cheap (curves fixed,
  1–2 free params), and curve changes are discrete — acceptable. If a future
  dataset makes this slow, a debounce or an opt-out "pause auto-fit" toggle can be
  added (not now — YAGNI).
- **Reference CIs:** selecting the un-optimised reference row still shows NaN CIs
  (df=0) — the known pre-existing `param_ci` limitation, unchanged.

## Related memory

- `mixdra-staged-fitting` — staged selection is required (the constraint this
  preserves); update once shipped.
- `excel-joint-fitting` — joint stays a per-model reporting polish.
