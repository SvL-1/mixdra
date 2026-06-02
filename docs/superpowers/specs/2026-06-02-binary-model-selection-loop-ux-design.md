# Binary tab: model selection as a joint-fit loop (UX spine)

**Date:** 2026-06-02
**Status:** Design approved; pending spec review → implementation plan

## Context

The automatic choice of interaction model is the centerpiece of this app, but
the current Binary tab UX buries it and contains two competing fitting
philosophies bolted together:

- **Stage 2 "Find best model"** (`analyse_mixture`) is **staged**: the curve
  parameters (`max`, `slope*`, `ec50*`) are estimated once from the
  single-compound rows and **frozen**, then only the interaction terms (`a`,
  `b`) are fitted per model, and reference/SA/DR/DL are compared by
  likelihood-ratio test. "Find best" is a small button tucked into the Stage 4
  results card.
- **Stage 3 "Optimize all parameters"** is the **joint** Excel-faithful fit, but
  it runs on **one** displayed model at a time as a manual post-selection polish
  (see `2026-06-02-binary-joint-fit-design.md`).

So "best" is currently decided on *frozen-curve* fits whose SSRs don't match the
Excel, while the Excel-style joint fit is a separate manual afterthought. That
split is the root of the "the flow is vague" complaint.

The original Excel workbook fit **every model fully jointly** (one SSR per
model, curves + interaction together) and compared models on those joint fits.
That is also the statistically standard nested comparison ("does the full model
beat the reduced model"). Our staged approach was a deliberate robustness choice
(joint curve+interaction fits on mixture data are numerically harder), **not**
something the comparison test requires.

**This design supersedes the "No change to model selection" non-goal of
`2026-06-02-binary-joint-fit-design.md`.** We move model selection itself to the
joint regime and make that loop the spine of the UI.

## Goal

Restructure the Binary tab around the model-selection loop the user described:

> seed the curves → for each interaction model, fully fit it and record its
> results → compare all models in a table → highlight the winner.

Concretely:

1. The selection loop fits **all four models jointly** (every parameter free),
   **seeded** from the single-compound curves (robust start, not a freeze).
2. Each model's full results are recorded (all params, joint SSR/deviance, df,
   n, R²) and shown **live, row by row**, as each fit completes.
3. The existing nested test selects and highlights the parsimonious winner.
4. The standalone "Optimize all" stage is **retired** — optimization now happens
   inside the loop for every model, dissolving the staged-vs-joint split.

## Non-goals

- **No new optimizer engine.** Reuse `fit_model()` (L-BFGS-B + Nelder-Mead
  fallback, `n_starts`/Thorough multistart, `time_limit`). "Joint" here means
  fitting all params per model, not a new solver.
- **No change to the per-point model math** (CA/IA × reference/SA/DR/DL).
- **No change to the ternary subproject.**
- **Not removing manual exploration.** The by-hand single-model controls
  (Autofit `a`/`b`, Simulate) survive, demoted to an "explore" affordance.

## The fitting loop (engine)

A joint variant of `analyse_mixture` (new mode or new function, e.g.
`analyse_mixture_joint()` or a `joint = TRUE` arg):

1. **Seed** the curves from the single-compound rows exactly as today
   (`fit_curve_from_singles`), or take a user-supplied `start`. These become
   *starting values*, no longer fixed.
2. **Warm-started nested chain** (important for soundness and robustness):
   - Fit **reference** jointly (curves free), seeded from the curve seed.
   - Seed **SA** from the reference joint fit + `a = 0`; fit jointly.
   - Seed **DR** and **DL** from the SA joint fit + `b = 0`; fit jointly.
   Warm-starting each child from its parent's joint optimum guarantees the child
   can do **no worse** than the parent (it can always reproduce the parent by
   driving the extra term to 0), so the child objective ≤ parent objective and
   the LR statistic stays ≥ 0. This both makes the nested comparison sound and
   sidesteps the "a more complex model accidentally fit worse" local-minimum
   trap.
3. **Compare** with the existing `lr_test` / `select_parsimonious`. These are
   **regime-agnostic** — they need only each fit's `objective` and `df`. Under
   the joint regime the df are the full counts (binary: reference = 5, SA = 6,
   DR/DL = 7), but the *differences* are still 1, so the test form is unchanged
   and becomes the standard unconditional nested comparison.

No change to `compare.R` is required for the test to work.

### Open decision — test statistic

`lr_test` currently uses the asymptotic LR statistic `n·log(SSRp/SSRc)` for
continuous data; the Excel reports an **F-test** (`F`, `P`). The selection
*logic* is identical in form either way. **Recommendation:** keep the existing LR
statistic for now (consistent with the binary-quantal deviance path and the
ternary subproject); offer an F-statistic swap as a small, separate refinement
if exact Excel-number parity is wanted. Flagged here so it's an explicit choice,
not an accident.

## User-facing design (3 stages)

```
Stage 1 · Single curves (seed)
  Two chemical panels (Autofit / Simulate), unchanged fitting.
  Reworded: these seed the joint fits — they are a starting point,
  not frozen truth. Loop workspace appears once both are seeded.

Stage 2 · Compare interaction models      ← THE HERO
  [ Fit & compare all models ]   alpha [0.05]
  ┌── Model comparison ───────────────────────────────────────────┐
  │ Model      a       b      SSR/Dev   df   F/χ²    p      ◄ best │
  │ reference  —       —      1240.3    5     —      —             │
  │ SA         1.55    —      1188.7    6    8.9   0.003          │   ← rows
  │ DR         1.62   0.21    1180.1    7    2.1   0.15           │   ← fill
  │ DL         1.49  -0.08    1185.0    7    1.0   0.31    ◄ best  │   ← live
  └────────────────────────────────────────────────────────────────┘
  Each row appears as that model's joint fit finishes; the winner
  is highlighted once the chain completes.

Stage 3 · Inspect a model
  Pick any row (defaults to the winner) → surface / isobole /
  obs-vs-pred / CIs / plain-language interaction help for that model.
  "Explore by hand" (collapsed): manual Autofit a/b, Simulate.
```

**Behavior:**

- One primary action ("Fit & compare all models") runs the whole loop. `alpha`
  lives next to it (moved out of the old results card). The `n_starts` /
  Thorough / `time_limit` knobs stay in the sidebar accordion and apply to the
  loop.
- The comparison **table is promoted to the centerpiece** of Stage 2 (params ·
  joint objective · df · statistic · p), best row highlighted — it directly
  answers "which interaction model wins."
- Rows fill **live** as each model's joint fit completes (see below).
- Stage 3 is the drill-down: selecting a model shows its diagnostics, CIs, and
  help. The manual single-model controls live here, collapsed, as optional
  exploration rather than the main path.

### What is removed

- The standalone **Stage 3 "Optimize all parameters"** card and its
  `refine_joint` button-flow — optimization is now intrinsic to the loop.
- The staged **"Find best model"** button buried in the old Stage 4 results
  card — its job is now the Stage 2 hero action, run jointly.

## Live row-by-row updates (Shiny mechanism)

A single blocking `analyse_mixture` call cannot paint the table between models,
because Shiny flushes outputs only after the observer returns. Two options:

- **Stepped state machine (recommended, no new deps):** hold a queue of
  remaining models in a `reactiveVal`. An observer fits the *next* model, writes
  its fit into `fits_store`, then re-schedules itself with
  `shiny::invalidateLater(0)` until the queue drains. Each cycle flushes the UI,
  so the table grows one row per model. `select_parsimonious` + the
  best-row highlight run when the queue empties (it needs all four). Progress via
  `withProgress` per step.
- **Async (`promises`/`future`):** offload fits to a worker and update on
  resolution. More moving parts and a new dependency; deferred unless the stepped
  approach proves too janky.

The warm-started chain (reference → SA → DR/DL) maps naturally onto the queue:
each step seeds from the previous step's result.

## App wiring (`R/app-binary.R`)

- Replace the Stage 2/3/4 structure with the 3-stage layout. Stage 1 panels and
  their `inject*` write-back are largely unchanged.
- Stage 2 observer drives the stepped loop into `fits_store` and stores the
  comparison/chosen in `last_compare` when complete.
- Stage 3 reads `fits_store()[[selected]]` for diagnostics/CIs/help; the
  "explore by hand" Autofit/Simulate observers stay but write into the same
  store (so a hand fit can replace a row).
- Existing invalidation (clear `fits_store` on curve/reference/response change)
  is unchanged and still correct.
- `refine_joint` / the bounds-grid helpers from the previous joint-fit feature
  may be **reused internally** as the per-model joint fit primitive, or retired
  if the loop calls `fit_model` directly with curves free.

## Testing / verification

1. **Joint analyse unit:** the new loop returns four fits with full df (5/6/7),
   a comparison frame, and a chosen model; child objective ≤ parent objective
   (warm-start monotonicity invariant).
2. **Comparison validity:** `lr_test`/`select_parsimonious` produce the expected
   selection given joint SSRs with full df (differences still 1).
3. **Synthetic recovery** (reuse `simulate_*` + recovery oracle): noise-free
   binary data from known params; the joint loop recovers params-in ≈
   params-out and selects the true model.
4. **Live update integration** (testServer): driving the loop fills `fits_store`
   incrementally (one model per cycle) and lands on the chosen model.
5. **Excel cross-check (acceptance):** real CPF/IMI binary data; the winning
   model's joint params land near the published Excel values (approximate — GRG
   vs. our optimizer).

## Risks / notes

- **4× full joint fits are slower** than 4 cheap `a`/`b` fits; `n_starts` /
  Thorough matter more. Warm-starting the chain mitigates both speed and local
  minima. The live table makes the cost visible rather than a frozen spinner.
- **Shared `max`:** the joint fits re-estimate a single shared `max` from the
  full data, resolving the staged "average of two maxes" heuristic for the
  reported result.
- **Excel parity is approximate** by design (different optimizer; CA bisection
  non-smoothness).

## Related memory

- `excel-joint-fitting.md` — Excel fits all params jointly; engine already
  supports joint via unfixing. **Update:** model *selection* now uses joint too.
- `mixdra-staged-fitting.md` — staged fitting; **the loop is no longer frozen**,
  staged curves are now only a seed. Needs updating once shipped.
- `mixdra-synthetic-recovery.md` — recovery-oracle test pattern reused.
- `binary-workbook-structure.md` — fixtures for the Excel cross-check.
