# Binary tab: final "Optimize all parameters" (joint) stage

**Date:** 2026-06-02
**Status:** Design approved; pending spec review → implementation plan

## Context

The Binary ("dual") Mixture tab fits mixture dose-response models in a **staged**
way: Stage 1 fits each chemical's single curve, the curves are frozen (with a
shared `max` = the average of the two single-fit maxes), then Stage 2 fits only
the interaction parameters (`a`, `b`) to the mixture rows, and compares
reference/SA/DR/DL via likelihood-ratio tests ("Find best model").

The original Excel workbook (`recieved/.../Model4Isoplane.xlsm`) instead fits
**all parameters jointly**: its fitting sheet ("CA - cont. data Adv C (3) free")
has Solver adjust Max + all slopes + all EC50s + all interaction terms
(`solver_adj = Q2:Q12`) against a single global SSR target (`solver_opt = Q20 =
SUM(R33:R416)`, minimize). It is a single-start GRG Nonlinear solve from
hand-entered starting values. Sam's earlier R attempt
(`MixTox_shiny_v2/modules/binary_module.R`) was also staged (and incomplete) and
explicitly noted being unable to recreate Excel's GRG solver.

So the per-point model math is faithful to the originals, but the app's
*parameter estimation* deliberately differs from the Excel. This feature closes
that gap **without removing the staged workflow the team values**.

## Goal

Add a final, optional stage to the Binary tab that **jointly optimizes every
parameter** of the currently displayed model — seeded from the staged fit — so
the user can reproduce the Excel's joint-fit objective. A per-parameter
lower/upper bounds grid lets the user constrain or pin any parameter (including
`a`/`b`), spanning the full spectrum from "freeze the curves, fit only `a`/`b`"
(today's staged behavior) to "free everything" (Excel-style).

## Non-goals

- **No change to model selection.** Stage 2's four-model comparison and
  "Find best model" stay staged (frozen curves), preserving the clean "only the
  interaction differs between models" interpretation. Optimize-all is a
  **post-selection polish on the single displayed model**.
- **No new optimizer engine.** We reuse the existing `fit_model()` machinery
  (L-BFGS-B + Nelder-Mead fallback, multistart via `n_starts`/Thorough,
  `time_limit` guard). We do **not** chase a bit-exact GRG clone; the CA model's
  bisection inner solver makes the objective only piecewise-smooth, which is why
  derivative-free/robust methods are the right choice. "Matching the Excel" means
  reaching the same minimum objective, not using the same-branded solver.
- **No change to the ternary subproject.**

## User-facing design

The Binary tab keeps its current layout. Stage 1 (single curves → freeze) and
Stage 2 (interaction model: Autofit / Simulate / Find best / comparison table)
are **unchanged**. A new card is added at the bottom of Stage 2, revealed once a
model is displayed:

```
┌─ Optimize all parameters (joint) ─────────────────────────────┐
│  Refines the displayed model by fitting every parameter at     │
│  once, seeded from the current fit. Fix any parameter by        │
│  setting its Lower = Upper.                                     │
│                                                                 │
│  Parameter   Meaning            Lower     Upper     Value(start)│
│  max         control response   [    ]    [    ]    786.40      │
│  slope1      …                  [    ]    [    ]    21.10       │
│  slope2      …                  [    ]    [    ]    1.54        │
│  ec501       …                  [    ]    [    ]    0.094       │
│  ec502       …                  [    ]    [    ]    0.334       │
│  a           interaction        [    ]    [    ]    1.546       │
│  b           (DR/DL only)       [    ]    [    ]    —           │
│                                                                 │
│  [ Optimize all params ]        SSR: 1240.3 → 1188.7  ✓ improved│
└────────────────────────────────────────────────────────────────┘
```

**Behavior:**

- The grid is pre-filled with the **displayed model's current parameter values**
  (the staged fit) as the starting point. Lower/Upper start blank = free.
- The `b` row is shown only for DR/DL (matching the existing param-grid logic).
- Clicking **Optimize all params** runs one joint optimization over all of that
  model's parameters, seeded from those values, honoring any bounds, using the
  existing **n_starts / Thorough / time_limit** controls.
- Setting a row's **Lower = Upper** pins that parameter at that value — so
  "freeze the curves, fit only `a`/`b`" is reproducible here, as is anything in
  between.
- A **before → after** objective readout (SSR for continuous, Deviance for
  quantal) shows the staged value vs. the joint value. Seeded from the staged
  fit, the joint result can only match or improve it — a built-in sanity check.
- The result becomes the displayed model's current fit, so **Stage 3
  diagnostics (surface, isobole, obs-vs-pred, CIs)** and the Stage 2 objective
  readout update to reflect it.

## Technical design

### Engine change (`R/fit.R`) — minimal

1. **Allow `a`/`b` in `lower`/`upper`.** Relax the validation that currently
   rejects bound names outside the base curve params, so the bound names may be
   any model parameter. Curve-param positivity defaults are unchanged; `a`/`b`
   still default to ±Inf when unbounded.
2. **Fixing stays via the `fixed` argument, resolved in the app layer.** L-BFGS-B
   errors when `lower == upper`, so a pinned parameter is NOT expressed as an
   equal bound to the optimizer. Instead the app translates each grid row:
   - Lower == Upper (both present and equal) → parameter goes into `fixed` at
     that value (and `start[param]` is set to it).
   - Otherwise → passes through as `lower`/`upper`.
   This reuses `fit_model`'s existing `fixed` machinery (which already supports
   fixing `a`/`b`, since they are in `spec$params`) and avoids the
   `lower == upper` error. No engine change is needed for fixing itself.

The joint fit is then:

```r
fit_model(df, reference, deviation, response,
          start    = <displayed model's params>,
          fixed    = <names where Lower == Upper>,
          lower    = <remaining named bounds>,
          upper    = <remaining named bounds>,
          n_starts = n_starts_eff(), time_limit = input$time_limit)
```

with nothing fixed by default.

### New helpers (pure, unit-testable)

- `collect_bounds_all(values, params = c("max","slope1","slope2","ec501","ec502","a","b"))`
  in `R/app-io.R`: like the existing `collect_bounds()` but over the full
  7-parameter set; returns `lower`/`upper` named numerics (blank/NA dropped).
- A split helper (in `R/app-io.R` or the engine) that, given the grid's
  lower/upper and the start vector, returns `list(fixed = <names>, start =
  <values incl. pinned>, lower = <rest>, upper = <rest>)`.
- Optionally a thin `refine_joint(fit, df, grid_inputs, n_starts, time_limit)`
  wrapper in the engine that performs the split and calls `fit_model`, returning
  the enriched fit. Keeps the Shiny observer thin.

### App wiring (`R/app-binary.R`)

- New `output$optimize_panel` UI card with the 7-row bounds/fix grid, the
  **Optimize all params** button, and the before→after readout. Reuses the
  existing `param_row()`/grid styling; the `b` row uses the same DR/DL
  conditional as the Stage 2 grid.
- New `observeEvent(input$optimize_all, …)`: captures the pre-refine objective,
  reads the grid via `collect_bounds_all()`, runs the split + `fit_model` inside
  `shiny::withProgress`, and on success **overwrites the displayed model's entry
  in `fits_store`** so all downstream outputs reflect it. Errors surface via
  `showNotification`, matching the Autofit/Simulate pattern.
- **Freeze invalidation is automatic:** the observer that clears `fits_store` on
  curve/reference/response change already covers the joint fit, since it lives in
  the same store. No new invalidation logic.

### Confidence intervals (`R/summary.R`)

`param_ci()` already computes the Hessian over **all** parameters in `fit$par`,
so a joint fit flows through it correctly (all params genuinely were free — which
is exactly the Hessian's assumption). One refinement: for any parameter the user
**pinned** (Lower == Upper), display its fixed value with **no interval** (it was
not estimated), instead of a misleading SE. The existing `output$cis` path is
otherwise unchanged.

## Testing / verification

1. **Engine units** (`tests/testthat/test-fit*.R`): `fit_model` accepts bounds on
   `a`/`b` (naming `a` in `lower` no longer errors); bounding and out-of-range
   start clamping work for an interaction parameter.
2. **Fix-routing unit** (`tests/testthat/test-app-io.R`): `collect_bounds_all()`
   + split helper turn Lower==Upper rows into `fixed` (at that value) and the
   rest into `lower`/`upper`.
3. **Synthetic recovery** (aligns with existing `simulate_*`/recovery oracle):
   noise-free binary data from known `max/slopes/ec50s/a`; joint optimize from a
   perturbed start recovers params-in ≈ params-out within tolerance.
4. **Monotonic-improvement invariant:** on the same data, joint objective ≤
   staged objective (joint is seeded from the staged fit).
5. **Excel cross-check (acceptance):** feed the real binary CPF/IMI data
   (`Raw_data_binary.xlsx`) and confirm the joint fit lands near the published
   Excel params (Max ≈ 786.4, a ≈ 1.546, etc.). Approximate (GRG vs. our
   optimizer) but should be close — the "did we match the Excel" check.
6. **App wiring** (`tests/testthat/test-app-modules.R`, testServer): Optimize-all
   updates the displayed model's `fits_store` entry; a curve/reference change
   clears it via the existing freeze invalidation.

## Risks / notes

- **Joint fits are harder than 2-parameter interaction fits.** The existing
  `n_starts`/Thorough/`time_limit` guards become more important; the new card
  shares those Stage 2 controls. Default `n_starts = 1` may land in a local
  minimum on difficult data — the before→after readout and Thorough toggle are
  the user's levers.
- **Averaged `max`** in the staged seed is unchanged; the joint fit naturally
  re-estimates a single shared `max` when `max` is left free, which resolves the
  heuristic for the polished result.
- **Excel cross-check is approximate**, not bit-exact, by design (different
  optimizer; CA bisection non-smoothness). Acceptance is "close to published
  values," not equality.

## Related memory

- `excel-joint-fitting.md` — verified Excel joint-fit method + that the engine
  already supports joint via unfixing.
- `mixdra-staged-fitting.md` — the staged approach being preserved.
- `mixdra-synthetic-recovery.md` — the recovery-oracle test pattern reused here.
