# (De)potentiation tab — design

Date: 2026-11-19. Requested by Sam van Loon (via Jelle), agreed in conversation.

## Intent

The app gets a top-level split between two kinds of experiment:

- **(De)potentiation** — a full dose–response curve of stressor 1, measured
  without stressor 2 and with stressor 2 at one, two or three levels. Stressor 2
  can be anything: a chemical, microplastics, heat stress. The only question
  the data can answer is *does stressor 2 change the toxicity of stressor 1, and
  if so, does it become more or less toxic* (and, with several levels, from
  which level on).
- **Mixture toxicity** — the existing full mixture analysis (CA/IA, S/A, DR/DL,
  ternary), unchanged.

Why: a user who has stressor 2 in only one or a few fixed levels has a
potentiation experiment, not a mixture dose–response experiment, and must not be
pushed through the mixture models (Sam's issue #13). This tab gives that data a
correct, deliberately limited analysis.

The tab is generic: it never assumes what stressor 1 or 2 is. Names come only
from what the user types or uploads, with "Stressor 1" / "Stressor 2" as
defaults.

The method is Sam's, as implemented in his draft workbook
`DRCs different ratios FBSA with IMI and CPF.xlsx`, sheet
`75-25 drc CPF on FBSA IMI`, three comparison blocks in columns L–T, V–AD and
AF–AN (labelled "25%", "75%" and "50%" there). The workbook is a draft: some
labels in it, including compound names, are unreliable, so nothing in the
package or app takes names from it. Results must equal the workbook, except
where the workbook has a formula error (see "Reference values").

## Goals

1. Navbar becomes **Introduction | Single Stressor | (De)potentiation |
   Mixture toxicity**. "Mixture toxicity" is the current "Multiple stressors"
   tab, renamed (user-facing labels only, as the #12 rename was).
2. A new **(De)potentiation** tab: upload, fit, plot, test, conclude (below).
3. A new exported engine function `analyse_potentiation()` that reproduces the
   workbook's numbers.
4. The workbook comparison is a validation oracle, visible in the validation
   report like the existing ones.

## Non-goals

- No mixture models, no curve for stressor 2, no toxic units.
- No test of the slope. Sam tests the EC50 only; slopes stay free per curve.
- No multiple-testing correction across levels (the workbook has none).
- No Lower/Upper/Fix grid in this tab in the first version (so the default
  bounds below are fixed in this tab for now; a grid can make them adjustable
  later, as on the other tabs).
- No automatic routing of uploads between tabs (#13's quality check is a later
  step; this tab is what it would route to).
- No change to the #15 workflow revision of the mixture tab.

## Anti-goals

- Numbers that differ from Sam's (corrected) workbook for the same data.
- Conclusions beyond "influence / no influence, more / less toxic, from level x".
- Stressor-specific names or assumptions anywhere in the tab, example or docs.

## Data format

One CSV, one row per observation:

| column | meaning |
|---|---|
| `Series` | free-text label of the curve, e.g. `Stressor 1`, `+ Stressor 2 level 1`, `Heat 30 °C` |
| `C1` | concentration of stressor 1 (≥ 0) |
| `Res` | continuous response — **or** `Exposed` + `Affected` for quantal |

How stressor 2 was applied (fixed concentration, fixed ratio, temperature) is
deliberately not encoded: it only labels the series. The **reference series**
(stressor 1 alone) is chosen in the UI and defaults to the first series in the
file. The other series are taken in file order, which the user is told should be
increasing stressor-2 level. At least two series; each needs enough rows to fit
three parameters (same minimum as `fit_single()`, > 3 rows).

Validation follows the existing `validate_upload()` style: required columns per
response type, numeric `C1 ≥ 0`, non-empty `Series`, quantal `0 ≤ Affected ≤
Exposed`.

### Controls — the workbook approach

In the workbook the 20 control observations (`C1 = 0`) are listed in both
curves of every comparison, so they are counted twice (n = 90 for 45 + 45
rows). We do the same, and say so:

- A series that has its own `C1 = 0` rows uses them (e.g. "with heat" has a
  heat-only control).
- A series without any `C1 = 0` rows **borrows the reference series' controls**,
  which are then counted in both curves of that comparison.
- The app states this plainly next to the results, per comparison, with n, e.g.
  *"Controls: the 20 control observations of Stressor 1 are used for both
  curves, as in the Excel workbook, so they are counted twice (n = 90). This is
  the approach used here; it makes the test slightly more liberal than counting
  them once."*

## Method

Model: the existing three-parameter log-logistic, `Y = max / (1 + (C/EC50)^b)`
(`ll3_predict()`), with `fit_single()`'s default bounds (max > 0, slope in
[1e-3, 50], EC50 > 0; quantal max ≤ 1). The slope bound of 50 is the app-wide
default, not specific to this tab.

For each non-reference series *k*, compared with the reference series *r*:

- **Run 1 (separate):** *r* and *k* each get their own max, slope, EC50.
  The two curves share nothing, so this is two independent fits; the reference
  curve is fitted once and reused for every comparison.
  `obj1 = obj(r) + obj(k)`.
- **Run 2 (shared EC50):** one EC50 for both curves; max and slope stay per
  curve. Five parameters, one joint fit. `obj2`.
- **Test:** likelihood ratio, df = 1.
  - continuous: objective = residual sum of squares,
    `X² = n · ln(obj2 / obj1)`, n = rows of *r* + rows of *k* (including any
    borrowed controls). This is the workbook's formula and `lr_test()`'s.
  - quantal: objective = binomial deviance (`obj_deviance()`),
    `X² = obj2 − obj1`, matching the mixture tab's binary convention.
  - `p = pchisq(X², 1, lower.tail = FALSE)`, α = 0.05.
  - X² is floored at 0: run 2 is nested in run 1, so a negative value is
    optimiser noise (the workbook's third block shows X² = −1.3e−5).
- **Direction:** `ratio = EC50_k / EC50_r` from run 1. Significant and
  ratio < 1 → *more toxic*; significant and ratio > 1 → *less toxic*;
  otherwise *no significant influence*.
- **Onset:** series are taken as increasing stressor-2 level in file order. An
  onset is named only when, from some series on, every later series is
  significant in the same direction and no earlier one is ("…more toxic from
  *Heat 30 °C* on"). Otherwise no onset is named and the summary lists the
  significant series per direction and the non-significant ones, and says the
  effect does not increase consistently with the level. With one non-reference
  series the summary is simply that comparison's verdict.
- **Extrapolated EC50:** when a series' run-1 EC50 lies above its highest
  tested `C1` (or below its lowest positive `C1`), its row and the summary carry
  a short caveat that the EC50 — and so the shift — is extrapolated. The verdict
  itself is not changed. (The existing `ec50_range_note()` talks about toxic
  units and Fix boxes, so this tab gets its own wording.)

### Optimisation

Use `optim` L-BFGS-B with `parscale`, as `fit_single()` does. A single seed
finds the optimum in all three reference blocks (the builder check confirmed
the level-3 optimum from 225 starts). As cheap insurance, start run 2 from
shared EC50 = run-1 EC50 of *r*, of *k*, and their geometric mean (max/slopes
from run 1), keep the lowest objective, and seed run 1's *k* fit from the run-2
solution too, so that obj1 ≤ obj2 holds.

Several optima sit on the slope bound of 50 (a very steep reference curve), but
lifting it to 1000 changes SS by < 0.01 and X² not at all (checked on all three
blocks).

## Reference values

The example and the oracle use generic labels. Mapping to the workbook:

| here | workbook block | columns (data rows 30–74 reference, 81–125 series) |
|---|---|---|
| Level 1 | "25%" | L–T (conc M, response N) |
| Level 2 | "50%" | AF–AN (conc AG, response AH) |
| Level 3 | "75%" | V–AD (conc W, response X) |

Re-fitted in R from the workbook data (scratch scripts, not in the repo):

| comparison | workbook SS₁ | workbook SS₂ | workbook X² | R |
|---|---|---|---|---|
| Level 1 | 301 550.972 (M25) | 341 687.917 (P25) | 11.2463 (M27), p = 0.0008 | identical |
| Level 2 | 393 465.224 (AG25) | 393 465.167 (AJ25) | −1.3e−5 (AG27), ≈ 0 | identical (X² ≈ 0) |
| Level 3 | 379 874.993 (W25) | 399 798.700 (Z25) | 4.6007 (W27) | SS₁ identical; SS₂ 411 358.37, X² 7.166, p = 0.0074 |

**The level-3 workbook SS₂ contains a formula error.** In run 2 of that block
only AC30 uses the block's own reference parameters (`$Z$17:$Z$19`); AC31:AC74
use `$P$17:$P$19`, the level-1 block's run-2 parameters. So SS₂ = 399 799 is
not a shared-EC50 fit of this block. Correcting AC31:AC74 to `$Z$17:$Z$19` and
re-solving gives SS₂ = 411 358.37, X² = 7.166, p = 0.0074 — same conclusion
(significant, less toxic), different numbers. The level-1 and level-2 blocks'
formulas are correct.

In all three blocks the reference rows are identical, and the 20 controls of
each series are the same values as the reference's. Reference EC50 ≈ 0.102;
level 1: 0.134 (less toxic, significant); level 2: 0.101 (no influence);
level 3: 0.361 (less toxic, significant, extrapolated: the highest stressor-1
dose in that series is 0.1).

## Engine API

```r
analyse_potentiation(df, reference = NULL,
                     response = c("continuous", "quantal"), alpha = 0.05)
```

`df`: columns `Series`, `C1`, and `Res` or `Exposed`/`Affected`. `reference`:
a series label, default the first in `df`. Returns a list with

- `reference`, `response`, `alpha`;
- `curves`: the run-1 fit per series (par max/slope/ec50, objective, n,
  convergence, `controls_borrowed`, `ec50_extrapolated`);
- `comparisons`: one row per non-reference series — `series`, `n`, `obj1`,
  `obj2`, `chi`, `df`, `p`, `ec50_reference`, `ec50_series`, `ec50_shared`,
  `ratio`, `significant`, `direction`, `controls_borrowed`,
  `ec50_extrapolated`;
- `shared`: the run-2 fit per comparison (for plotting on demand);
- `summary`: the conclusion text (one or two sentences).

Internal helpers are the builder's call.

## The tab

```
(De)potentiation
 sidebar:  Stressor 1 name + unit, Stressor 2 name   (default "Stressor 1/2";
                                                      label axis and summary)
           Response type (Continuous / Quantal)
           Download template · Upload CSV · upload errors
           Reference series (select, default first)
 main:     Plot — stressor-1 response vs concentration, one colour per series,
                  observed points + run-1 fitted curve; log x-axis with the
                  controls handled as plot_dose_response() does
           Table — per series: max, slope, EC50, EC50 ratio vs reference,
                   extrapolation caveat
           Tests — per comparison: n, SS/deviance run 1 and 2, X², p, verdict
           Controls note (see above)
           Conclusion — for the example data, roughly: "Stressor 1 is
                        significantly less toxic with Level 1 (EC50 0.102 →
                        0.134, p = 0.0008) and Level 3 (EC50 0.102 → 0.361,
                        p = 0.0074, extrapolated); no significant influence
                        with Level 2. The effect does not increase
                        consistently with the level."
```

An example dataset (series `Stressor 1`, `Level 1`, `Level 2`, `Level 3`, in
that order) is loaded until the user uploads a file, as the mixture tab does.
In the example the controls are listed once, under `Stressor 1`, so the
borrowing note is visible and the numbers still equal the workbook.

Introduction tab: a short paragraph on which tab to use — one stressor at a few
fixed levels → (De)potentiation; dose–response of every stressor and of their
mixtures → Mixture toxicity.

## Validation

- Fixture `tests/testthat/fixtures/potentiation/potentiation_example_continuous.csv`
  (generic labels, controls listed once). The workbook itself is not committed:
  it is a draft with known errors, and the same observations are already public
  in `tests/testthat/fixtures/binary/mps_cpf_fbsa/`. Provenance (file, sheet,
  cells, and the formula correction) goes in the oracle rows.
- `inst/validation/oracles.csv` rows for SS₁, SS₂ and X² of the three
  comparisons, sourced to the workbook cells above; relative tolerance, tight
  (level 1 agrees to 1e−8). Level 2's X² is a "≈ 0" check.
- Level 3: SS₁ is a normal workbook match. SS₂ and X² are asserted against the
  **corrected** workbook values (411 358.37 / 7.166), with the source naming
  the correction (AC31:AC74 `$P$17:$P$19` → `$Z$17:$Z$19`). Because the
  corrected value was re-solved here rather than by Sam's Solver, the report
  must not count it as an independent match: record it as a regression pin
  (the report's "=" kind) with a note, unless Sam re-solves the corrected
  workbook, in which case it becomes a match.
- Regenerate `VALIDATION.md` (CI checks it).
- Unit tests: recovery from simulated data (no shift → not significant; EC50
  shifted → significant, right direction), borrowed vs own controls, quantal
  path, X² floor, reference choice, onset/summary cases (monotone, not
  monotone, mixed directions, single series), extrapolation caveat, validation
  errors. Run filtered tests only.

## Constraints

- Base R + existing dependencies only; app modules follow the existing
  `*_ui()/*_server()` pattern in `R/app-*.R`.
- The rename touches user-facing strings only, as #12 did; internal
  `campaign_*` names stay.

## Left to the builder

Internal helper structure, plot styling, table layout, how the example dataset
is stored, exact wording of help texts within the meaning above.
