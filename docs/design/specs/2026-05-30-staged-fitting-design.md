# Staged fitting: fix single-compound curves, then fit interaction params

**Date:** 2026-05-30
**Status:** Approved, ready for implementation plan

## Problem

When `analyse_mixture()` fits the binary or ternary models, it currently uses the
single-compound curve fits only as *starting seeds* (`seed_from_singles()`) and
then re-fits **all** parameters jointly to the full mixture dataset — the curve
parameters (`max`, `slope*`, `ec50*`) and the interaction parameters (`a`, `b`,
`b1/b2/b3`) move together.

This is not the intended workflow. The correct procedure (per Sam, and as
described in the MixTox method paper, `papers/1-s2.0-S0304389425034132-main.pdf`)
is **staged**:

1. Fit `max`, `slope*`, `ec50*` to the **single-compound data only**.
2. **Fix** those curve parameters.
3. Fit **only** the interaction parameters (`a`, `b`, …) to the mixture data.

Besides being the correct estimator, this fixes the large majority of parameters
before the expensive mixture fit, so it is also faster.

## Scope

- Applies to **both** binary (2-chem) and ternary (3-chem) analyses.
- Staged fitting becomes the **only** path in `analyse_mixture()` (no opt-out
  flag). `fit_model()` keeps its existing low-level `start`/`fixed` flexibility,
  so a joint fit is still reachable directly if ever needed.
- **Out of scope:** which deviation models to fit for ternary (SA-only vs also
  DR/DL). That stays with the on-hold ternary subproject. For ternary, DR/DL are
  not part of the intended workflow, but this change does not remove them — it
  only changes *how* whatever models are fit get fit.

## Key facts that make this clean

- On single-compound rows (every chemical but one at 0), every mixture model —
  reference, SA, DR, DL, for both CA and IA — collapses **exactly** to the
  independent 3-parameter log-logistic curve. The interaction params (`a`, `b`)
  never appear. So the curve params are fully identified by the single-compound
  rows, and stage 1 is reference-vs-CA/IA-agnostic.
- Model selection is unaffected: `lr_test()` and `select_parsimonious()` use the
  **difference** in free-parameter counts. The base params are shared across all
  four models and fixed identically, so they cancel. There is no AIC/BIC in the
  codebase, so nothing depends on the absolute `df`.

## Design

### Stage 1 — `fit_curve_from_singles()` (new internal helper)

```
fit_curve_from_singles(df, reference, response, lower, upper, n_starts, time_limit)
  -> named numeric vector (max, slope1..n, ec50..n)
```

- Select control + single-compound rows: `rowSums(conc_cols > 0) <= 1`.
- Fit the `reference` deviation on that subset via `fit_model()`, forwarding
  `lower/upper/n_starts/time_limit`. This yields **one shared `max`** plus
  per-compound `slope`/`ec50`.
- Use the caller's `reference` (CA/IA) — irrelevant on single rows, but keeps one
  code path.
- **Fallback:** if the single-compound subset is too small to fit (a chemical
  with < 4 distinct single-compound points) or the fit errors, fall back to the
  existing `seed_from_singles()` heuristic for that estimate. Never errors.

### Stage 2 — interaction fit (per deviation), in `analyse_mixture()`

For each deviation in `reference, SA, DR, DL`:

```
fit_model(df_full, reference, dev, response,
          start = c(stage1_base, interaction_defaults),
          fixed = base_param_names,            # max, slope*, ec50*
          lower = lower, upper = upper,
          n_starts = n_starts, time_limit = time_limit)
```

- Only interaction params (`a`, `b`, `b1/b2/b3`) are free.
- `reference` deviation has no interaction params → 0 free params → returns the
  fixed prediction with `df = 0`.
- `base_param_names` come from `model_spec(reference, "reference", n_chem)$params`
  (i.e. `setdiff(spec$params, spec$extra)`), so the same set is fixed in every
  model.

### Interaction-parameter start values (in `fit_model()`)

Change the default initialisation of the interaction params so that, when not
supplied in `start`:

- `a` starts at **0** (already the case).
- every `b`-family param (`b`, `b1`, `b2`, `b3`) starts at **1** (currently 0).

Today `fit_model()` zero-initialises the whole vector; this introduces a
per-name default for the `b`-family. Multistart jitter and `parscale` for the
interaction params continue to work (additive jitter around the start; `parscale`
floor of 1).

### `df` semantics

Each fit's reported `df` now counts only its free interaction params (stage-2
free params). LR tests use differences so are unchanged; `result_table()` will
display interaction-param counts. Document this in the `analyse_mixture()` /
`fit_model()` docs.

## Components changed

- `R/analyse.R` — `analyse_mixture()` becomes a two-stage orchestrator; add
  `fit_curve_from_singles()`; keep `seed_from_singles()` as fallback.
- `R/fit.R` — interaction-param default start values (`a=0`, `b*=1`).
- Docs/man pages regenerated for the above.

## Testing

New tests:

- Stage 1 on single-compound rows recovers the same curve params as a direct
  `fit_single()` per compound (shared `max` consistent with the marginals).
- Staged `analyse_mixture()` holds the base params **exactly equal** across all
  four fits (they are the stage-1 values, untouched by stage 2).
- Interaction start values: `b` starts at 1, `a` at 0 (assert via a fit with
  `n_starts = 1` and a 0/low-iteration budget, or by inspecting the assembled
  start vector).
- Fallback path: a design with a sparse single-compound series still returns a
  usable curve estimate without erroring.

Update existing expectations:

- `test-validation-binary`, `test-nchem`, `test-summary` — re-baseline to the
  staged results.

## Non-goals / YAGNI

- No `fit = 'staged'|'joint'` strategy flag on `analyse_mixture()`.
- No change to ternary model selection (SA-only vs DR/DL).
- No new public API surface.

## Outcome / as-built (2026-05-31)

Implemented on branch `mixdra-engine` (pushed). Full suite green: 205 pass / 0 fail.

**Scientific result — staged fitting changes real-data conclusions, and that is correct.**
On the binary MPs+CPF workbook data:

| | joint (old / workbook) | staged (now) |
|---|---|---|
| continuous reference SS | 1,633,770 | 2,224,216 |
| continuous chosen | DL | DL (unchanged) |
| continuous CA-vs-S/A χ² (p) | ~13.2 (3e-4) | 5.0 (0.025) |
| quantal reference deviance | 184.4 | 231.4 |
| quantal chosen | reference | **DR** |

The quantal case flips `reference`→`DR` (robust across seeds; S/A p≈4e-8, DR-vs-S/A
p≈3e-5). **Why it's correct, not a bug:** the Excel/MixTox workbook (and the old
joint code that reproduced it) lets the *mixture* data reshape the single-compound
curves, so interaction is absorbed into bent curves and hidden. Staged fitting
fixes the curves from the single-compound data only, so a genuine dose-ratio
interaction surfaces as `b` instead. Investigated and ruled out the marginal-fit
method as the cause: deviance-fit vs least-squares-fit marginals are essentially
identical (e.g. ec50₁ 4.635 vs 4.650); the divergence is purely that no
singles-only curve can reach the workbook's 184 mixture deviance (any such curve
gives ~231) — only a mixture-informed curve can, confirming the workbook is not
staged. The MixTox Shiny source confirms it fits curve params on treatment/mixture
data and never fits `a`/`b` at all (placeholders).

**Test design under staged fitting.** On noiseless synthetic data the
reference/SA objectives are ~1e-15, which makes the likelihood-ratio *model
selection* numerically degenerate (it spuriously picks a deviation with a≈±1e-8).
So the noiseless tests (`test-staged`, `test-nchem` ternary, `test-analyse`)
assert *parameter values* (`a≈0`, base params fixed-identical across fits, df
counts) rather than `chosen`. Model selection is exercised on the real-data
validation fixtures and on a small-noise known-interaction test (a=5 recovered).

**Post-review hardening (commit `1f3378d`).** Bound-name validation was being
skipped on the all-fixed fast path (so the staged reference fit silently ignored
an illegal bound the other fits reject); moved the check above the fast path.
Added a sparse-marginal fallback test and an all-fixed-bounds-error test.

**Man-page note.** `@keywords internal` functions still get `.Rd` files (just no
NAMESPACE export); `init_start_par.Rd` and `fit_curve_from_singles.Rd` exist and
are correct.

**Deferred follow-ups (not in scope here; surfaced to the user):**
- `param_ci()` is not staged-aware — it Hessians over the fixed base params, so it
  returns `NA` CIs for the staged reference fit. Pre-existing; needs a doc note or
  a staged-aware variant.
- A user-supplied *partial* `start` to `analyse_mixture()` (the fixed base) is not
  validated for completeness; a partial vector would silently leave some curve
  params free. The default path (`start = NULL`) is always safe.
- The generic `analyse_mixture()` still offers DR/DL for 3-chem input. The real
  ternary workflow uses the separate `analyse_ternary()` (Advanced S/A only), so
  this is a cosmetic cleanup left to the ternary subproject.
