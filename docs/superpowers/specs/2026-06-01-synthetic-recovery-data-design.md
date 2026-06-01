# Synthetic Recovery Data — Design

**Date:** 2026-06-01
**Status:** Approved (pending spec review)
**Branch:** mixdra-engine

## Purpose

Provide an exported, deterministic **forward generator** that produces mixture
dose-response datasets from *known* parameters, so the engine's fitting and
model-selection can be verified by recovery: feed the engine data generated from
a known model and confirm it recovers the parameters that produced it.

Priorities, in the user's order:

1. **Parameter recovery (primary).** At zero noise, params in == params out
   within optimiser tolerance. A true correctness oracle, not a regression check.
2. **Simplicity.** A test with no unknowns and no randomness — easy to run,
   trivial to reason about.
3. **Noise robustness (deferred).** Adjustable noise to characterise *when*
   recovery/selection degrades. Ships after recovery is proven.
4. **Public verifiability.** Exported functions make this a reproducible,
   inspectable method suitable for GitHub.

### Why this engine especially benefits

The engine is an inverse problem: it recovers log-logistic curves
(`max`, `slope*`, `ec50*`) plus an interaction deviation (CA/IA reference +
SA/DR/DL with `a`, `b`). The natural oracle for an inverse problem is the forward
model run with known inputs.

The current validation oracle is weak: tests assert magic numbers
(`res$fits$reference$objective == 2224215.9`) that detect regressions but cannot
confirm correctness, and which intentionally differ from the published workbook
so there is no external ground truth. Recovery testing supplies that ground
truth.

## Architecture

One model definition, used both directions — the generator goes forward, the
fitter goes inverse, and they share a single predict path so they cannot drift.

```
                 ┌─────────────────────────────┐
   known par ───▶│  mixture_predict(df, par,    │───▶ μ (noise-free expected)
   design grid   │     reference, deviation)    │      │
                 └─────────────────────────────┘      │
                         ▲  (shared)                   ▼
                         │                      ┌──────────────┐
   fit_model() ──────────┘                      │ noise layer  │──▶ observed data
   (inverse)                                    └──────────────┘     │
        ▲                                                            │
        └──────────────── analyse_mixture() ◀───────────────────────┘
                          recovers par   →  compare to known par
```

The fitter's entire forward evaluation today is (fit.R:59-72):

```r
spec <- model_spec(reference, deviation, n_chem)
conc <- setNames(as.list(df[conc_cols]), paste0("c", seq_along(conc_cols)))
do.call(spec$fn, c(conc, as.list(p_full)))
```

`model_spec` already dispatches every `reference × deviation × n_chem` combo, so
binary and ternary fall out of the same path with no special-casing.

## Components

### New files

| File | Contents |
|------|----------|
| `R/predict-mixture.R` | `mixture_predict(df, par, reference, deviation)` — the shared forward helper, extracted from `fit_model`. |
| `R/simulate.R` | Exported `simulate_mixture()`, `simulate_single()`, `mixture_design()`. |
| `tests/testthat/test-simulate.R` | Generator unit tests. |
| `tests/testthat/test-recovery.R` | Round-trip recovery + selection across the model grid (the payoff). |

Plus `man/` pages and NAMESPACE entries for the three exported functions.

### `mixture_predict(df, par, reference, deviation)` — internal

- Infers `n_chem` from the `C1..Cn` columns present in `df`.
- Calls `model_spec(reference, deviation, n_chem)` and applies `spec$fn` over the
  concentration columns with the supplied full parameter vector.
- Returns a numeric vector of noise-free expected responses on the model's
  natural scale (continuous = response value; binary = probability in (0,1)).
- `fit_model` is refactored so its inline `predict_with` becomes a thin closure
  over `mixture_predict`. **Acceptance:** the full existing test suite passes
  unchanged after the refactor (no behaviour change to `fit_model`).

### `simulate_mixture()` — exported

```r
simulate_mixture(par, reference = "CA", deviation = "reference",
                 response = c("continuous", "binary"),
                 design = NULL,            # default: mixture_design(par, ...)
                 cv = 0, group_size = Inf, # the two noise knobs
                 reps = 1, seed = NULL)
```

- `par` — named full parameter vector matching
  `model_spec(reference, deviation, n_chem)$params`, e.g.
  `c(max, slope1, slope2, ec501, ec502, a, b)`. Validated against the spec's
  `params` (error on missing/extra names). `n_chem` inferred from the number of
  `slope*` entries.
- `design` — a data frame of `C1..Cn`. When `NULL`, built by `mixture_design()`.
- Computes `mu <- mixture_predict(design, par, reference, deviation)`, then
  applies the noise layer.
- Returns a data frame the engine consumes unchanged: `C1..Cn` plus `Res`
  (continuous) or `Exposed`/`Affected` (binary). With `reps > 1` each design row
  is repeated `reps` times before the noise layer.

### `simulate_single(curve, conc = NULL, response, cv, group_size, seed)` — exported

- `curve = c(max, slope, ec50)`. Convenience for the 1-chemical case.
- Returns `C1` + `Res` (or `Exposed`/`Affected`) for `analyse_single`/`fit_single`.
- Default `conc`: geometric ladder around `ec50` (see design grid).

### `mixture_design(par, ratios = NULL, n_per_ray = 7, ...)` — exported

Takes `par` so concentrations can be placed around each chemical's `ec50`
(otherwise `slope`/`ec50` are not identifiable from the data):

- 1 control row (all concentrations 0).
- Per-chemical **single-compound ladder**: geometric, `ec50ᵢ · 2^(-3:3)`, all
  other chemicals 0. These rows feed the staged curve recovery.
- **Mixture rays**: fixed-ratio rays, each a geometric dose ladder.
  - Binary default ratios (in toxic units): 1:1, 1:3, 3:1.
  - Ternary default: centroid 1:1:1 plus the three edge midpoints.
- Exported and overridable so a grid can be inspected or hand-crafted.

## Noise layer

Two knobs, one per response type:

- **Continuous:** `Res = μ + rnorm(n, 0, cv * abs(μ))`. Relative (CV) Gaussian,
  portable across any `max` scale. `cv = 0` → `Res = μ` exactly.
- **Quantal:** `group_size` is the knob.
  - `group_size = Inf` (default) → **exact**: `Exposed = 100`,
    `Affected = 100 · μ` (fractional). `binlik` (objective.R:18) uses only
    `pi = affected/exposed` and logs — no factorials — so fractional `Affected`
    is mathematically exact and `Affected/Exposed = μ` recovers with zero
    quantisation error.
  - `group_size = N` finite → `Affected ~ Binomial(N, μ)`. Small N = more
    scatter; large N → exact. Statistically honest binomial noise.

## Determinism

`seed` seeds the RNG locally for the call. With `cv = 0` and `group_size = Inf`
the output is fully deterministic regardless of `seed` — the headline recovery
test contains no randomness at all.

## Tests

### `test-simulate.R` — generator unit tests

- `mixture_design` produces the expected control + single ladders + mixture rays
  for binary and ternary; concentrations span each `ec50`.
- At `cv = 0` / `group_size = Inf`, the generated `Res` / `Affected/Exposed`
  equal `mixture_predict` exactly.
- Noise scaling sanity: over many seeded continuous draws, residual SD ≈ `cv·|μ|`;
  over many seeded binomial draws, mean `Affected/Exposed ≈ μ`.
- `par` validation rejects mismatched names.

### `test-recovery.R` — the payoff

For `reference ∈ {CA, IA}` × `deviation ∈ {reference, SA, DR, DL}` ×
`response ∈ {continuous, binary}`:

1. Generate zero-noise data from a known `par` (interaction params chosen clearly
   non-trivial and well-separated from sibling models).
2. Run `analyse_mixture`.
3. Assert:
   - curve params (`max`, `slope*`, `ec50*`) recovered within ~1e-3 relative;
   - interaction params (`a`, `b`) recovered within tolerance;
   - `chosen == injected deviation` (selection works).

Ternary: **reference-level recovery only** — CA/IA reference and the SA/DR/DL
deviations where implemented. ASA (Advanced S/A) is skipped: its fitting is on
hold, so its round-trip cannot be validated yet.

## Out of scope (deferred)

- **Noise sweep (goal #3).** A characterisation of recovery/selection vs
  increasing `cv` / decreasing `group_size`. Added on top once exact recovery is
  trusted. The noise knobs themselves ship now; the sweep harness/report does not.
- **Vignette.** A public runnable recovery/sweep document. Easy to add later;
  not in this iteration.
- **Ternary ASA recovery.** Blocked on the Advanced S/A deviation-model hold.

## Caveats (explicit, not hidden)

1. **Shared-model circularity.** Because the generator and fitter share
   `mixture_predict`, a bug *inside the model formula itself* is not caught by
   recovery — both sides would be wrong identically. The workbook fixtures remain
   the external truth-guard. Recovery tests the optimiser, objective, staging,
   and selection — not the formula.
2. **Too-easy data.** Model-generated noise-free data is, by construction, the
   easiest possible case. It complements, never replaces, the real workbook
   fixtures.
3. **Identifiability.** DR and DL siblings can be ambiguous for some `par`;
   selection tests use clearly-separated injected interactions. A tiny injected
   interaction may (correctly) cause the engine to select a simpler model.
4. **Fractional `Affected`.** The exact quantal representation yields non-integer
   `Affected`, which is mathematically fine for `binlik` but cosmetically unusual
   if such a frame is printed or plotted as raw counts.
