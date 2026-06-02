# Model-unification migration: wire production through `mix_response()`

Date: 2026-06-02
Status: ✅ IMPLEMENTED (2026-06-02) — merged to `main` (branch `model-unification-migration`,
commits `20be962`..`033cc9d`). Net −1,539 lines in `R/`; 764 assertions green incl. the 576
equivalence contract; one ~5e-8 tolerance re-baseline; ASA kept separate. See plan
`docs/superpowers/plans/2026-06-02-model-unification-migration.md`.
Supersedes the "Step 2" left undone by `2026-06-01-model-unification-design.md`.

## Problem

Commit `cd10b00` introduced `R/mix-response.R`, a single n-agnostic predictor that
provably replaces the 16 verbatim per-`(reference × deviation × n_chem)` model
functions. Step 1 (write the unified predictor and prove equivalence) is done and
pinned by `tests/testthat/test-mixture-predict.R` (576 assertions). Step 2 — wiring
production to it and deleting the originals — was never completed.

As a result the package ships **two complete implementations** of the binary/ternary
model family:

| Source | Lines | Status |
| --- | --- | --- |
| `R/models-binary.R` (8 fns + Vectorize wrappers) | 408 | live in production |
| `R/models-ternary.R` (8 fns + Vectorize wrappers) | 1184 | live in production |
| `R/mix-response.R` (the replacement) | 81 | unused except by its own test |

`registry.R::model_spec()` resolves the old functions by string key
(`get("ca_dl_bi_vec")`); `mix_response()` is referenced only by `test-mixture-predict.R`
and design docs. The two are kept in lockstep by hand. This is the largest
maintenance liability in the codebase.

## Goal

Make `mix_response()` the single production predictor for the 16 SA/DR/DL/reference
models, remove the ~1,592 duplicated lines from the shipped package, and change **zero**
downstream API behavior (beyond accepted last-digit numeric shifts). The ternary
Advanced-S/A (ASA) model is explicitly **out of scope** and stays as-is.

## Decisions (locked)

1. **Old functions become a test-only oracle.** Relocate all 16 verbatim functions
   and their `Vectorize` wrappers into `tests/testthat/helper-legacy-models.R`.
   `test-mixture-predict.R` keeps comparing `mix_response()` against them dynamically,
   so the live equivalence proof survives — but the shipped package sheds the lines.
2. **Accept re-baselining.** Keep `mix_response()`'s 200-iteration / `1e-10` bisection.
   Fitted `a`/`b` and reported SSR/R²/F/P on the reference data shift in the last
   digits (agreement to ~`1e-4`, not bit-identical). Re-run reference fits, update any
   golden-value assertions, note the shift in the commit message.
3. **Defer ASA.** This migration touches only the 16. `registry.R` keeps resolving
   `ca_asa_tri_vec` for the `ASA` branch unchanged. Folding ASA into `mix_response()`
   is a follow-up once the 16→1 switch is proven in production.

## The single change point: `registry.R`

`model_spec()` is the only seam. Every production consumer obtains its predictor from
`model_spec(...)$fn`:

- `fit.R:59,68` (`.mixture_eval`)
- `predict-mixture.R:36,42` (`mixture_predict` / `.mixture_eval`)
- `plot-data.R:72,74` (`predict_grid`)
- `summary.R:14,22` (`param_ci`, via `do.call(spec$fn, args)` **by name**)
- `simulate.R` (via `mixture_predict`)

Only what `$fn` resolves to changes:

- `reference` / `SA` / `DR` / `DL` → an **adapter closure** wrapping `mix_response()`.
- `ASA` → unchanged, `get("ca_asa_tri_vec")`.

`$params`, `$extra`, `$parent` are untouched. No consumer changes.

## The adapter

The call contract that must be preserved (from `.mixture_eval` and `param_ci`):

```r
do.call(fn, c(conc, as.list(par)))
# conc = list(c1=, c2=[, c3=])   — concentration vectors, one per chemical
# par  = named flat vector: max, slope1, slope2[, slope3], ec501, ec502[/ec50_*], a, b[1..3]
```

The adapter accepts that exact by-name signature and repacks per row into
`mix_response(concs, max, slopes, ec50s, reference, deviation, a, b)`:

```r
make_adapter <- function(reference, deviation, n_chem) {
  slope_names <- if (n_chem == 2) c("slope1","slope2")
                 else            c("slope1","slope2","slope3")
  ec50_names  <- if (n_chem == 2) c("ec501","ec502")
                 else            c("ec50_1","ec50_2","ec50_3")
  cc_names    <- paste0("c", seq_len(n_chem))

  function(...) {
    A      <- list(...)
    concs  <- lapply(cc_names, function(k) A[[k]])      # one vector per chemical
    slopes <- vapply(slope_names, function(k) A[[k]], numeric(1))
    ec50s  <- vapply(ec50_names,  function(k) A[[k]], numeric(1))
    a      <- if (is.null(A$a)) 0 else A$a
    b      <- adapter_bmap(deviation, n_chem, A)
    n_pts  <- length(concs[[1]])
    vapply(seq_len(n_pts), function(i) {
      mix_response(vapply(concs, `[[`, numeric(1), i),
                   A$max, slopes, ec50s, reference, deviation, a, b)
    }, numeric(1))
  }
}
```

### `b` mapping (`adapter_bmap`)

`mix_response()` expects `b[active]` as a per-chemical vector for DR and a scalar for
DL. The flat `par` carries different shapes; the mapping (already pinned by
`test-mixture-predict.R`):

| deviation | n_chem | flat `par` carries | passed to `mix_response` as |
| --- | --- | --- | --- |
| reference | any | (no b) | `0` |
| SA | any | (no b) | `0` |
| DR | 2 | scalar `b` | `c(b, 0)` |
| DR | 3 | `b1, b2, b3` | `c(b1, b2, b3)` |
| DL | 2 or 3 | scalar `b` | scalar `b` (passthrough) |

The binary-DR `c(b, 0)` form matches the original `(a + b * z1)`, where `b` weights only
chemical 1's dose fraction.

## Files touched

- **`R/registry.R`** — add `make_adapter` + `adapter_bmap` (or place them in
  `R/mix-response.R` next to the predictor); change the `model_spec()` resolution so the
  16 branches return `make_adapter(...)` and ASA keeps `get(...)`.
- **`R/models-binary.R`** — deleted (moved to test helper).
- **`R/models-ternary.R`** — deleted (moved to test helper).
- **`R/models-ternary-asa.R`** — unchanged.
- **`tests/testthat/helper-legacy-models.R`** — new; the 16 functions + Vectorize wrappers.
- **`tests/testthat/test-mixture-predict.R`** — unchanged behavior; now compares against
  the helper-loaded oracle instead of package functions.
- **Golden-value tests** (e.g. reference-fit assertions) — re-baselined where last-digit
  shifts occur.
- **`man/`** — regenerate; the 16 `.Rd` pages for the deleted internal functions go away
  (they were `@keywords internal`, not exported, so no NAMESPACE change).

## Risks and how they are handled

- **`b`-mapping correctness** — the one fiddly transform. Covered directly by
  `test-mixture-predict.R`'s DR/DL cases; verify those pass before deleting anything.
- **Bisection numerics** — accepted (decision 2); re-baseline.
- **Skipped edge cases** — `test-mixture-predict.R` skips rows where an original returns
  `NULL`/`NaN` (e.g. mixed-sign slopes). The package only fits all-positive or
  all-negative slope regimes, so this is theoretical, but it means mixed-sign behavior is
  *unverified* rather than *proven-equivalent*. Note, do not block on it.
- **roxygen drift** — regenerate `man/` after deletion so no orphan `.Rd` pages remain.

## Verification gates

1. `test-mixture-predict.R` (576 assertions) green against the relocated oracle.
2. Filtered runs of the model / fit / predict-mixture / simulate / summary / plot tests.
3. Reference-fit re-baseline: re-run, confirm changes are last-digit only, update goldens.
4. Smoke check: `model_spec` dispatch returns finite predictions across a binary and a
   ternary grid for every `(reference, deviation)` combination, ASA included.

## Net effect

~1,592 lines removed from the shipped package, one ~30-line adapter added, no downstream
API change, ASA untouched, the equivalence oracle preserved as a test helper.
