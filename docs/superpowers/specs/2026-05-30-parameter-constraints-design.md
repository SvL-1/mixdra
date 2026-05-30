# Parameter constraints for mixdra fitting — design

Date: 2026-05-30
Status: approved (pending user review of this spec)
Related: [parameter-constraints-decisions.md](../../parameter-constraints-decisions.md)

## Problem

`fit_model()` derives its box constraints from the starting values via fixed
multipliers:

```r
lower[base] <- pmax(1e-8, theta0[base] * lower_frac)   # default 0.1x
upper[base] <- pmax(lower[base] * 1.01, theta0[base] * upper_mult)  # default 10x
```

Consequences:

- The feasible region depends on the seed, not on domain meaning. A parameter
  can be silently pinned to an artificial bound (e.g. EC50 Chem 1 = 100 in the
  first binary example was stuck at `seed * 10`, not fitted).
- There is no way to express a real biological/chemical constraint (e.g. a
  concentration cannot exceed 100%, or a slope should not exceed a plausible
  maximum).
- Deviation parameters (`a`, `b`, …) are fully unconstrained, which is correct
  for them but is incidental rather than intentional.

The original MixTox workflow constrains parameters to domain-meaningful ranges;
the engine should let the user do the same explicitly.

## Goals

1. Let the user set hard `[lower, upper]` bounds **per parameter**, by name.
2. Default to **positivity only** for base curve params when the user supplies
   no bound. Remove the seed-multiplier scheme.
3. Keep deviation parameters (`a`, `b`, `b1/b2/b3`) **unconstrained** and
   **not user-overridable**, to protect downstream interaction analysis (e.g.
   the concentration at which an interaction switches synergistic ↔ antagonistic).
4. Preserve the binary correctness rule: `max ≤ 1`.
5. Document the new interface in roxygen and a user-facing vignette.

## Non-goals

- No automatic/domain-aware *default* bounds beyond positivity (the user
  rejected smart defaults — constraints are manual).
- No change to the optimiser, model functions, multi-start strategy, or
  `analyse_mixture`'s model-selection logic.
- No new constraint types (only box bounds + the existing `fixed` pinning).

## Design

### API

`fit_model()` loses `lower_frac` / `upper_mult` and gains `lower` / `upper`:

```r
fit_model(df, reference, deviation = "reference",
          response = c("continuous", "binary"),
          start, fixed = character(0),
          lower = NULL, upper = NULL, n_starts = 1)
```

- `lower`, `upper`: optional named numeric vectors keyed by base parameter name
  (`"max"`, `"slope1"`, `"slope2"`, `"ec501"`, `"ec502"`, and the ternary
  `"slope3"`, `"ec50_1"…"ec50_3"`). Only base curve params are accepted.
- Names not present fall back to the defaults below. Unknown names (typos, or a
  deviation-param name like `"a"`) raise a clear error so silent no-ops can't
  hide a mistake.

### Bound resolution

For each **free** parameter:

| Group | Default lower | Default upper | User-overridable |
|---|---|---|---|
| Base curve params (`max`, `slope*`, `ec50*`) | `1e-8` | `Inf` | yes, via `lower`/`upper` |
| Deviation params (`spec$extra`: `a`, `b`, `b1/b2/b3`) | `-Inf` | `+Inf` | no |

Binary correctness exception: when `response == "binary"` and `max` is free,
`upper["max"]` defaults to `1`. A user-supplied `upper["max"]` may *tighten* it
(e.g. 0.98) but a value > 1 is clamped to 1 with a warning, since the model
yields probabilities > 1 above that.

Resolution order, per base param `p`:

```
lo[p] <- if (p in names(lower)) lower[[p]] else 1e-8
hi[p] <- if (p in names(upper)) upper[[p]] else (1 if binary & p=="max" else Inf)
hi[p] <- min(hi[p], 1) if binary & p=="max"   # with a warning if it bit
# validate lo[p] < hi[p], and start[p] within [lo, hi] (clamp start with a warning if not)
```

Deviation params keep `-Inf/+Inf` regardless of `lower`/`upper`.

### Touch points in `fit.R`

- Replace the `lower`/`upper` construction block (current lines ~60–68) with the
  resolution above. `base <- setdiff(free, spec$extra)` already isolates the base
  params, so deviation params are left at `±Inf` by construction.
- **Enforce bounds in the objective.** L-BFGS-B respects `lower`/`upper`
  natively, but the Nelder-Mead *fallback* does not — so a parameter could
  escape its bound whenever L-BFGS-B fails to converge and NM takes over.
  `obj_free` therefore returns `1e12` for any infeasible `theta`, making the box
  constraints bite regardless of which optimiser is active.
- Multi-start jitter (lines ~96–108) is unchanged: it clamps base params to
  `[lower, upper]` (a no-op when `upper = Inf`) and jitters deviation params
  additively as before.
- `parscale` logic unchanged.

### `analyse_mixture()`

Add `lower = NULL`, `upper = NULL` and forward them to every `fit_model()` call,
so one bound set applies across reference/SA/DR/DL (they share base params).

### Validation & errors

- `lower`/`upper` names must be a subset of the active base params → else
  `stop()` listing the offending names and the valid set.
- Each resolved `lower < upper` → else `stop()`.
- `start` outside `[lower, upper]` → clamp into range with a `warning()` (so the
  optimiser starts feasible).
- Binary `upper["max"] > 1` → clamp to 1 with a `warning()`.

## Documentation

### roxygen (regenerate `man/` via `devtools::document()`)

- `fit_model`: document `lower`/`upper` (named vectors, base params only,
  positivity default, binary `max ≤ 1`, deviation params always unconstrained).
  Remove `lower_frac`/`upper_mult` from `@param`.
- `analyse_mixture`: document the forwarded `lower`/`upper`.
- Update `DESCRIPTION` roxygen version stamp as in prior commits.

### Vignette (new)

Create `vignettes/constraints.Rmd` (and add `knitr`, `rmarkdown` to `Suggests`
plus `VignetteBuilder: knitr` in `DESCRIPTION`). Contents:

1. Why constrain — biological/chemical meaning (concentration ≤ 100%, plausible
   slope ceiling), and why `a`/`b` are deliberately left free.
2. Worked binary example reproducing the EC50 Chem 1 case: fit unconstrained,
   show the parameter sitting at a high value, then `upper = c(ec501 = 100)`
   and `fixed = "max"` with `start["max"] = 0.98`, and compare.
3. Reference table of the base parameter names per model size.

## Testing

- `lower`/`upper` named by base param produce the expected `optim` bounds
  (introspect via a small fit where the bound bites).
- Unknown / deviation-param names in `lower`/`upper` raise an error.
- Binary `max` upper defaults to 1; an explicit `upper["max"] = 0.98` is honoured;
  `upper["max"] = 1.5` clamps to 1 with a warning.
- Deviation params remain unconstrained even if (erroneously) named.
- A constrained `ec501 = c(upper = 100)` fit lands at ≤ 100 and is reproducible.
- Removing `lower_frac`/`upper_mult` doesn't break existing `analyse_mixture`
  default runs (positivity-only base bounds still converge on the example data).

## Migration

`lower_frac`/`upper_mult` are removed (not deprecated) since the package is
pre-release (`0.0.0.9000`) and they were internal tuning knobs. Existing callers
relying on defaults are unaffected except that base params are no longer capped
at `10×seed`; if a previous fit depended on that cap to avoid a runaway
parameter, the user now sets an explicit `upper`.
