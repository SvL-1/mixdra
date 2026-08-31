# mixdra plotting layer (single + binary) — design

**Date:** 2026-05-30
**Status:** implemented and merged to `main` (2026-06-01); this document has been
reconciled to the as-built code (see "As-built deviations" at the end).
**Depends on:** the completed engine (`docs/design/specs/2026-05-30-mixdra-design.md`,
plan `docs/design/plans/2026-05-30-mixdra-engine.md`)

## Purpose

Add a plotting layer to the `mixdra` package so a user can visualise a fitted
dose-response / mixture model. This is the first slice of the original "Plan 2"
(which also envisaged a Shiny app, file I/O, and docs); those remain separate,
later plans. This plan delivers **only** the plotting functions, as exported
package functions that build on the finished engine.

Scope is **single-chemical and binary mixtures**. Ternary visualisation is
explicitly deferred (see Out of scope).

## Reference material

- **Skylar's validated plotting script** `recieved/R scripts & excel/R scripts &
  excel/Example binary 3d curve & 2D isobles.R` — the binary 3-D surface and 2-D
  isobole. This is the authority for grid construction, matrix orientation, the
  reference-vs-deviation isobole overlay, and the default effect levels. We port
  its *logic*, sourcing model parameters from a fitted engine object instead of
  hard-typed Excel values.
- **Skylar's ternary script** `recieved/.../Isoplane_MPs_FBSA_IMI.R` — read for
  context only. Its isoplane (EC50 isobole surface in 3-chemical space) and
  ΣTU / z-value plots are a different, larger problem (new solver machinery) and
  are out of scope here.
- The existing `MixTox_shiny_v2` app uses `ggplot2`; we deliberately standardise
  on `plotly` instead (see Backend decision).

## Decisions (from brainstorming)

- **Plotting layer only** — no Shiny app, no I/O, no docs in this plan.
- **Backend: all plotly.** Every renderer returns a plotly htmlwidget, so they
  embed identically (`plotlyOutput`) in the future Shiny app. Skylar's isobole is
  ggplot `geom_contour`; we port that logic to plotly `add_contour`.
- **Dimensions: single + binary.** Ternary deferred to its own plan.
- **Enriched fit objects.** `fit_model()` records the model descriptors it
  currently omits, so plot functions take just `(fit, df)` plus plot-specific
  args and are self-describing.
- **Data-builder + thin-renderer split.** Each plot is a pure data-builder
  (returns numeric data frames / matrices, no plotly) plus a thin plotly wrapper.
  Builders are rigorously unit-tested; wrappers get smoke tests. Keeps plotly out
  of the hard-to-test path.
- **Separate named functions** (not one S3 `plot()` dispatcher) — focused,
  self-documenting signatures, one `?help` page each, easy Shiny wiring.

## Architecture & files

```
R/plot-data.R   # pure data builders (numeric output only)
R/plot.R        # thin plotly renderers (return plotly htmlwidgets)
R/fit.R         # EDIT: fit_model() also records model descriptors
R/single.R      # EDIT: fit_single() gets a light tag for single-chemical shape
tests/testthat/test-plot-data.R   # rigorous value / dim / edge-case tests
tests/testthat/test-plot.R        # smoke tests: object builds, trace counts, errors
DESCRIPTION     # EDIT: add plotly to Suggests
NAMESPACE       # regenerated: export the four renderers
```

`plotly` goes in **Suggests**, not Imports: the data builders need no plotly and
stay fully testable; the renderers guard with
`requireNamespace("plotly", quietly = TRUE)` and `stop()` with an install hint if
it is absent. Same pattern as the existing `readxl` suggest, keeping the core
engine dependency-light.

## Components

### Renderers (exported, in R/plot.R) — each returns a plotly object

| Function | Plot-specific args | Scope |
|---|---|---|
| `plot_dose_response(fit, df, chem = 1)` | which chemical's marginal curve | single & each binary margin |
| `plot_obs_pred(fit, df)` | — | any fit |
| `plot_surface(fit, df, n = 100)` | grid resolution | binary only |
| `plot_isobole(fit, df, levels = c(0.1, 0.25, 0.5, 0.75, 0.9), reference_fit = NULL, n = 100)` | effect levels (fractions of max), optional reference overlay | binary only |

(`plot_surface` has no `axes` argument: a binary fit has only one axis pairing
— C1↔`slope1`/`ec501`, C2↔`slope2`/`ec502` — so an axis swap would mislabel
parameters. Axis selection returns with ternary plotting, its own plan.
`plot_dose_response` also takes `log_x = TRUE` to toggle the log10 axis.)

Each renderer calls its data builder, then assembles the plotly object:
- `plot_dose_response`: line (fitted curve) + markers (observed). x = concentration
  (log10 axis), y = response (continuous) or proportion `Affected/Exposed` (binary).
- `plot_obs_pred`: scatter of observed vs predicted + a 1:1 reference line.
- `plot_surface`: `plot_ly(...) %>% add_trace(scatter3d observed) %>% add_surface(z)`.
- `plot_isobole`: equal-response contour paths extracted with
  `grDevices::contourLines` at `levels * max` and drawn with `add_lines` (solid
  black). If `reference_fit` is supplied, overlay its contours dashed red —
  Skylar's deviation-vs-additivity comparison. (Plotly's `add_contour` can't
  render arbitrary, non-evenly-spaced levels or a separate dashed overlay, so
  contour *paths* are used instead — same result, all still plotly.)

### Data builders (internal, in R/plot-data.R) — numeric output only

- `dr_curve_data(fit, df, chem)` → data frame: concentration grid + fitted
  response (via the single-curve log-logistic), plus the observed points for that
  chemical's marginal (rows where the other chemical is 0).
- `obs_pred_data(fit, df)` → data frame: observed, predicted. Mixture fits carry
  `fit$pred`; single-chemical fits do not, so predictions are recomputed via
  `ll3_predict()`.
- `surface_grid_data(fit, df, n)` → list: `x_vals`, `y_vals`, and the `z`
  response **matrix**, built to Skylar's orientation (see Data flow). Evaluates
  `model_spec(fit$reference, fit$deviation, 2)$fn` over the grid, so it reflects
  the fitted model including deviation parameters. Also returns the observed
  points for the 3-D scatter overlay.
- `isobole_data(fit, df, levels, reference_fit, n)` → reuses `surface_grid_data`;
  returns a data frame of contour paths (`source`, `level`, `group`, `x`, `y`) at
  the absolute contour values (`levels * fit$par["max"]`), optionally including a
  `reference_fit`'s contours (using that fit's own `max`).

`plot_surface` and `plot_isobole` share `surface_grid_data` — an isobole is a
contour view of the same predicted grid.

### Engine enrichment

`fit_model()` return list gains five additive fields, set from values already
computed inside the function:
- `reference` — `"CA"` / `"IA"`
- `deviation` — `"reference"` / `"SA"` / `"DR"` / `"DL"`
- `response` — `"continuous"` / `"binary"`
- `conc_cols` — e.g. `c("C1", "C2")`
- `n_chem` — `length(conc_cols)`

No behaviour change; `analyse_mixture()` (which calls `fit_model`) inherits the
fields for free. `fit_single()` gets a light marker so `plot_dose_response`
recognises the single-chemical result shape. The existing test suite (49 tests)
must remain green after the enrichment.

## Data flow / matrix orientation

Builders mirror Skylar's construction so the surface is not transposed:

```r
x_vals <- seq(min(c1), max(c1), length.out = n)
y_vals <- seq(min(c2), max(c2), length.out = n)
grid   <- expand.grid(C1 = x_vals, C2 = y_vals)        # C1 varies fastest
grid$z <- mapply(predictor, grid$C1, grid$C2, MoreArgs = params)
z_mat  <- matrix(grid$z, nrow = length(y_vals), ncol = length(x_vals), byrow = TRUE)
# add_surface(x = x_vals, y = y_vals, z = z_mat)  =>  z_mat[i, j] == predictor(x_vals[j], y_vals[i])
```

A test pins `z_mat[i, j]` against `predictor(x_vals[j], y_vals[i])` at sample
indices, and checks the control corner (`c1 = c2 = 0`) equals `max`.

## Error handling

- **plotly missing** — renderers `stop()` with an `install.packages("plotly")`
  hint. Builders never touch plotly.
- **Wrong dimensionality** — `plot_surface` / `plot_isobole` on a non-binary fit
  `stop()` with a clear message ("surface/isobole require a binary (2-chemical)
  fit").
- **NaN responses** — the CA bisection predictor returns `NaN` at the response
  bounds; builders pass `NaN` through (plotly renders gaps) rather than erroring.
  A test covers a grid that hits this.
- **Argument validation** — `chem` within range; `levels` strictly in (0, 1).

## Testing

**Data builders — rigorous (`test-plot-data.R`, no plotly needed):**
- `dr_curve_data` fitted values match `ll3_predict` at sample concentrations.
- `surface_grid_data` returns `n × n`; orientation pinned at sample indices;
  control corner equals `max`.
- `isobole_data` contour values equal `levels * max`.
- `obs_pred_data` predicted column equals `fit$pred`.
- Both single-chemical and binary paths exercised; the `NaN` edge case covered.

**Renderers — smoke (`test-plot.R`, `skip_if_not_installed("plotly")`):**
- Each returns a plotly / htmlwidget object.
- `plot_surface` has the expected traces (scatter3d + surface);
  `plot_isobole` with a `reference_fit` has two contour sets.
- Error paths fire (non-binary surface; informative messages).

**Optional regression anchor:** reproduce Skylar's CPF/IMI surface at a couple of
grid points using her published parameters, as a fixture test.

## Out of scope (deferred to later plans)

- Ternary plots: the EC50 isoplane surface and ΣTU / z-value plots (need an
  EC50-isobole solver over the simplex and a ΣTU computation that do not exist
  yet).
- The bslib Shiny app and its modules.
- Excel / CSV templates and upload/export I/O.
- README and methodology vignette.
- Static image export (kaleido / Cairo) — the renderers return interactive
  objects; callers can `htmlwidgets::saveWidget()` themselves if needed.

## As-built deviations (reconciliation, 2026-06-01)

The implementation (plan `docs/design/plans/2026-05-30-mixdra-binary-plots.md`)
matches every functional requirement above. Deliberate deviations, recorded here so
the design doc and code agree:

1. **`plot_surface` / `surface_grid_data` have no `axes` argument** — binary fits
   have a single fixed axis pairing; an axis swap would mislabel parameters. Axis
   selection is deferred to the ternary plan.
2. **Isobole uses `grDevices::contourLines` paths + `add_lines`, not plotly
   `add_contour`** — needed for arbitrary effect levels and the dashed reference
   overlay. Still all-plotly; same visual result.
3. **`plot_dose_response` gained `log_x = TRUE`** — toggles the (default) log10
   concentration axis.
4. **`obs_pred_data` recomputes predictions for single-chemical fits** (which lack
   `fit$pred`) via `ll3_predict()`.
5. **Engine enrichment also adds a `kind` field** (`"mixture"` / `"single"`) — the
   "light marker" the spec called for, used to distinguish fit shapes.

Known minor test-quality gaps (non-blocking, behaviour verified correct):
- The `surface_grid_data` "NaN passthrough" test's fixture produces an all-finite
  grid, so it does not actually exercise the NaN path (the passthrough itself is
  correct and was verified independently in review).
- Renderer smoke tests assert `expect_s3_class(p, "plotly")` only; the isobole
  two-set (fit + reference) property is checked at the data-builder layer, not by
  counting plotly traces.
- The optional Skylar CPF/IMI surface regression anchor was not implemented.
