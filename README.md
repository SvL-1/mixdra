# mixdra

**Mixture dose–response analysis in R — the Jonker et al. (2005) method, made reproducible.**

`mixdra` fits single-chemical and **binary/ternary mixture** dose–response models, detects and
quantifies how mixtures *deviate* from additivity (synergism, antagonism, dose-ratio- and
dose-level-dependent interactions), and tells you **which interaction pattern the data actually
support** via likelihood-ratio testing. It ships with a point-and-click Shiny app for
non-programmers and a clean scriptable API for everyone else.

It replaces a manual Excel/VBA + GRG-Solver workflow — the kind that lives in one spreadsheet on
one laptop and is almost impossible to reproduce or peer-review — with tested, version-controlled,
open-source code.

---

## Why this exists

The standard reference for mixture-toxicity interaction analysis — **Jonker, Svendsen, Bedaux,
Bongers & Kammenga (2005)**, *Environmental Toxicology and Chemistry* 24(10):2701–2713
([doi:10.1897/04-431R.1](https://doi.org/10.1897/04-431R.1)) — is widely cited, but in practice
it is run through a hand-driven Excel workbook: the model equations live in VBA, and the fitting
is done by clicking Excel's Solver, one model at a time, by hand.

That means every analysis is:

- **hard to reproduce** — the result depends on which cells were selected and which Solver run was
  accepted;
- **hard to trust** — there is no test suite, no version history, no record of the starting values;
- **hard to share** — it is one `.xlsm` file, not a tool.

`mixdra` takes the *same science* and makes it a proper instrument: a pure computational engine
(no spreadsheet, no GUI required) with a regression test suite pinned against the original Excel
results, wrapped in a modern Shiny dashboard.

## What it does

For **1, 2, or 3 chemicals**, `mixdra`:

- fits the per-chemical **log-logistic dose–response curves**;
- fits the two standard additivity **reference models** — **Concentration Addition (CA)** and
  **Independent Action (IA)**;
- fits the **deviation models** layered on top of each reference:
  - **S/A** — overall synergism / antagonism,
  - **DR** — dose-ratio-dependent deviation (binary),
  - **DL** — dose-level-dependent deviation (binary),
  - **Advanced S/A** — per-corner interaction for ternary mixtures;
- runs **likelihood-ratio tests** between each nested pair and flags the **most parsimonious
  model** — the simplest description the data don't reject;
- handles both **continuous** responses (sum-of-squares objective) and **binary / quantal**
  responses (binomial likelihood, *exposed* / *affected* counts).

## Features

- **One engine, three scales.** Single, binary, and ternary mixtures all run through a single
  n-agnostic predictor — no copy-pasted per-case model code to drift out of sync.
- **Automatic fit-and-compare.** `analyse_mixture()` fits the whole reference → S/A → {DR, DL}
  chain, builds the LR-test comparison table, and selects the winning model in one call.
- **Staged maximum-likelihood fitting.** Curve parameters are identified from the single-compound
  rows, then only the interaction parameters are fit to the mixture data — robust, fast, and
  faithful to how the method is meant to work.
- **A real optimiser, not a spreadsheet.** Base-R `optim` (L-BFGS-B with a Nelder-Mead fallback
  for the non-smooth CA bisection surface), **multi-start** seeding, hard parameter **bounds**,
  and the ability to **fix / pin** individual parameters — plus a "joint refine all parameters"
  mode that mirrors the Excel single-SSR fit when you want a digit-for-digit comparison.
- **Rich diagnostics & plots.** Dose–response curves, observed-vs-predicted, 3-D response
  surfaces, 2-D isoboles/contours, ternary **isoplanes** and **σ–TU** curves.
- **Validated.** A `testthat` suite checks the engine against the original Excel workbooks *and*
  against synthetic data with known parameters (params-in == params-out recovery on noise-free
  data).
- **Lean.** The engine depends only on `stats` and `numDeriv`; the Shiny UI stack is optional
  (Suggested) and only needed if you launch the app.

## Installation

```r
# install.packages("remotes")
remotes::install_github("SvL-1/mixdra")
```

To use the interactive app, also install the (Suggested) UI packages:

```r
install.packages(c("shiny", "bslib", "plotly", "DT"))
```

## Quick start

### The app (no coding)

```r
mixdra::run_app()
```

Opens a `bslib` dashboard with **Single Chemical** and **Binary Mixture** workflows: upload your
data (or use the built-in example), fit the curves, and step through fit-and-compare with live
plots and an editable data table.

### Scripted analysis

```r
library(mixdra)

# `df` has C1, C2 concentration columns plus either `Res` (continuous)
# or `Exposed` / `Affected` (quantal).
result <- analyse_mixture(df, reference = "CA", response = "continuous")

result$chosen        # name of the selected model, e.g. "CA_DR"
result$comparison    # likelihood-ratio test table (χ², df, p-value vs parent)
result$fits          # every fitted model, with parameters and fit statistics
```

Single chemicals and ternary mixtures use the matching entry points:

```r
analyse_single(conc, resp)               # one log-logistic curve
analyse_ternary(df, reference = "CA")    # 3-chemical fit + Advanced S/A
```

Need synthetic data to sanity-check a design before you run the experiment?

```r
sim <- simulate_mixture(reference = "CA", deviation = "DR", ...)  # generate, then recover
```

## Outputs

Every analysis returns the parameter block you'd expect from the paper and the Excel workbook —
fitted `max`, per-chemical slopes and EC50s, the interaction parameters (`a`, `b`/`b_i`),
the objective value (SS or −2·logL), and the LR-test comparison — plus tidy result tables
(`result_table()`, `ternary_effect_table()`) and the full plotting family
(`plot_dose_response()`, `plot_obs_pred()`, `plot_surface()`, `plot_isobole()`,
`plot_isoplane()`, `plot_sigma_tu()`).

## The method, briefly

Each chemical is described by a three-parameter log-logistic curve. **CA** combines chemicals on a
shared toxic-unit (concentration) scale; **IA** combines them on a probabilistic (effect) scale.
A single interaction term is then added to the reference model to capture a specific *pattern* of
non-additivity, and nested models are compared by likelihood ratio so the reported interaction is
the one the data justify — not the most flattering one. See Jonker et al. (2005) for the full
derivation.

## Validation & testing

The engine is checked two ways:

1. **Against the source workbooks** — `mixdra` is expected to select the same best model and land
   on parameters consistent with the published Excel/VBA + Solver results (fixtures extracted from
   the original binary and ternary workbooks).
2. **Against the truth** — synthetic datasets are generated with known parameters and noise-free,
   then re-fit; the engine must recover the generating parameters.

```r
# run the suite (or a single file)
devtools::test()
devtools::test(filter = "recovery")
```

## Scope

**In scope (v1):** single + binary + ternary mixtures; CA and IA references; S/A, DR, DL, and
Advanced S/A deviations; continuous and quantal responses; run-locally Shiny app.

**Out of scope (for now):** hosted/click-a-link web deployment; more than 3 chemicals;
non-monotonic (J-shaped / hormesis) curves.

## Citation

If you use `mixdra`, please cite the underlying methodology:

> Jonker, M.J., Svendsen, C., Bedaux, J.J.M., Bongers, M. & Kammenga, J.E. (2005). Significance
> testing of synergistic/antagonistic, dose level–dependent, or dose ratio–dependent effects in
> mixture dose–response analysis. *Environmental Toxicology and Chemistry*, 24(10), 2701–2713.
> https://doi.org/10.1897/04-431R.1

## License

GPL-3.
