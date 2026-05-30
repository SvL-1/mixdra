# mixdra — Mixture Dose–Response Analysis: Design Spec

**Date:** 2026-05-30
**Status:** Approved design (pre-implementation)

## 1. Purpose

`mixdra` makes the Jonker et al. (2005) mixture-toxicity methodology publicly available to
other scientists as an R package with a Shiny dashboard. It fits single-chemical
dose–response curves and mixture models for **1, 2, or 3 chemicals**, detects and quantifies
**synergism/antagonism (S/A), dose-ratio-dependent (DR), and dose-level-dependent (DL)**
deviations from the **concentration addition (CA)** and **independent action (IA)** reference
models, and identifies which interaction pattern best describes the data via likelihood-ratio
testing.

It replaces a manual Excel/VBA + GRG-Solver workflow (used in Sam's published paper) and an
incomplete work-in-progress Shiny app (`MixTox_shiny_v2`).

### Background / prior art reviewed
- **Jonker et al. 2005**, *Environ. Toxicol. Chem.* 24(10):2701–2713 — the methodology
  (single log-logistic curve Eq. 14; CA Eq. 2/3; IA Eq. 4/5; deviation functions S/A Eq. 7,
  DR Eq. 8/11, DL Eq. 12/13; ML fitting Eq. 16; likelihood-ratio comparison → Table 2).
- **`MixTox Model_binary_MPs_CPF.xlsm`** — binary workbook. VBA contains only the worksheet
  model functions (`CA`, `IA`, `CA_SA`, `CA_DR`, `CA_DL`, `IA_*`) and likelihood functions
  (`BinLik`, `TotalBinLik`). Optimisation is done by Excel's built-in **GRG Nonlinear Solver**,
  driven manually — there is no solver code in the VBA.
- **`FBSA CPF IMI ternary -simplified.xls`** — ternary workbook (validation fixture).
- **`MixTox_shiny_v2`** — existing R Shiny app. `functions/model_functions.R` is a faithful,
  complete R port of the VBA model functions for **both binary and ternary**, all models. The
  fitting/UI layers are incomplete (crude 2-step SBPLX least-squares only, continuous only,
  several noted bugs, tangled logic).

### Key reframing
The scientific math is small and already ported. The real engineering is: (1) a reliable
constrained ML optimiser, (2) clean glue + model-comparison statistics, (3) the UI, and (4)
plots. None of this requires matching Excel digit-for-digit.

## 2. Decisions (locked)

| Decision | Choice |
|---|---|
| Language / platform | **R + Shiny**, shipped as a proper R package on GitHub |
| Distribution | **Run locally** (`shiny::runGitHub(...)` or `install_github` + `mixdra::run_app()`). No server, no hosting, no webR/WASM constraints. |
| Scope | **Single + binary + ternary**, one shared engine |
| Fitting | **Auto-fit & compare**: fit Reference + S/A + DR + DL, run LR tests, flag most parsimonious; **advanced panel** for manual start values + fixing parameters |
| Reproducibility | **Scientifically equivalent** to Excel/paper (same best model, parameters within CIs); sanity-checked against the two Excel workbooks. Robust modern optimiser is acceptable. |
| Response types | **Continuous** (sum-of-squares objective) **and binary/quantal** (binomial objective, Eq. 16) |
| Package name | **`mixdra`** (`mixtox` is taken on CRAN in the same field) |
| Licence | GPL-3 (proposed; common in R ecosystem) |

## 3. Architecture — engine / UI separation

The current app tangles fitting logic into Shiny modules. `mixdra` separates a **pure
computational engine** (no Shiny) from a **thin Shiny layer**, so the engine is independently
testable, scriptable, and trustworthy.

```
mixdra/
├─ DESCRIPTION, NAMESPACE         # R package metadata
├─ R/                             # ENGINE (no Shiny):
│   ├─ models-single.R           #   log-logistic single dose–response
│   ├─ models-ca.R               #   CA reference + CA_SA / CA_DR / CA_DL (binary & ternary)
│   ├─ models-ia.R               #   IA reference + IA_SA / IA_DR / IA_DL (binary & ternary)
│   ├─ objective.R               #   SS (continuous) + binomial (binary) objectives
│   ├─ fit.R                     #   staged ML fitter, multi-start, bounds, fixed params
│   ├─ compare.R                 #   LR tests, df, p-values, most-parsimonious selection, CIs
│   ├─ plots.R                   #   DR curve, obs-vs-pred, 3D surface, isobole/contour
│   └─ io.R                      #   template read/write, result export
├─ inst/app/                      # SHINY app (run via mixdra::run_app())
│   ├─ app.R, ui.R, server.R
│   └─ modules/ (setup, single, binary, ternary, summary)
├─ inst/templates/                # downloadable Excel/CSV templates
├─ inst/extdata/                  # the 2 Excel datasets as validation fixtures
├─ tests/testthat/                # regression tests vs Excel results
├─ vignettes/                     # methodology recap (Jonker) + how-to
└─ README.md                      # install + quick start + screenshots
```

The model functions are adapted from the existing `model_functions.R` (already validated
against the VBA), refactored for clarity and consistency.

## 4. The fitting engine

### 4.1 Staged maximum-likelihood (mirrors paper & Excel)
1. **Single curves** — fit log-logistic per chemical → `max`, `slope`, `EC50`. Seeds mixtures.
   Replaces the `drc` dependency with a base-`optim` fit.
2. **Reference model** — refine on mixture data. CA uses the bisection solver from
   `model_functions.R`; IA is closed-form.
3. **Deviation models** — add S/A (`a`), DR (`a`, `b`/`b_i`), DL (`a`, `b_DL`) on top of the
   reference.

### 4.2 Objectives (swappable by response type)
- **Continuous:** minimise sum of squared residuals (≡ Gaussian ML).
- **Binary/quantal:** minimise the binomial objective of Eq. 16 (the `BinLik`/`TotalBinLik`
  logic), using exposed/affected counts.

### 4.3 Optimiser
- Backbone: **`bbmle::mle2`** — purpose-built for ML; yields `logLik`/AIC, profile-likelihood
  **95% CIs**, and `anova()` LR tests that map 1:1 onto the Table 2 output.
- Engine: `optim` **L-BFGS-B** (box constraints) with **Nelder-Mead** as the derivative-free
  fallback (robust to the non-smooth CA bisection inner solver).
- **Multi-start** from perturbed single-fit seeds to avoid local minima.
- **`DEoptim`** as a heavy-duty global fallback for stubborn fits.
- Parameter bounds and sign constraints per Table 1 (e.g. EC50 > 0; slope sign fixed by
  endpoint direction).

### 4.4 Model comparison
For each selected reference (CA and/or IA), fit {Reference, S/A, DR, DL}, then:
- compute χ² = difference in −2·logL between nested models, df = number of extra parameters,
  and the LR-test p-value;
- flag the **most parsimonious** model (simplest model not significantly improved upon),
  reproducing the paper's Table 2 logic.

## 5. UI / workflow (Shiny, `bslib`)

Modern `bslib` dashboard (replaces the unmaintained-feel `shinydashboard`). One tab per stage:

- **Setup** — experiment metadata; chemical names & units; **response-type toggle
  (continuous / binary)**; **download template** and upload data; editable data table.
  - Continuous template columns: `C1, C2, C3, Response`.
  - Binary template columns: `C1, C2, C3, N_exposed, N_affected`.
  - (Columns scale to the number of chemicals: 1, 2, or 3.)
- **Single / Binary / Ternary** — each runs **auto-fit & compare** and shows results +
  plots, with an **Advanced** panel to override start values and fix parameters. Mixture
  stages reuse fitted parameters from the previous stage as starting values.
- **Summary / Export** — collected parameters, comparison tables, and plots; export to
  Excel/PDF.

## 6. Outputs & plots

- **Parameter + comparison block** (the Excel B10:P20 region): for CA and IA, the
  Ref / S/A / DR / DL fitted parameters (μmax, slopes, EC50s, `a`, `b`), SS or −2·logL, χ²,
  df, p-value, 95% CIs, and the most-parsimonious flag.
- **Residuals** (observed vs predicted).
- **Plots:**
  - single dose–response curves (log-x), observed points + fitted curve;
  - observed-vs-predicted diagnostic;
  - **3-D response surface** (`plotly`);
  - **2-D isobole / contour** plots (cf. Jonker Fig. 1).

## 7. Validation & deliverables

- **`testthat` regression tests**: load the two Excel datasets from `inst/extdata/`, run the
  engine, assert it selects the published best model and lands within tolerance / within the
  reported confidence intervals.
- **README**: install + run instructions, quick start, screenshots.
- **Methodology vignette**: concise recap of Jonker et al. + worked example.
- **Licence**: GPL-3.

## 8. Dependencies (lean, maintained)
`shiny`, `bslib`, `bbmle`, `DEoptim`, `plotly`, `DT`, `readxl`/`writexl` (or `openxlsx`),
`numDeriv`. Drops from the old app: `drc`, `nloptr`, `Rsolnp`, `rhandsontable`, `Cairo`,
`shinydashboard`.

## 9. Out of scope (v1)
- Hosted/click-a-link web deployment (shinylive/webR) — possible later; not designed for now.
- >3 chemicals.
- Non-monotonic (J-shaped / hormesis) dose–response curves.

## 10. Open items
- The referenced email ("[Pasted text #1]") did not reach the assistant; revisit if it
  contains additional constraints.
- Confirm GPL-3 vs MIT licence with the authors.
- Confirm exact slope-sign convention and parameter bounds against the workbooks during
  implementation.
