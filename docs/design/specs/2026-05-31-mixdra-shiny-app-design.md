# mixdra Shiny App (single + binary) — Design

**Date:** 2026-05-31
**Status:** Approved design — ready for an implementation plan.
**Plan slot:** This is "Plan 4 — Shiny app" from the Plan-2 split, brought forward
ahead of the ternary subproject at the user's request.

Builds on the finished engine
(`docs/design/specs/2026-05-30-mixdra-design.md`,
`docs/design/plans/2026-05-30-mixdra-engine.md`) and the validated plotting
layer (`docs/design/specs/2026-05-30-mixdra-binary-plots-design.md`). The app
adds **no model math and no plotting math** — it is a thin, testable UI shell over
the package's exported API.

---

## 1. Goal & scope

Provide an interactive Shiny front-end for the proven part of the old
`MixTox_shiny_v2` workflow — **single-chemical** and **binary-mixture**
dose-response analysis — rebuilt on the `mixdra` engine and the `plotly` plotting
layer, packaged so it launches with `run_app()`.

**In scope (v1):** three stages — **Introduction**, **Single Chemical**, **Binary
Mixture** — modernised from `shinydashboard` to `bslib`, with `plotly` plots and
the `mixdra` fitting engine.

**Out of scope (own later plans):**
- Ternary stage (its engine subproject is on hold;
  `docs/design/specs/2026-05-30-ternary-subproject-design.md`).
- A Summary/aggregate-results tab.
- Rich file I/O — Excel templates, result-file export, multi-format upload — stays
  with **Plan 5 (file I/O)**. v1 includes only a minimal CSV-template download and
  fixed-format CSV upload (see §4).
- Asynchronous fitting (`ExtendedTask`/`promises`). v1 fits synchronously (§6).

### Decisions captured during brainstorming (2026-05-31)

| # | Decision |
|---|---|
| Scope | Faithful core: Intro + Single + Binary only. Ternary/Summary deferred. |
| Data input | **File upload only** — no in-app paste/editing of values. |
| File format | **CSV only** (`.csv`). |
| Column handling | **Fixed column schema** + a **downloadable template** the user fills in. No column-mapping UI, no auto-detect. |
| Fit controls | **Expose everything, default everything** — progressive disclosure: primary controls up front, an Advanced panel for the rest, all pre-filled with engine defaults. |
| Layout | `bslib` **top navbar + per-stage sidebar** (`page_navbar` + `layout_sidebar`). |
| Single response | Single tab supports continuous and quantal (quantal fit as proportions). |
| Binary plots | Default to the **chosen** model; a selector lets the user view any of the four fits; the isobole always overlays the reference. |

---

## 2. Architecture

A package-internal Shiny app: modules live in `R/` (so they are unit-testable with
`shiny::testServer`), the assembled app lives in `inst/app/`, and a single exported
`run_app()` launches it.

```
R/app-run.R       run_app(): guards on Suggests, launches inst/app. Exported + documented.
R/app-intro.R     intro_ui(id) / intro_server(id, meta)
R/app-single.R    single_ui(id) / single_server(id, meta)
R/app-binary.R    binary_ui(id) / binary_server(id, meta)
R/app-io.R        shared, plotly-free, Shiny-free helpers:
                    template_csv(stage, response)  -> character/data.frame template
                    validate_upload(df, stage, response) -> ok | structured error
                    read_upload(path)              -> data.frame
inst/app/app.R    bslib::page_navbar() wiring the three modules + a shared `meta` reactiveValues
tests/testthat/test-app-io.R     pure-helper unit tests (no Shiny)
tests/testthat/test-app-modules.R testServer tests for each module server
tests/testthat/test-app-smoke.R   shinytest2 app-builds smoke test (skip_if_not_installed)
DESCRIPTION       add shiny, bslib, plotly, DT to Suggests
```

- **UI shell:** `bslib::page_navbar` with three `nav_panel`s. Each analysis panel is
  a `bslib::layout_sidebar` — controls in the sidebar, `plotly`/`DT` outputs in the
  main area. New stages later (Ternary, Summary) are just more `nav_panel`s.
- **Shared state:** one `reactiveValues` object, `meta`, holding the Introduction
  tab's experiment ID, species, endpoint, chemical-1 name, chemical-2 name, and
  concentration unit. Single and Binary read it for plot titles and axis labels.
  **No fitted results flow between stages** (see §3).
- **Dependency isolation:** `shiny`, `bslib`, `plotly`, `DT` are all **Suggests**.
  `run_app()` checks for them and stops with an install hint if any are missing.
  The engine and plotting layer therefore remain installable and `R CMD
  check`-clean without the UI stack.

---

## 3. State & data flow

```
Introduction ──(meta: names, unit, endpoint, …)──▶ Single Chemical
                                              └────▶ Binary Mixture
```

- `meta` carries **labels only**. Stages do not seed each other.
- **Why no Single→Binary seeding (unlike the old app):** `analyse_mixture(df, …)`
  seeds its own starting values from the binary dataset's single-chemical rows
  (the rows where every other chemical is 0), via `seed_from_singles()`. The Binary
  stage is therefore self-contained: its uploaded CSV holds the full dataset
  (single-chemical series + mixture rows), exactly as the validated fixture does.
  The Single stage is for standalone per-chemical inspection.

---

## 4. Upload, templates, validation

Fixed schema per stage and response type. The template download writes a CSV with
the correct headers plus a couple of example rows; the user fills it in Excel and
uploads it back.

| Stage | Response | Required columns |
|---|---|---|
| Single | continuous | `Conc, Res` |
| Single | quantal | `Conc, Affected, Exposed` |
| Binary | continuous | `C1, C2, Res` |
| Binary | quantal | `C1, C2, Affected, Exposed` |

- The **Binary** template includes example single-chemical rows (`C2 = 0` series and
  `C1 = 0` series) plus mixture rows, with a header comment noting the
  single-chemical rows are required for seeding.
- The app maps stage column names to the engine's internal convention before
  fitting (Single's `Conc` → the engine's single-curve input; Binary's `C1/C2`
  are already the engine convention).
- `validate_upload()` checks, before the engine is ever called:
  - all required columns present (report missing **and** unexpected columns);
  - all values numeric;
  - concentrations ≥ 0;
  - quantal: `Affected` and `Exposed` non-negative integers with `Affected ≤
    Exposed`;
  - enough rows to fit (single: ≥ 4 distinct concentrations).
- Failures render as a clear inline message; the Fit button stays disabled until a
  valid file is loaded.

Richer templates (Excel, pre-populated metadata) and result-file export are **Plan
5**, not v1.

---

## 5. Stage specifications

### 5.1 Introduction
- Inputs: Experiment ID, Species, Endpoint, Chemical 1 name, Chemical 2 name,
  Concentration unit. (Chemical 3 omitted — ternary is out of scope.)
- Writes each value to `meta`. No fitting, no outputs beyond the form.

### 5.2 Single Chemical
- **Sidebar:** response-type selector (continuous / quantal); `Download template`;
  file upload; `Fit` (enabled only after a valid upload). No Advanced panel in v1 —
  `fit_single` minimises sum-of-squares with positivity bounds and takes no
  reference/deviation/constraint arguments, so there is nothing further to expose.
- **Engine call:** `analyse_single(df)` (wraps `fit_single`). The fit object is
  enriched (`kind = "single"`) so the plot builders accept it directly.
- **Main area:**
  - dose-response curve — `plot_dose_response(fit, df)`;
  - observed-vs-predicted — `plot_obs_pred(fit, df)`;
  - parameter table (`DT`): max, slope, EC50, n, SSR.
- Axis/title labels from `meta` (chemical name + unit + endpoint), falling back to
  generic labels when blank.

### 5.3 Binary Mixture
- **Sidebar (primary):** response-type (continuous / quantal); reference model
  (CA / IA); `Download template`; file upload; **Thorough fit** toggle; `Fit`.
- **Sidebar (Advanced accordion, collapsed, defaults pre-filled):**
  - `n_starts` (numeric; the Thorough-fit toggle sets 1 ↔ 20);
  - `alpha` (default 0.05);
  - `time_limit` seconds per model (default 30);
  - per-parameter **lower / upper / fix** for `max, slope1, slope2, ec501, ec502`
    (binary `max` upper pre-filled at 1; user may fix it, e.g. to 0.98). `a`/`b`
    deviation parameters stay unconstrained per the constraints design.
- **Engine call:**
  `analyse_mixture(df, reference, response, n_starts, alpha, lower, upper, time_limit)`.
- **Main area:**
  - **Results table** (`result_table(res)`) — parameters for reference/SA/DR/DL +
    objective + df; the chosen model highlighted.
  - **Model comparison** (`res$comparison`) — LR χ², df, p vs each parent.
  - **Confidence intervals** (`param_ci`) for the selected fit.
  - **Plots**, driven by a "model to display" selector (defaults to `res$chosen`):
    dose-response per chemical (`plot_dose_response`, chem 1 / chem 2),
    observed-vs-predicted (`plot_obs_pred`), **3-D surface** (`plot_surface`), and
    **2-D isobole** (`plot_isobole`, `reference_fit = res$fits$reference`).

---

## 6. Long-running fits

- Default `n_starts = 1` keeps the normal `Fit` responsive (seconds).
- **Thorough fit** sets multi-start (`n_starts = 20`), explicitly opt-in, with an
  inline "this may take several minutes" warning.
- Each fit runs inside `shiny::withProgress` with a spinner + notification so the UI
  never appears frozen. `time_limit` bounds each model's runtime.
- Fitting is **synchronous** in v1. Asynchronous execution
  (`ExtendedTask`/`promises`/`future`) is a documented future enhancement, not built
  now (the app targets local single-user use).

---

## 7. Error handling

- **Upload/validation:** friendly inline messages (missing/unexpected columns,
  non-numeric, out-of-range, `Affected > Exposed`, too few points). Fit disabled
  until valid.
- **Fit failures:** wrapped in `tryCatch`; non-convergence or engine errors surface
  as a `showNotification` reporting the convergence status, without crashing the
  session. Plots/tables that depend on a fit are guarded with `req()`.
- **Missing dependencies:** `run_app()` stops early with an `install.packages(...)`
  hint listing the missing Suggests.

---

## 8. Testing strategy

- **Pure helpers** (`template_csv`, `validate_upload`, `read_upload`) — plain
  `testthat` unit tests; no Shiny runtime. These carry the bulk of the coverage.
- **Module servers** — `shiny::testServer` drives inputs and asserts reactive
  outputs and that the correct engine/plot functions are invoked. Headless and fast.
- **Smoke** — one `shinytest2` test that the app object builds / launches,
  `skip_if_not_installed("shinytest2")` and friends.
- The existing engine + plotting suite (54 tests) must stay green; all UI deps are
  `Suggests` + guarded, so `R CMD check` stays clean.

---

## 9. Dependencies

Added to `DESCRIPTION` **Suggests** (guarded, never `Imports`): `shiny`, `bslib`,
`plotly` (already present), `DT`, and for tests `shinytest2`. No new hard
dependencies; the engine remains lean.

---

## 10. Open items for the implementation plan

- Exact `bslib` theme / version pin and minimum Shiny version.
- Final template example rows (use trimmed real values from the binary fixture so
  the example is scientifically sensible).
