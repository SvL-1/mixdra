# Design: mixdra methodology vignette (Plan 6)

**Date:** 2026-06-04
**Status:** Approved (brainstorming) — ready for implementation plan
**Roadmap:** Plan 6 (docs) of the mixdra plan split. README is done; this is the
remaining "methodology vignette" deliverable. Builds on the finished engine, plotting
layer, app, and I/O — no engine changes.

## Goal

Ship a single HTML package vignette that explains the Jonker et al. (2005) mixture
dose–response method *as mixdra implements it*, taught through one real worked example
end-to-end. It is the deep, narrative companion to the README (which remains the
"why + quick reference"). A reader who finishes it should understand both **what the
method does** and **how to drive the package to do it**.

## Shape

- **One vignette:** `vignettes/methodology.Rmd`, e.g. titled
  *"The Jonker method, worked end-to-end with mixdra."*
- **Worked-example-led, theory inline:** lead with a dataset and walk through it; each
  piece of theory (curve equation, CA vs IA, the deviation terms, staged fitting,
  LR testing) is introduced at the moment it is first needed, not front-loaded.
- **Binary spine + ternary coda:** a binary mixture carries the full method narrative
  (CA/IA references + S/A/DR/DL deviations + LR model selection); a short closing
  section extends to ternary + Advanced S/A.
- Consumes the **public API only** — no internals, no engine changes.

## Worked-example data

Both bundled datasets already exist in `inst/extdata/`:

- Binary spine: `binary_ca_cpf_imi_fbsa_continuous.csv` (145 rows; columns `C1`, `C2`,
  `Res` — two concentrations + a continuous reproduction response).
- Ternary coda: `ternary_ca_fbsa_cpf_imi_continuous.csv` (420 rows; `C1`, `C2`, `C3`,
  `Res`).

Loaded in the vignette via `system.file("extdata", ..., package = "mixdra")` so the
example is reproducible from an installed package. Exact chemical identities/labels for
the two binary columns are confirmed from the data + provenance during implementation
(F. candida soil-reproduction assay: CPF / MPs / IMI / FBSA).

## Content outline (the narrative)

1. **Setup** — load `mixdra`; load the bundled binary example via `system.file()`;
   show its structure (concentration columns + continuous `Res`). One paragraph framing:
   this is the Excel/VBA + Solver workflow, made reproducible.
2. **The dose–response curve** — the 3-parameter log-logistic model; fit each compound
   from its single-compound rows (`analyse_single()`); plot the two curves
   (`plot_dose_response()`). Theory inline: the curve equation, meaning of `max`, slope,
   EC50.
3. **Reference models** — CA vs IA: concentration (toxic-unit) scale vs probabilistic
   (effect) scale; what each assumes; fit both as the additivity baselines.
4. **Deviation models** — S/A (overall synergism/antagonism), DR (dose-ratio-dependent),
   DL (dose-level-dependent) layered on a reference. The interaction-term equations shown
   inline (Jonker 2005 forms). **Staged fitting** explained: curve parameters are fixed
   from the single-compound fits, then only the interaction parameter(s) are fit to the
   mixture rows.
5. **Model selection** — `analyse_mixture()` fits the reference → S/A → {DR, DL} chain;
   walk through `result$comparison` (LR test: χ², df, p vs parent) and how the most
   parsimonious adequate model becomes `result$chosen`.
6. **Interpreting & visualising** — `result_table()` for the parameter block, then the
   diagnostic plots: `plot_obs_pred()`, `plot_surface()`, `plot_isobole()`.
7. **Ternary coda** — extend to three chemicals: Advanced S/A (per-corner A1–A3 + the
   three-way A4); the **overall-vs-per-ratio "averaging-out"** insight (synergy and
   antagonism on different ratios can cancel in a single overall fit); `analyse_ternary()`
   + `ternary_effect_table()`; one `plot_isoplane()` (and/or `plot_sigma_tu()`) figure.
8. **References** — Jonker et al. (2005); van Loon et al. (2025). Cross-link the README.

## Build mechanics — the pre-compute pattern

The engine fits take minutes (multi-start), plots are interactive plotly, Rtools is not
present here, and CRAN/CI cap vignette build time. So the vignette is **pre-computed**:
the slow code runs once, locally, and the shipped vignette is static.

- **Source of truth:** `vignettes/methodology.Rmd.orig` — fully runnable, real code,
  real fits. Uses enough `n_starts` and a fixed `set.seed()` so DR/DL can actually be
  selected (slowness is acceptable; it runs once). Listed in `.Rbuildignore` so it is
  not shipped/built.
- **Precompute script:** `vignettes/precompute.R` knits `.Rmd.orig` →
  `vignettes/methodology.Rmd` (e.g. `knitr::knit("methodology.Rmd.orig",
  "methodology.Rmd")`), with computed results and figure paths baked in. `R CMD build`
  then re-runs nothing slow: it just assembles the already-knitted `.Rmd`. No Rtools,
  no plotly, no fitting on the build path. (Same pattern as `data.table` / `sf`.)
- **Plots → static PNG (decided).** During precompute, each plotly object is rendered to
  PNG once via the established headless path (`htmlwidgets::saveWidget()` →
  `webshot2::webshot()`), saved under `vignettes/` (e.g. `vignettes/figure-*.png`), and
  embedded with `knitr::include_graphics()`. Static PNG is chosen over embedded
  interactive widgets: far lighter, CRAN-safe, no JS dependency shipped in the package.
- **DESCRIPTION:** add `knitr` and `rmarkdown` to `Suggests`; add `VignetteBuilder: knitr`.
  `plotly`, `webshot2`, and `DT` are only needed at precompute time (already Suggests /
  installed locally), not at build time.

### Repository layout after this work

```
vignettes/
  methodology.Rmd.orig     # runnable source (in .Rbuildignore)
  methodology.Rmd          # knitted, committed, shipped
  precompute.R             # knit .orig -> .Rmd, render figures
  figure-*.png             # committed static figures
.Rbuildignore              # add: ^vignettes/methodology\.Rmd\.orig$, ^vignettes/precompute\.R$
DESCRIPTION                # + knitr, rmarkdown (Suggests); + VignetteBuilder: knitr
```

## Acceptance / testing

- **Acceptance gate:** `vignettes/precompute.R` runs to completion locally and produces a
  `methodology.Rmd` that knits to HTML without error, with all figures present. Because
  the source code is real and executed, a clean precompute proves the documented workflow
  matches actual package behaviour.
- **No new unit tests.** The vignette uses only the public API, which is already covered
  by the existing `testthat` suite (engine, plotting, ternary, recovery). Adding
  engine-level tests here would duplicate that coverage.
- The committed artefacts are `methodology.Rmd` + the `figure-*.png` files; `.Rmd.orig`
  and `precompute.R` are committed for reproducibility but excluded from the build.

## Out of scope

- Multiple/separate vignettes (getting-started, theory-only, per-scale) — one vignette only.
- A `pkgdown` site.
- Re-running fits at `R CMD build` time (explicitly avoided via pre-compute).
- Interactive plotly widgets embedded in the vignette (static PNG instead).
- Any engine, plotting, app, or I/O code changes.

## Environment notes (carried from prior work)

- Non-interactive R here needs `.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))`;
  RStudio finds packages automatically.
- `R CMD check` / `rcmdcheck` cannot fully run here (no Rtools); run the final
  `devtools::check()` / build on a machine with Rtools.
- Headless plotly→PNG: `saveWidget()` to HTML then `webshot2::webshot(html, png,
  vwidth=, vheight=, delay=2)`; `plotly`, `webshot2` (+ chromote) are in the user library.
- This terminal can drop/double streamed stdout on long runs — write markers to a file
  and Read them; background long knits and Monitor for a completion marker.
