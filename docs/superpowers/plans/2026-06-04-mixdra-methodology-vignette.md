# mixdra Methodology Vignette Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship one pre-computed HTML package vignette, `vignettes/methodology.Rmd`, that teaches the Jonker et al. (2005) mixture method as mixdra implements it, through one worked binary example end-to-end plus a short ternary coda.

**Architecture:** The vignette is *pre-computed* (the `data.table`/`sf` pattern). A runnable `vignettes/methodology.Rmd.orig` holds the real code and real (slow) fits. `vignettes/precompute.R` runs `knitr::knit()` on it ONCE to produce the committed, static `vignettes/methodology.Rmd` (executable `{r}` chunks become inert markdown, so `R CMD build` re-runs nothing) plus committed static PNG figures. Plots (plotly) are rendered to PNG during precompute via `htmlwidgets::saveWidget()` → `webshot2::webshot()`. Content tasks author one section at a time into the `.Rmd.orig` and verify that section's R code runs standalone (fast); the full precompute is a single final assembly task.

**Tech Stack:** R, knitr, rmarkdown, html_vignette; mixdra public API (`analyse_single`, `analyse_mixture`, `result_table`, `plot_*`, `analyse_ternary`, `ternary_effect_table`); plotly + webshot2 + htmlwidgets (precompute-time only).

**Spec:** `docs/superpowers/specs/2026-06-04-mixdra-methodology-vignette-design.md`

---

## Environment preamble (read once)

Every non-interactive R invocation in this plan MUST set the user library first. Use this exact prefix inside R scripts:

```r
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
library(devtools)
```

For ad-hoc verification, prefer writing the package state with `devtools::load_all(".")` (no Rtools/compile needed — mixdra is pure R) over installing. The bundled example data is read with `system.file("extdata", "<file>.csv", package = "mixdra")` only AFTER `load_all()` (which makes `inst/extdata` visible). During development you may also read it directly: `read.csv("inst/extdata/binary_ca_cpf_imi_fbsa_continuous.csv")`.

**Terminal caveat:** this terminal drops/doubles streamed stdout on long runs. For any step that runs fits, have the R script `writeLines()` its key results AND a final marker line to a file, then Read that file — do not rely on streamed console output. Long runs auto-background; Monitor for the marker line.

---

## File structure

- `vignettes/methodology.Rmd.orig` — **Create.** Runnable source. Holds YAML, a hidden setup chunk (knitr opts + the `save_png()` figure helper), all prose, and all executable `{r}` chunks. In `.Rbuildignore`.
- `vignettes/precompute.R` — **Create.** Knits `.Rmd.orig` → `methodology.Rmd`, then renders to HTML to verify, writing a DONE marker. In `.Rbuildignore`.
- `vignettes/methodology.Rmd` — **Generated & committed.** The static, shipped vignette (output of precompute). Do NOT hand-edit.
- `vignettes/figure-*.png` — **Generated & committed.** Static figures referenced by the vignette. Ship in the build.
- `DESCRIPTION` — **Modify.** Add `knitr`, `rmarkdown` to `Suggests`; add `VignetteBuilder: knitr`.
- `.Rbuildignore` — **Modify.** Exclude `methodology.Rmd.orig` and `precompute.R` (but NOT the `.Rmd` or PNGs).

Authoring order: the `.Rmd.orig` is built up section by section (Tasks 2–8). Each section is a self-contained block of prose + chunks appended to the file. The full knit happens only in Task 9.

---

### Task 1: Vignette infrastructure + precompute scaffold

**Files:**
- Modify: `DESCRIPTION`
- Modify: `.Rbuildignore`
- Create: `vignettes/methodology.Rmd.orig`
- Create: `vignettes/precompute.R`

- [ ] **Step 1: Add vignette build deps to DESCRIPTION**

In `DESCRIPTION`, add `knitr` and `rmarkdown` to the `Suggests:` block (alphabetical order) and add a `VignetteBuilder` field. The `Suggests:` block becomes:

```
Suggests:
    bslib,
    DT,
    knitr,
    plotly,
    readxl,
    rmarkdown,
    shiny,
    shinytest2,
    testthat (>= 3.0.0)
VignetteBuilder: knitr
```

(Place `VignetteBuilder: knitr` immediately after the `Suggests:` block, before `Config/testthat/edition`.)

- [ ] **Step 2: Exclude the source + driver from the build**

Append to `.Rbuildignore`:

```
^vignettes/methodology\.Rmd\.orig$
^vignettes/precompute\.R$
```

- [ ] **Step 3: Create the skeleton `vignettes/methodology.Rmd.orig`**

Create the file with the YAML, hidden setup chunk (knitr options + `save_png()` helper), and a single trivial chunk so the scaffold knits before any heavy content exists:

````markdown
---
title: "The Jonker method, worked end-to-end with mixdra"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The Jonker method, worked end-to-end with mixdra}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figure-",
  echo = TRUE
)

# Render a plotly object to a static PNG under vignettes/ and return the file
# name, so chunks can do: knitr::include_graphics(save_png(p, "name.png")).
# Static PNGs keep the shipped vignette light and free of a JS dependency.
save_png <- function(widget, file, vwidth = 700, vheight = 500) {
  html <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(widget, html, selfcontained = TRUE)
  webshot2::webshot(html, file, vwidth = vwidth, vheight = vheight, delay = 2)
  file
}

library(mixdra)
```

# Introduction

_Scaffold — replaced in later tasks._

```{r scaffold}
1 + 1
```
````

- [ ] **Step 4: Create `vignettes/precompute.R`**

```r
# Pre-compute the methodology vignette: knit the runnable .Rmd.orig into the
# static, shipped methodology.Rmd (executable chunks become inert markdown, so
# R CMD build re-runs nothing), then render to HTML to confirm it builds.
# Run from the package root:  Rscript vignettes/precompute.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))

old <- setwd("vignettes")
on.exit(setwd(old), add = TRUE)

knitr::knit("methodology.Rmd.orig", output = "methodology.Rmd")
rmarkdown::render("methodology.Rmd", quiet = TRUE)  # writes methodology.html

writeLines("PRECOMPUTE_DONE", "precompute.marker")
```

- [ ] **Step 5: Run the scaffold precompute and verify it builds**

Run: `Rscript vignettes/precompute.R`
Then Read `vignettes/precompute.marker` and confirm it contains `PRECOMPUTE_DONE`, and confirm `vignettes/methodology.Rmd` and `vignettes/methodology.html` now exist (`ls vignettes`).
Expected: marker present; both output files exist; no error from knit/render.

- [ ] **Step 6: Commit**

```bash
git add DESCRIPTION .Rbuildignore vignettes/methodology.Rmd.orig vignettes/precompute.R vignettes/methodology.Rmd vignettes/methodology.html
git rm --cached vignettes/precompute.marker 2>/dev/null; echo "vignettes/precompute.marker" >> .gitignore
git add .gitignore
git commit -m "docs(vignette): scaffold pre-computed methodology vignette"
```

---

### Task 2: Setup & data-loading section

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Replace the scaffold Introduction with the real intro + data load**

Replace everything from `# Introduction` to the end of the file with:

````markdown
# Introduction

The interaction-analysis method of Jonker et al. (2005) is usually run by hand in
an Excel/VBA workbook with the Solver add-in: one model fitted per click, no test
suite, no record of the starting values. `mixdra` reproduces the same science as a
tested R package. This vignette walks the whole method through one real dataset, so
you see both *what each step does* and *how to drive the package to do it*.

We use a bundled binary mixture: a *Folsomia candida* soil reproduction assay with
two co-applied chemicals. Columns `C1` and `C2` are the two concentrations (mg/kg);
`Res` is the continuous reproduction response (offspring count).

```{r load-data}
path <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv",
                    package = "mixdra")
df <- read.csv(path)
str(df)
head(df)
```

The design crosses the two chemicals: rows with `C2 == 0` are the chemical-1-only
series, rows with `C1 == 0` are the chemical-2-only series, the `C1 == 0 & C2 == 0`
rows are untreated controls, and the remainder are genuine mixtures.
````

- [ ] **Step 2: Verify the data-load code runs standalone**

Run this verification script (write results to a file to dodge the streamed-stdout caveat):

```r
# file: /tmp/verify-t2.R  (or any scratch path)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
path <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra")
df <- read.csv(path)
writeLines(c(
  paste("nrow:", nrow(df)),
  paste("cols:", paste(names(df), collapse = ",")),
  "T2_DONE"
), "verify-t2.out")
```

Run: `Rscript /tmp/verify-t2.R` then Read `verify-t2.out`.
Expected: `nrow: 145`, `cols: C1,C2,Res`, `T2_DONE`. Delete `verify-t2.out` after.

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): intro + data-loading section"
```

---

### Task 3: Single-compound dose-response curves section

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Append the single-curve section**

Append to `vignettes/methodology.Rmd.orig`:

````markdown
# Step 1 — the single-compound curves

Each chemical alone is described by a three-parameter log-logistic curve

$$ y(c) \;=\; \frac{\mathrm{max}}{1 + (c / \mathrm{EC50})^{\,\mathrm{slope}}}, $$

where `max` is the control-level response, `EC50` is the concentration giving half
that response, and `slope` controls steepness. `analyse_single()` fits this curve
from a one-chemical data frame (`C1` + `Res`). We fit each chemical from its own
single-compound rows.

```{r single-fits}
s1 <- subset(df, C2 == 0)                              # chemical-1-only + controls
s2 <- subset(df, C1 == 0)                              # chemical-2-only + controls

fit1 <- analyse_single(s1)
fit2 <- analyse_single(data.frame(C1 = s2$C2, Res = s2$Res))

rbind(`chemical 1` = fit1$par, `chemical 2` = fit2$par)
```

These per-chemical `max`, `slope`, and `EC50` values are exactly what the staged
mixture fit (Step 3) will hold fixed: the curves come from the single-compound data,
never from the mixture rows.
````

(`analyse_single()` returns a `fit_single()` result whose `$par` is the named
`max`/`slope`/`ec50` vector. If the field name differs, adjust the readout in
Step 2 — verify before committing.)

- [ ] **Step 2: Verify the single-fit code runs and inspect the return shape**

```r
# file: /tmp/verify-t3.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
df <- read.csv(system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra"))
s1 <- subset(df, C2 == 0); s2 <- subset(df, C1 == 0)
fit1 <- analyse_single(s1)
fit2 <- analyse_single(data.frame(C1 = s2$C2, Res = s2$Res))
sink("verify-t3.out")
cat("names(fit1):", paste(names(fit1), collapse = ","), "\n")
cat("fit1$par:\n"); print(fit1$par)
cat("fit2$par:\n"); print(fit2$par)
cat("T3_DONE\n")
sink()
```

Run: `Rscript /tmp/verify-t3.R` then Read `verify-t3.out`.
Expected: both fits return a named parameter vector (`max`, `slope`, `ec50`); `T3_DONE` present. If `$par` is not the right accessor, note the correct one and fix the chunk's readout line. Delete `verify-t3.out` after.

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): single-compound curves section"
```

---

### Task 4: Reference models (CA vs IA) section

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Append the reference-models section**

Append to `vignettes/methodology.Rmd.orig`:

````markdown
# Step 2 — the additivity reference models

Before asking whether two chemicals *interact*, we need a baseline for "no
interaction". The method offers two:

- **Concentration Addition (CA)** assumes the chemicals act through the same
  mechanism and add up on the *concentration* (toxic-unit) scale: a mixture is as
  toxic as the sum of each chemical's potency-scaled dose. CA predicts the mixture
  response by solving, at each mixture point, for the equi-effective single dose.
- **Independent Action (IA)** assumes independent mechanisms and combines on the
  *effect* (probability) scale: the chemicals act on different sites, so their
  fractional effects multiply.

Neither has extra free parameters beyond the single-compound curves — they are pure
predictions. `mixdra` fits both as the additivity baselines; departures from them
are what the deviation models (Step 3) capture. This worked example uses CA as the
reference throughout (it is the reference behind the source workbook); IA is fit the
same way by passing `reference = "IA"`.
````

- [ ] **Step 2: No code to verify (prose-only task)**

This section introduces no new code chunk; the references are fit by
`analyse_mixture()` in Task 5. Confirm the file still reads coherently (the next
task wires CA in).

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): CA vs IA reference-models section"
```

---

### Task 5: Deviation models, staged fitting, and the analyse_mixture run

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Append the deviation-models + fit section**

Append to `vignettes/methodology.Rmd.orig`:

````markdown
# Step 3 — deviation models and staged fitting

A *deviation* model adds a single interaction term on top of a reference, each
capturing a different *pattern* of non-additivity (Jonker et al. 2005, Table 1):

- **S/A** — one parameter `a`: overall synergism (`a > 0`) or antagonism (`a < 0`),
  the same everywhere in the mixture.
- **DR** — dose-ratio-dependent: the interaction changes with the *ratio* of the two
  chemicals (parameters `a` and `b`).
- **DL** — dose-level-dependent: the interaction changes with the overall *dose
  level* / effect magnitude (parameters `a` and `b`).

`mixdra` fits these in **stages**, which is the whole point of the method: the curve
parameters are estimated from the single-compound rows and then *held fixed*, and
only the interaction parameter(s) are fitted to the full mixture data. The mixture
points can therefore never distort the single-chemical potencies.

`analyse_mixture()` does the entire chain in one call — it fits the reference and
each deviation, then compares them (Step 4). We raise `n_starts` so the multi-start
optimiser reliably escapes the local minima of the non-smooth CA surface.

```{r fit-mixture}
set.seed(1)
res <- analyse_mixture(df, reference = "CA", response = "continuous",
                       n_starts = 12)
names(res)
names(res$fits)   # reference, SA, DR, DL
```
````

- [ ] **Step 2: Verify the full mixture fit runs end-to-end (SLOW — background it)**

This step runs the real multi-start fit (minutes). Run it in the background and Monitor for the marker.

```r
# file: /tmp/verify-t5.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
df <- read.csv(system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra"))
set.seed(1)
res <- analyse_mixture(df, reference = "CA", response = "continuous", n_starts = 12)
sink("verify-t5.out")
cat("names(res):", paste(names(res), collapse = ","), "\n")
cat("fits:", paste(names(res$fits), collapse = ","), "\n")
cat("chosen:", res$chosen, "\n")
cat("T5_DONE\n")
sink()
saveRDS(res, "verify-t5.rds")   # reuse in later verification to avoid refitting
```

Run (background): `Rscript /tmp/verify-t5.R`
Monitor `verify-t5.out` for `T5_DONE`, then Read it.
Expected: `names(res)` includes `fits,comparison,chosen,reference,response`; `fits` are `reference,SA,DR,DL`; `chosen` is one of those names. Keep `verify-t5.rds` for Tasks 6–7. Delete `verify-t5.out`.

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): deviation models + staged fitting section"
```

---

### Task 6: Model-selection (likelihood-ratio) section

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Append the model-selection section**

Append to `vignettes/methodology.Rmd.orig`:

````markdown
# Step 4 — which interaction does the data support?

Each deviation model nests its reference (set the interaction parameters to zero and
you recover CA), so nested pairs can be compared by a **likelihood-ratio test**.
`analyse_mixture()` builds the comparison table and applies the parsimony rule:
prefer the simplest model the data do not reject.

```{r comparison}
res$comparison
res$chosen
```

Read the table top-down: each row tests a model against its parent (the change in
objective against the degrees of freedom it spent). A small *p*-value means the
extra interaction parameter bought a real improvement; a large one means it did not,
and the simpler parent stands. `res$chosen` is the model this rule selects — the
interaction pattern the data actually justify, not the most flattering one.
````

- [ ] **Step 2: Verify the comparison readout from the saved fit**

```r
# file: /tmp/verify-t6.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
res <- readRDS("verify-t5.rds")
sink("verify-t6.out")
cat("comparison columns:", paste(names(res$comparison), collapse = ","), "\n")
print(res$comparison)
cat("chosen:", res$chosen, "\n")
cat("T6_DONE\n")
sink()
```

Run: `Rscript /tmp/verify-t6.R` then Read `verify-t6.out`.
Expected: `res$comparison` is a data frame with one row per non-root model (LR statistic, df, p-value vs parent); `T6_DONE` present. Confirm the prose's description of the columns matches the actual column names; adjust the prose wording if needed. Delete `verify-t6.out`.

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): LR-test model-selection section"
```

---

### Task 7: Interpreting & visualising section (binary plots)

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Append the interpret + plots section**

Append to `vignettes/methodology.Rmd.orig`. Each plot chunk renders to a committed PNG via `save_png()` and embeds it with `include_graphics()`; the plotly object is built from the *chosen* fit (`res$fits[[res$chosen]]`) with the reference fit overlaid on the isobole.

````markdown
# Step 5 — reading and visualising the result

`result_table()` turns the selected fit into the parameter block you would read off
the paper or workbook: fitted `max`, per-chemical slopes and EC50s, the interaction
parameter(s), and the fit statistic.

```{r result-table}
result_table(res)
```

The fitted marginal curve for each chemical, with the observed points overlaid:

```{r plot-curve-1, results = "hide"}
fit <- res$fits[[res$chosen]]
knitr::include_graphics(save_png(plot_dose_response(fit, df, chem = 1),
                                 "figure-curve1.png"))
```

```{r plot-curve-2, results = "hide"}
knitr::include_graphics(save_png(plot_dose_response(fit, df, chem = 2),
                                 "figure-curve2.png"))
```

Observed-versus-predicted is the quickest global check of fit — points should hug
the dashed 1:1 line:

```{r plot-obs-pred, results = "hide"}
knitr::include_graphics(save_png(plot_obs_pred(fit, df), "figure-obs-pred.png"))
```

The full fitted response surface over the (C1, C2) plane, with the observed cloud:

```{r plot-surface, results = "hide"}
knitr::include_graphics(save_png(plot_surface(fit, df), "figure-surface.png"))
```

Finally the isoboles — equal-response contours. Overlaying the CA reference
(dashed) on the chosen model (solid) makes any departure from additivity visible: if
the solid contours bow away from the reference, the mixture is more (or less) potent
than additivity predicts.

```{r plot-isobole, results = "hide"}
knitr::include_graphics(save_png(
  plot_isobole(fit, df, reference_fit = res$fits$reference),
  "figure-isobole.png"))
```
````

- [ ] **Step 2: Verify each plot object builds and renders to PNG**

```r
# file: /tmp/verify-t7.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
df <- read.csv(system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra"))
res <- readRDS("verify-t5.rds")
fit <- res$fits[[res$chosen]]
save_png <- function(widget, file, vwidth = 700, vheight = 500) {
  html <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(widget, html, selfcontained = TRUE)
  webshot2::webshot(html, file, vwidth = vwidth, vheight = vheight, delay = 2)
  file
}
print(result_table(res))
save_png(plot_dose_response(fit, df, chem = 1), "t7-curve1.png")
save_png(plot_dose_response(fit, df, chem = 2), "t7-curve2.png")
save_png(plot_obs_pred(fit, df), "t7-obspred.png")
save_png(plot_surface(fit, df), "t7-surface.png")
save_png(plot_isobole(fit, df, reference_fit = res$fits$reference), "t7-isobole.png")
writeLines("T7_DONE", "verify-t7.out")
```

Run: `Rscript /tmp/verify-t7.R` then Read `verify-t7.out` and confirm the five `t7-*.png` files exist and are non-empty (`ls -l t7-*.png`).
Expected: `result_table(res)` prints a parameter block; all five PNGs created; `T7_DONE` present. Delete the `t7-*.png` scratch files and `verify-t7.out` after (the real figures are generated in Task 9).

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): interpret + binary plots section"
```

---

### Task 8: Ternary coda section

**Files:**
- Modify: `vignettes/methodology.Rmd.orig`

- [ ] **Step 1: Append the ternary coda**

Append to `vignettes/methodology.Rmd.orig`:

````markdown
# Coda — three chemicals and Advanced S/A

The same machinery scales to three chemicals. With three components there are more
ways to interact, so the S/A term is generalised to **Advanced S/A**: a per-corner
interaction with pairwise parameters `A1`, `A2`, `A3` (one per chemical pair) plus a
three-way `A4` (all three present together).

```{r load-ternary}
tpath <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                     package = "mixdra")
tdf <- read.csv(tpath)
str(tdf)
```

`analyse_ternary()` fits the same staged way: curves from the singles, pairwise terms
from the binary rows, and the three-way `A4` from the ternary rows. Crucially it fits
`A4` two ways — an **overall** `A4` across all ternary rows, and an **individual**
`A4` per mixture ratio.

```{r fit-ternary}
set.seed(1)
tern <- analyse_ternary(tdf, reference = "CA", n_starts = 12)
```

`ternary_effect_table()` contrasts them. Comparing each ratio's own `A4` against the
overall `A4` reveals *averaging-out*: synergy on one ratio and antagonism on another
can cancel in a single overall fit, so the per-ratio readout shows interaction the
pooled fit hides.

```{r ternary-effect-table}
ternary_effect_table(tern, tdf)
```

The EC50 **isoplane** is the three-chemical analogue of the isobole — the surface of
equi-effective mixtures across the concentration simplex:

```{r plot-isoplane, results = "hide"}
knitr::include_graphics(save_png(plot_isoplane(tern, tdf), "figure-isoplane.png"))
```

For the full theory — the reference-model equations, the deviation terms, and the
significance testing — see Jonker et al. (2005); for the microplastic/PFAS/pesticide
study this method was built for, see van Loon et al. (2025).

# References

- Jonker, M.J., Svendsen, C., Bedaux, J.J.M., Bongers, M. & Kammenga, J.E. (2005).
  Significance testing of synergistic/antagonistic, dose level–dependent, or dose
  ratio–dependent effects in mixture dose–response analysis. *Environmental
  Toxicology and Chemistry*, 24(10), 2701–2713.
  <https://doi.org/10.1897/04-431R.1>
- van Loon, S., Xie, G., Svendsen, C., Kraak, M.H.S., de Jeu, L., Schut, N.C.,
  Sprokkereef, E., Hurley, R., van Wezel, A.P. & van Gestel, C.A.M. (2025).
  Microplastics and PFAS as ubiquitous pollutants affect potencies of highly toxic
  chemicals in mixtures. *Journal of Hazardous Materials*, 500, 140493.
  <https://doi.org/10.1016/j.jhazmat.2025.140493>
````

- [ ] **Step 2: Verify the ternary fit + table + isoplane (SLOW — background it)**

```r
# file: /tmp/verify-t8.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
tdf <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra"))
set.seed(1)
tern <- analyse_ternary(tdf, reference = "CA", n_starts = 12)
et <- ternary_effect_table(tern, tdf)
save_png <- function(widget, file, vwidth = 700, vheight = 500) {
  html <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(widget, html, selfcontained = TRUE)
  webshot2::webshot(html, file, vwidth = vwidth, vheight = vheight, delay = 2)
  file
}
save_png(plot_isoplane(tern, tdf), "t8-isoplane.png")
sink("verify-t8.out")
cat("effect-table columns:", paste(names(et), collapse = ","), "\n")
print(et)
cat("T8_DONE\n")
sink()
```

Run (background): `Rscript /tmp/verify-t8.R`
Monitor `verify-t8.out` for `T8_DONE`, then Read it and confirm `t8-isoplane.png` exists.
Expected: `ternary_effect_table()` returns rows per ratio (`ratio,C1,C2,C3,pred_CA,pred_SA,pred_ASA,a4_effect`); isoplane PNG created; `T8_DONE`. Delete `t8-isoplane.png` and `verify-t8.out` after.

- [ ] **Step 3: Commit**

```bash
git add vignettes/methodology.Rmd.orig
git commit -m "docs(vignette): ternary coda + references section"
```

---

### Task 9: Final precompute, HTML verification, and ship

**Files:**
- Generated: `vignettes/methodology.Rmd`, `vignettes/figure-*.png`

- [ ] **Step 1: Run the full precompute (SLOW — background it)**

This runs every chunk for real (both 12-start fits + all figures) and produces the shipped `methodology.Rmd` + PNGs, then renders HTML to confirm it builds. Add a `set.seed` is already inside the chunks. Run in the background and Monitor.

Run (background): `Rscript vignettes/precompute.R`
Monitor `vignettes/precompute.marker` for `PRECOMPUTE_DONE`.
Expected: marker present; `vignettes/methodology.Rmd`, `vignettes/methodology.html`, and `vignettes/figure-curve1.png`, `figure-curve2.png`, `figure-obs-pred.png`, `figure-surface.png`, `figure-isobole.png`, `figure-isoplane.png` all exist.

- [ ] **Step 2: Sanity-check the generated `methodology.Rmd` is static**

Confirm the knit turned executable chunks inert. Grep the generated file:

Run: `grep -c '```{r' vignettes/methodology.Rmd` → Expected: `0` (no live `{r}` chunks remain; knit replaced them with plain ` ```r ` blocks). `grep -c 'figure-' vignettes/methodology.Rmd` → Expected: `6` (the six embedded figures).

If live `{r}` chunks remain, the build would re-run fits — stop and fix the precompute (it must use `knitr::knit`, not a copy).

- [ ] **Step 3: Verify the vignette is registered and the package builds the vignette index**

```r
# file: /tmp/verify-t9.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
sink("verify-t9.out")
cat("html exists:", file.exists("vignettes/methodology.html"), "\n")
cat("rmd static (no live chunks):",
    !any(grepl("```\\{r", readLines("vignettes/methodology.Rmd"))), "\n")
cat("T9_DONE\n")
sink()
```

Run: `Rscript /tmp/verify-t9.R` then Read `verify-t9.out`.
Expected: `html exists: TRUE`, `rmd static (no live chunks): TRUE`, `T9_DONE`. Delete `verify-t9.out`.
(Full `devtools::check()` / `R CMD build` to validate vignette packaging cannot run here — no Rtools; note this and run it on an Rtools machine before release.)

- [ ] **Step 4: Open the rendered HTML and eyeball it**

Open `vignettes/methodology.html` in a browser (or read it) and confirm: all six figures display, the parameter/comparison/effect tables rendered, equations render, and the narrative reads top-to-bottom without scaffold leftovers. Fix any prose/figure issues in `.Rmd.orig` and re-run Step 1 if needed.

- [ ] **Step 5: Commit the generated artefacts**

```bash
git add vignettes/methodology.Rmd vignettes/figure-curve1.png vignettes/figure-curve2.png vignettes/figure-obs-pred.png vignettes/figure-surface.png vignettes/figure-isobole.png vignettes/figure-isoplane.png
git commit -m "docs(vignette): pre-computed methodology vignette + figures"
```

(`vignettes/methodology.html` is a build/verification artefact — `R CMD build`
regenerates it into `inst/doc/`. It is git-ignored, NOT committed.)

- [ ] **Step 6: Update the README to point at the vignette**

In `README.md`, under the "The method, briefly" section, add a closing line:

```markdown
For a full worked walkthrough — fitting curves, comparing models, and reading the
results on a real dataset — see the **methodology vignette**
(`vignette("methodology", package = "mixdra")`).
```

Commit:

```bash
git add README.md
git commit -m "docs(readme): link to the methodology vignette"
```

---

## Self-review notes (for the implementer)

- **Spec coverage:** §Shape → Tasks 2–8 (binary spine + ternary coda, theory inline); §Content outline items 1–8 → Tasks 2,3,4,5,6,7,8 respectively; §Build mechanics (pre-compute, `.Rmd.orig`→`.Rmd`, static PNG, DESCRIPTION/.Rbuildignore) → Tasks 1 & 9; §Acceptance (precompute runs + knits to HTML, no new unit tests) → Task 9.
- **Unverified API accessors to confirm during implementation (flagged in-task):** `analyse_single()$par` field name (Task 3 Step 2); `res$comparison` column names for the prose (Task 6 Step 2). Both have an explicit verify-and-adjust step.
- **Plot-object/`fit` contract:** all binary `plot_*` take a single enriched fit (`res$fits[[res$chosen]]` or `res$fits$reference`), confirmed against `R/plot.R`; ternary plots take the `analyse_ternary()` result, confirmed against `R/plot.R` and `R/ternary-isoplane.R`.
- **Build-time safety:** Task 9 Step 2 explicitly asserts no live `{r}` chunks survive in the shipped `.Rmd`, which is the whole point of the pre-compute pattern.
