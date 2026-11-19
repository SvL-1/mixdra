# (De)potentiation Tab Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new top-level "(De)potentiation" tab that fits stressor 1's curve per series of stressor 2, tests each series' EC50 shift against the reference by likelihood ratio (Sam's workbook method), and states whether stressor 2 makes stressor 1 more or less toxic; the "Multiple stressors" tab is renamed "Mixture toxicity".

**Architecture:** A pure engine function `analyse_potentiation()` (new `R/potentiation.R`) does all fitting and testing with base-R `optim`; a pure summary builder turns its comparisons into the conclusion text; a plot data builder + plotly renderer draw the curves; a new Shiny module `R/app-potentiation.R` wires upload → engine → plot/tables, following the existing `*_ui()/*_server()` pattern. The workbook comparison enters the table-driven validation report.

**Tech Stack:** R package `mixdra`, base R `stats::optim` (L-BFGS-B), testthat 3, shiny + bslib + plotly + DT (Suggested), roxygen2.

**Spec:** `docs/design/specs/2026-11-19-depotentiation-tab-design.md` — read it first; this plan argues from it.

## Global Constraints

- Generic naming everywhere: no compound names (CPF, MPs, FBSA, IMI) in code, UI, example data, tests or oracle labels for this feature. Defaults are "Stressor 1" / "Stressor 2"; example series are `Stressor 1`, `Level 1`, `Level 2`, `Level 3`.
- Model `Y = max / (1 + (C/EC50)^b)` via existing `ll3_predict()`. Bounds: max in [1e-8, Inf] (quantal: upper 1), slope in [1e-3, 50], ec50 in [1e-8, Inf].
- Continuous: objective = residual SS; `X² = n · ln(obj2/obj1)`. Quantal: objective = `obj_deviance()` with predictions clamped to [1e-8, 1 − 1e-8]; `X² = obj2 − obj1`. df = 1. `p = pchisq(X², 1, lower.tail = FALSE)`. X² floored at 0. α default 0.05.
- Controls: a series with no `C1 == 0` rows borrows the reference series' `C1 == 0` rows (counted in both curves, as in the workbook).
- Base R + existing dependencies only. Internal `campaign_*` names stay; the rename touches user-facing strings only.
- Every non-interactive R call: prefix `.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))`. Run R through the **PowerShell** tool (`Rscript -e "..."`); the Bash tool's Rscript output is unreliable here. If output looks truncated, have R `writeLines()` to a file in the scratchpad and Read it.
- Tests: run **filtered** only — `devtools::test(filter = "<name>")`, never the full suite.
- Git: commit on branch `depotentiation-tab`. **Never** run HEAD-moving git (checkout, switch, reset, rebase, stash) — the working directory is shared with other agents. No `Co-Authored-By` trailer in commit messages.
- After adding/altering roxygen tags: `devtools::document()` and commit `NAMESPACE`/`man/` changes in the same task.

## Review Focus

1. **Series column read as factor/number** — `read.csv` gives character; labels like `1`, `2` would read as integers. Expect: `Series` is coerced with `as.character()` before use, and numeric-looking labels work. (Test in Task 3.)
2. **Reference series with no controls** — if the reference has no `C1 == 0` rows and another series lacks them too, nothing can be borrowed. Expect: the fit proceeds without controls (no error), `controls_borrowed = FALSE`. (Test in Task 1.)
3. **A series with all responses equal or a flat curve** — EC50 runs to a bound. Expect: no error; the fit returns, the comparison is reported (likely non-significant), and the EC50 is flagged extrapolated. (Test in Task 1.)
4. **Reference choice not first in file** — the user picks another series as reference. Expect: comparisons are all other series in file order; the summary's onset logic uses that order. (Test in Task 1.)
5. **Whitespace/case in Series labels** (`"Level 1 "` vs `"Level 1"`) — Expect: labels are trimmed with `trimws()` so stray spaces don't split a series. (Test in Task 3.)

---

## File structure

| File | Responsibility |
|---|---|
| `R/potentiation.R` (create) | `analyse_potentiation()` + internal fitting helpers |
| `R/potentiation-summary.R` (create) | conclusion text: `potentiation_summary()`, number formatting |
| `R/plot-data.R`, `R/plot.R` (modify) | `potentiation_curve_data()`, `plot_potentiation()` |
| `R/app-io.R` (modify) | `"potentiation"` stage in `upload_schema()`, `template_df()`, `validate_upload()` |
| `R/app-potentiation.R` (create) | `potentiation_ui()`, `potentiation_server()` |
| `R/app-run.R`, `R/app-intro.R` (modify) | navbar + rename + "which tab" paragraph |
| `R/validation-report.R`, `inst/validation/oracles.csv` (modify) | register dataset, `comparison` stage, oracle rows |
| `tests/testthat/fixtures/potentiation/` (create) | `potentiation_example_continuous.csv`, `extract_potentiation_fixture.R` |
| `inst/extdata/potentiation_example_continuous.csv` (create) | same CSV, bundled example for the app |
| tests | `test-potentiation.R`, `test-potentiation-summary.R`, additions to `test-app-io.R`, `test-plot-data.R`, `test-app-run.R`, `test-app-modules.R`, `test-validation-potentiation.R` |

---

### Task 1: Example fixture + engine `analyse_potentiation()`

**Files:**
- Create: `tests/testthat/fixtures/potentiation/extract_potentiation_fixture.R`, `tests/testthat/fixtures/potentiation/potentiation_example_continuous.csv`, `inst/extdata/potentiation_example_continuous.csv` (identical copy)
- Create: `R/potentiation.R`
- Test: `tests/testthat/test-potentiation.R`

**Interfaces:**
- Produces:
  - `analyse_potentiation(df, reference = NULL, response = c("continuous", "quantal"), alpha = 0.05, stressor1 = "Stressor 1")` → list with `reference` (chr), `response`, `alpha`, `series` (chr, file order, trimmed), `curves` (named list per series: `par` = named `c(max, slope, ec50)`, `objective`, `n`, `convergence`, `controls_borrowed` (lgl), `ec50_extrapolated` (lgl), `data` = the rows used incl. borrowed controls), `comparisons` (data.frame, one row per non-reference series in file order: `series, n, obj1, obj2, chi, df, p, ec50_reference, ec50_series, ec50_shared, ratio, significant, direction, controls_borrowed, ec50_extrapolated`; `direction` ∈ `"more toxic"`, `"less toxic"`, `"no significant influence"`), `shared` (named list per comparison: `par` = named `c(max_reference, slope_reference, max_series, slope_series, ec50)`, `objective`), `summary` (chr, from Task 2's `potentiation_summary()`; until Task 2 lands, `NA_character_`).
  - Exported (roxygen `@export`). `df` columns: `Series`, `C1`, and `Res` or `Exposed`/`Affected`.

- [ ] **Step 1: Write the extraction script and produce the fixture**

`extract_potentiation_fixture.R` takes the workbook path as its first argument (the workbook is NOT committed), reads sheet `75-25 drc CPF on FBSA IMI` with `readxl`, and writes the CSV (`Series,C1,Res`, `row.names = FALSE`). Layout (see spec "Reference values"): reference = `M30:N74`; Level 1 = `M101:N125`; Level 2 = `AG101:AH125`; Level 3 = `W101:X125` (the non-control rows of each series). It must `stopifnot()` that (a) `W30:X74` and `AG30:AH74` equal `M30:N74`, and (b) each series' control rows (`M81:N100`, `AG81:AH100`, `W81:X100`) equal the reference's controls `N30:N49`. Header comment: provenance (Sam's draft workbook, sheet, cells), and that labels are deliberately generic because the workbook's names are unreliable. Run it on `C:/Users/jelle/Downloads/DRCs different ratios FBSA with IMI and CPF.xlsx`; result: 120 rows (45 `Stressor 1` incl. 20 controls, 25 per level, no controls under the levels). Copy to `inst/extdata/`.

- [ ] **Step 2: Write the failing tests** in `tests/testthat/test-potentiation.R`

```r
ex <- function() read.csv(test_path("fixtures", "potentiation",
                                    "potentiation_example_continuous.csv"))

test_that("workbook level 1: separate vs shared EC50 reproduces Sam's sheet", {
  r <- analyse_potentiation(ex())
  c1 <- r$comparisons[r$comparisons$series == "Level 1", ]
  expect_equal(c1$n, 90)
  expect_equal(c1$obj1, 301550.972, tolerance = 1e-6)
  expect_equal(c1$obj2, 341687.917, tolerance = 1e-6)
  expect_equal(c1$chi, 11.2463, tolerance = 1e-4)
  expect_equal(c1$p, 0.000798, tolerance = 1e-2)
  expect_equal(c1$ec50_series, 0.13417, tolerance = 1e-3)
  expect_true(c1$significant); expect_equal(c1$direction, "less toxic")
  expect_true(c1$controls_borrowed)
})

test_that("workbook level 2: no influence, X2 floored at ~0", {
  c2 <- analyse_potentiation(ex())$comparisons
  c2 <- c2[c2$series == "Level 2", ]
  expect_equal(c2$obj1, 393465.224, tolerance = 1e-6)
  expect_gte(c2$chi, 0); expect_lt(c2$chi, 1e-3)
  expect_false(c2$significant); expect_equal(c2$direction, "no significant influence")
})

test_that("workbook level 3: corrected shared fit, extrapolated EC50", {
  c3 <- analyse_potentiation(ex())$comparisons
  c3 <- c3[c3$series == "Level 3", ]
  expect_equal(c3$obj1, 379874.993, tolerance = 1e-6)
  expect_equal(c3$obj2, 411358.37, tolerance = 1e-6)   # NOT the workbook's 399798.70 (cell error)
  expect_equal(c3$chi, 7.166, tolerance = 1e-3)
  expect_equal(c3$ec50_series, 0.3607, tolerance = 1e-3)
  expect_true(c3$significant); expect_equal(c3$direction, "less toxic")
  expect_true(c3$ec50_extrapolated)                    # max C1 in Level 3 is 0.1
})
```

Plus, each in its own `test_that()`:
- `"series and comparisons keep file order"`: `r$series == c("Stressor 1","Level 1","Level 2","Level 3")`, `r$comparisons$series == c("Level 1","Level 2","Level 3")`, `r$reference == "Stressor 1"`.
- `"recovery: no shift is not significant, a shift is detected with the right direction"`: build data with `simulate_single(c(max = 100, slope = 2, ec50 = 1), conc = c(0, 2^(-3:3)), cv = 0.05, seed = 1)` replicated 4× (bind with different seeds) as reference; series `"same"` from the same curve (different seed); series `"potent"` with `ec50 = 0.5`; series `"weak"` with `ec50 = 2`. Assert `same` not significant; `potent` significant and `"more toxic"`, `ratio < 1`; `weak` `"less toxic"`.
- `"a series with its own controls does not borrow"`: give `"potent"` its own `C1 == 0` rows → `controls_borrowed == FALSE` and `n` = its rows + reference rows.
- `"quantal path uses deviance"`: `simulate_single(..., response = "binary", group_size = 20, seed = )` for reference (ec50 1) and a series (ec50 0.4), rename to `Exposed/Affected` as returned; `analyse_potentiation(d, response = "quantal")`: `chi == pmax(0, obj2 - obj1)` exactly, series `"more toxic"`, all `curves[[i]]$par[["max"]] <= 1`.
- `"reference can be any series"` (Review Focus 4): `reference = "Level 2"` → `comparisons$series == c("Stressor 1","Level 1","Level 3")`.
- `"no controls anywhere still fits"` (Review Focus 2): drop all `C1 == 0` rows from the example → no error, all `controls_borrowed == FALSE`.
- `"flat series does not error and is flagged"` (Review Focus 3): a series whose `Res` is constant 100 at all doses → returns, `significant` is logical, `ec50_extrapolated` TRUE.
- `"obj1 <= obj2 always"`: for every comparison in every test result above.
- `"errors on bad input"`: unknown `reference` label → error mentioning the label; fewer than two series → error; a series with ≤ 3 rows (after borrowing) → error naming the series.

- [ ] **Step 3: Run to verify they fail**

Run: `Rscript -e ".libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths())); devtools::test(filter = 'potentiation$')"` (PowerShell)
Expected: FAIL — `could not find function "analyse_potentiation"`.

- [ ] **Step 4: Implement `R/potentiation.R`**

Internal helpers (names fixed so Task 2/3/4 can rely on the result, bodies are yours):
- `pot_objective(d, max, slope, ec50, response)` → numeric scalar (SS or deviance as in Global Constraints).
- `pot_fit_curve(d, response, starts)` → list(`par`, `objective`, `convergence`): 3-parameter fit; `starts` = list of named start vectors; default start as `fit_single()` (max = max response, slope 1, ec50 = median positive C1), clamped into bounds; best objective wins. Use `optim(method = "L-BFGS-B", control = list(parscale = pmax(abs(start), 1e-8), factr = 1e-9, maxit = 1000))` like `fit_single()`.
- `pot_fit_shared(d_ref, d_k, fit_ref, fit_k, response)` → list(`par`, `objective`, `convergence`): 5-parameter fit, seeds: shared ec50 = `fit_ref` ec50, `fit_k` ec50, and their geometric mean; max/slopes from run 1.
- After the shared fit, refit the series curve also from the shared solution's (`max_series`, `slope_series`, `ec50`) and keep the better, so `obj1 <= obj2`.
- `ec50_extrapolated`: ec50 > max(C1) or ec50 < min(C1[C1 > 0]) of that series' own rows.
- `direction`: as spec (ratio from run 1, only if `p < alpha`).
- `summary`: `NA_character_` for now (Task 2 fills it).
- Validation of `df` inside the function: as listed in the error test; `Series` via `trimws(as.character())`.

- [ ] **Step 5: Run to verify they pass** — same command. Expected: all pass. Also `devtools::document()`.

- [ ] **Step 6: Commit**

```bash
git add R/potentiation.R tests/testthat/test-potentiation.R tests/testthat/fixtures/potentiation inst/extdata/potentiation_example_continuous.csv NAMESPACE man/analyse_potentiation.Rd
git commit -m "feat: analyse_potentiation() — per-series EC50 shift test (Sam's workbook method)"
```

---

### Task 2: Conclusion text `potentiation_summary()`

**Files:**
- Create: `R/potentiation-summary.R`
- Modify: `R/potentiation.R` (fill `summary`)
- Test: `tests/testthat/test-potentiation-summary.R`

**Interfaces:**
- Consumes: `analyse_potentiation()` result (Task 1).
- Produces: `potentiation_summary(res, stressor1 = "Stressor 1")` → single character string (internal, `@keywords internal`). `format_p(p)` → `"p < 0.0001"` when `p < 1e-4`, else `paste("p =", format(signif(p, 2), scientific = FALSE))`. `format_num(x)` → `format(signif(x, 3), scientific = FALSE)`.

Exact templates (S = `stressor1`; per-significant item `"<series> (EC50 <ref> → <series ec50>, <format_p>[, extrapolated])"`):
- one comparison, significant: `"<S> is significantly <more|less> toxic with <item>."`
- one comparison, not significant: `"No significant influence of <series> on the toxicity of <S> (<format_p>)."`
- several, onset (spec definition): `"<S> is significantly <dir> toxic from <onset series> on: <item>; <item>."` (items = all significant series, file order)
- several, none significant: `"No significant influence on the toxicity of <S> at any level."`
- several, no onset: per direction present, `"<S> is significantly <dir> toxic with <item> and <item>"`, joined with `"; "`; then `"; no significant influence with <a>, <b>"` if any; end `". The effect does not increase consistently with the level."`

- [ ] **Step 1: Write the failing tests** — build `res` objects by hand (a minimal list with `series`, `reference`, `comparisons` data.frame) so each case is exact:
  - `"workbook example reads as the spec's example"`: on `analyse_potentiation(ex())`: summary contains `"less toxic with Level 1 (EC50"`, `"Level 3 (EC50"`, `"extrapolated"`, `"no significant influence with Level 2"`, and ends with `"The effect does not increase consistently with the level."`; does not contain `" from "`.
  - `"onset named when every later level is significant"`: sig pattern F, T, T same direction → contains `"from Level 2 on"`.
  - `"mixed directions list both and give no onset"`: Level 1 more toxic, Level 2 less toxic → contains `"more toxic with Level 1"` and `"less toxic with Level 2"`, no `" from "`.
  - `"single comparison"`: both templates exactly, e.g. `"No significant influence of Level 1 on the toxicity of Stressor 1 (p = 0.42)."`.
  - `"none significant"`: exact string.
  - `"stressor name is used"`: `stressor1 = "Copper"` → starts with `"Copper "`.
  - `format_p(0.000798) == "p = 0.0008"`, `format_p(0.00743) == "p = 0.0074"`, `format_p(5e-6) == "p < 0.0001"`.
- [ ] **Step 2: Run** `devtools::test(filter = 'potentiation-summary')` → FAIL (function missing).
- [ ] **Step 3: Implement** `potentiation_summary()`, `format_p()`, `format_num()`; set `res$summary <- potentiation_summary(res, stressor1)` at the end of `analyse_potentiation()`.
- [ ] **Step 4: Run** both filters (`potentiation-summary`, `potentiation$`) → PASS.
- [ ] **Step 5: Commit** `feat: conclusion text for the (de)potentiation analysis`.

---

### Task 3: Upload format for the `"potentiation"` stage

**Files:**
- Modify: `R/app-io.R` (`upload_schema()`, `template_df()`, `validate_upload()`)
- Test: `tests/testthat/test-app-io.R`

**Interfaces:**
- Produces: stage `"potentiation"`: schema continuous `c("Series", "C1", "Res")`, quantal `c("Series", "C1", "Affected", "Exposed")`. `template_df("potentiation", response)` → three series `"Stressor 1"`, `"Level 1"`, `"Level 2"`, each `C1 = c(0, 0.1, 0.3, 1, 3, 10)` with an illustrative declining response (levels steeper). `validate_upload(df, "potentiation", response)` → character errors.

- [ ] **Step 1: Write failing tests** in `test-app-io.R`:
  - schema for both responses (exact vectors above).
  - template round-trips: `validate_upload(template_df("potentiation", r), "potentiation", r)` is `character(0)` for both `r`.
  - `Series` is not reported as non-numeric; `C1`/`Res` still are.
  - errors (each asserted with `expect_match`): missing `Series`; empty/NA `Series` label (`"Series labels must not be empty."`); fewer than 2 series (`"at least two series"`); a series with fewer than 4 distinct `C1` values counting borrowed controls is NOT checked here (engine errors), but a series with fewer than 3 distinct positive `C1` values gives `"Series '<label>' has only <k> positive concentrations"`.
  - Review Focus 1: `Series = c(1, 1, …, 2, 2, …)` (numeric) validates without error, and `analyse_potentiation()` on it runs (labels `"1"`, `"2"`).
  - Review Focus 5: labels `"Level 1 "` and `"Level 1"` in one file count as one series (validation sees 2 series for a file with `"Stressor 1"` + those two spellings).
- [ ] **Step 2: Run** `devtools::test(filter = 'app-io')` → new tests FAIL.
- [ ] **Step 3: Implement.** Exclude `"Series"` from the numeric check; trim/`as.character` labels for counting.
- [ ] **Step 4: Run** → PASS.
- [ ] **Step 5: Commit** `feat(app): upload schema, template and validation for (de)potentiation data`.

---

### Task 4: Plot — curves per series

**Files:**
- Modify: `R/plot-data.R`, `R/plot.R`
- Test: `tests/testthat/test-plot-data.R`, `tests/testthat/test-plot.R`

**Interfaces:**
- Consumes: `analyse_potentiation()` result (`curves[[s]]$par`, `curves[[s]]$data`, `series`).
- Produces: `potentiation_curve_data(res)` (internal) → list(`curves` = data.frame `series, conc, response` on a shared log grid, 200 points per series, from the control offset (when any control) or the lowest positive dose to the highest dose over all series; `observed` = data.frame `series, conc, response, control` (lgl) with each series' own rows **without** borrowed controls (controls shown once, under the series they belong to); `control_x` = one decade below the lowest positive dose over all series, `NA` when no controls). `plot_potentiation(res, log_x = TRUE, xlab = "Concentration", ylab = "response")` (exported) → plotly object: one colour per series (points + line, same colour, legend = series label), controls as open circles at `control_x` with the `"0 (control)"` annotation, exactly as `plot_dose_response()` does.

- [ ] **Step 1: Failing tests**: on the example result — `unique(curves$series)` equals `res$series`; 200 rows per series; `min(curves$conc) == control_x`; `observed` has 120 rows (borrowed controls not duplicated); `sum(observed$control) == 20`; `plot_potentiation()` returns class `"plotly"` (skip if plotly missing).
- [ ] **Step 2: Run** filters `plot-data`, `plot$` → FAIL.
- [ ] **Step 3: Implement**; reuse the `obs_response()` convention for quantal (proportion).
- [ ] **Step 4: Run** → PASS; `devtools::document()`.
- [ ] **Step 5: Commit** `feat: plot the stressor-1 curve per series of stressor 2`.

---

### Task 5: The tab, navbar and rename

**Files:**
- Create: `R/app-potentiation.R`
- Modify: `R/app-run.R`, `R/app-intro.R`, `README.md:23`, `NEWS.md` (top, development section)
- Test: `tests/testthat/test-app-run.R`, `tests/testthat/test-app-modules.R`

**Interfaces:**
- Consumes: `analyse_potentiation()`, `potentiation_summary()`, `plot_potentiation()`, `template_df()/validate_upload()/read_upload()` (Tasks 1–4).
- Produces: `potentiation_ui(id)`, `potentiation_server(id, meta)` (internal). Navbar: `nav_panel("Introduction")`, `"Single Stressor"`, `"(De)potentiation"` → `potentiation_ui("potentiation")`, `"Mixture toxicity"` → `campaign_ui("campaign")`.

UI per spec "The tab": sidebar = `textInput` Stressor 1 name, Stressor 1 unit, Stressor 2 name; `radioButtons` response (Continuous/Quantal, values `continuous`/`quantal`); `downloadButton` template (`potentiation_<response>_template.csv`); `fileInput` CSV; errors; `selectInput` reference series (choices = series in file order, default first). Main = cards: plot (`plotly::plotlyOutput`), per-series table (max, slope, EC50, ratio vs reference, "extrapolated" flag), tests table (n, objective run 1/run 2 labelled "SS" or "deviance", X², p, verdict), controls note, conclusion. The bundled example (`inst/extdata/potentiation_example_continuous.csv`) is used until a file is uploaded, with a line saying so. Analysis re-runs on data/reference/response change (no Fit button).

Controls note (exact, per comparison that borrowed; S = stressor-1 name, k = number of borrowed control rows): `"Controls: the <k> control observations of <reference> are used for both curves of the comparison with <series>, as in the Excel workbook, so they are counted twice (n = <n>). This is the approach used here; it makes the test slightly more liberal than counting them once."` When nothing was borrowed: `"Controls: each series uses its own control observations."`

Introduction paragraph (add to the Welcome card): `"Which tab? If stressor 2 was applied at only one or a few fixed levels (for example with and without heat stress), use (De)potentiation: it tests whether stressor 2 changes the toxicity of stressor 1. If every stressor has its own dose-response series and their mixtures were dosed too, use Mixture toxicity."` and change "Multiple stressors tabs" → "Mixture toxicity tabs" in the existing sentence.

- [ ] **Step 1: Failing tests**: `test-app-run.R` — rename the test to `"app_ui exposes Introduction, Single Stressor, (De)potentiation and Mixture toxicity"`; assert all four labels present (`fixed = TRUE` for the parenthesis), `"Multiple stressors"` absent, existing `expect_false`s kept. `test-app-modules.R` — `potentiation_ui("p")` builds a tag containing `"Reference series"`; `shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), { ... })` with no upload: the analysis reactive returns a result whose `comparisons$series == c("Level 1","Level 2","Level 3")` (use the module's reactive name you define, e.g. `result()`).
- [ ] **Step 2: Run** filters `app-run`, `app-modules` → FAIL.
- [ ] **Step 3: Implement** the module, navbar, intro text; README line 23 "Multiple stressors tab" → "Mixture toxicity tab"; NEWS entry under development: the new tab (one paragraph: method, Excel-equal controls approach, generic) and the rename.
- [ ] **Step 4: Run** → PASS.
- [ ] **Step 5: Launch for a look** (controller does this, not a subagent): `Rscript -e ".libPaths(...); devtools::load_all(); run_app(port = 7788, launch.browser = FALSE)"` in the background; tell the user to open http://127.0.0.1:7788.
- [ ] **Step 6: Commit** `feat(app): (De)potentiation tab; rename Multiple stressors to Mixture toxicity`.

---

### Task 6: Validation report

**Files:**
- Modify: `R/validation-report.R` (`validation_datasets()`, `validation_extract()`), `inst/validation/oracles.csv`
- Create: `tests/testthat/test-validation-potentiation.R`
- Regenerate: `VALIDATION.md`, README badge block (via `tools/validation-report.R`)

**Interfaces:**
- Consumes: `analyse_potentiation()`; fixture path `potentiation/potentiation_example_continuous.csv`.
- Produces: dataset id `potentiation_example_continuous` (no `experiment`; title `"(De)potentiation example: stressor 1 with stressor 2 at three levels, continuous response"`; endpoint `"continuous response"`; `fit = function(df) analyse_potentiation(df)`). `validation_extract()` gains stage `"comparison"` where `parameter` is `"<series>:<field>"`, field a column of `res$comparisons` (e.g. `"Level 1:obj1"`).

Oracle rows (source_type `workbook`, source = `"Sam's draft workbook 'DRCs different ratios…', sheet 75-25…, cell <X>"`):

| parameter | reference | expected | tolerance | kind |
|---|---|---|---|---|
| `Level 1:obj1` | 301550.972 | 301550.972 | 1e-6 | match (M25) |
| `Level 1:obj2` | 341687.917 | 341687.917 | 1e-6 | match (P25) |
| `Level 1:chi` | 11.2463 | 11.2463 | 1e-4 | match (M27) |
| `Level 2:obj1` | 393465.224 | 393465.224 | 1e-6 | match (AG25) |
| `Level 2:obj2` | 393465.167 | 393465.167 | 1e-6 | match (AJ25) |
| `Level 2:chi` | -1.30e-5 | NA | NA | bounded (AG27; note: ≈ 0, asserted `0 <= chi < 1e-3` in the test) |
| `Level 3:obj1` | 379874.993 | 379874.993 | 1e-6 | match (W25) |
| `Level 3:obj2` | 399798.700 | 411358.37 | 1e-6 | divergent (Z25; note names the cell error AC31:AC74 `$P$17:$P$19` → `$Z$17:$Z$19`) |
| `Level 3:chi` | 4.6007 | 7.166 | 1e-3 | divergent (W27; same note) |

All rows: `stage = comparison`. Labels generic, e.g. `"Level 1 — residual SS, separate curves (run 1)"`.

- [ ] **Step 1: Failing test** `test-validation-potentiation.R`: `expect_matches_reference("potentiation_example_continuous")` and the bounded check on `Level 2:chi`.
- [ ] **Step 2: Run** filter `validation-potentiation` → FAIL (unregistered dataset).
- [ ] **Step 3: Implement** registration, `comparison` extraction, oracle rows.
- [ ] **Step 4: Run** → PASS.
- [ ] **Step 5: Regenerate the report** (controller, background — it refits the ternary, a few minutes): `Rscript tools/validation-report.R validation-report.html`; check the new section renders, the score line counts the new ✓ rows, and the ≠ footnote reads correctly.
- [ ] **Step 6: Commit** `feat: (de)potentiation workbook comparison in the validation report` (incl. `VALIDATION.md`, `README.md`).
