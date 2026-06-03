# Ternary app tab (T3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `Ternary Mixture` tab to the mixdra Shiny app that reuses the binary tab's UX skeleton (sidebar + three stage-cards + table-as-hub + click-to-inspect) and drives the existing staged Advanced-S/A engine.

**Architecture:** Three interactive single-curve panels (Stage 1) feed a frozen base into `analyse_ternary(base = …)` (a small new engine arg). Stage 2 is the per-ratio A4 hub table (Overall row + one row per ternary ratio); Stage 3 inspects the isoplane / ΣTU plots plus the selected ratio's effect-size readout. CA + continuous only; Quantal/IA radios shown but disabled.

**Tech Stack:** R, Shiny, bslib, plotly, DT, testthat. Reuses `curve_fit_ui/server`, `analyse_ternary`, `ternary_effect_table`, `plot_isoplane`, `plot_sigma_tu`.

**Spec:** `docs/superpowers/specs/2026-06-03-ternary-app-tab-design.md`

**Test convention (project memory):** run FILTERED tests only, e.g. `devtools::test(filter = "fit-ternary-asa")`. Never the full suite.

---

## File Structure

- `R/fit-ternary-asa.R` — MODIFY: add `base =` arg to `fit_ternary_asa()` + `analyse_ternary()`.
- `R/app-io.R` — MODIFY: `upload_schema`, `template_df`, `validate_upload` gain a "ternary" stage; add `marginal_df3()`, `assemble_curve_params3()`.
- `R/app-ternary.R` — CREATE: `ternary_ui()`, `ternary_server()`, plus the `disable_radio_js()` helper.
- `R/app-run.R` — MODIFY: add the nav panel + server wiring.
- `R/app-intro.R` — MODIFY: add a `chem3` name field.
- `inst/extdata/ternary_ca_fbsa_cpf_imi_continuous.csv` — CREATE (copy of the validated fixture).
- `tests/testthat/test-fit-ternary-asa.R` — MODIFY: `base =` tests.
- `tests/testthat/test-app-io.R` — MODIFY: ternary I/O tests.
- `tests/testthat/test-app-modules.R` — MODIFY: ternary module smoke test + bundled-example test.
- `tests/testthat/test-app-run.R` — MODIFY: assert the Ternary panel is wired.

---

## Task 1: Engine `base =` argument

**Files:**
- Modify: `R/fit-ternary-asa.R` (`fit_ternary_asa`, `analyse_ternary`)
- Test: `tests/testthat/test-fit-ternary-asa.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-fit-ternary-asa.R`:

```r
test_that("analyse_ternary uses a supplied base verbatim and skips the Stage-1 fit", {
  skip_on_cran()
  df <- utils::read.csv(testthat::test_path(
    "fixtures", "ternary", "fbsa_cpf_imi", "ternary_fbsa_cpf_imi_continuous.csv"))
  base <- c(max = 872.2, slope1 = 4.674, slope2 = 13, slope3 = 3.629,
            ec50_1 = 0.1275, ec50_2 = 5.58, ec50_3 = 0.575)
  res <- analyse_ternary(df, reference = "CA", response = "continuous", base = base)
  # base passed through unchanged (not re-fit from the singles)
  expect_equal(res$base[["max"]], 872.2)
  expect_equal(res$base[["ec50_1"]], 0.1275)
  expect_null(res$fits$base)                 # no Stage-1 fit object when base supplied
  expect_setequal(names(res$pairwise), c("A1", "A2", "A3"))
  expect_true(is.numeric(res$A4_overall))
  expect_gt(nrow(res$individual), 0)
})

test_that("analyse_ternary errors when the supplied base is missing a param", {
  df <- data.frame(C1 = c(0, 1, 1), C2 = c(0, 1, 0), C3 = c(0, 1, 1), Res = c(100, 20, 30))
  expect_error(analyse_ternary(df, base = c(max = 800, slope1 = 4)), "missing")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-fit-ternary-asa.R', filter='base')"`
Expected: FAIL — `unused argument (base = …)`.

- [ ] **Step 3: Implement the `base =` arg**

In `R/fit-ternary-asa.R`, change the `fit_ternary_asa()` signature and Stage-1 block.

Signature:

```r
fit_ternary_asa <- function(df, reference = "CA", response = "continuous",
                            base = NULL, lower = NULL, upper = NULL,
                            n_starts = 1, time_limit = 30) {
```

Replace the current Stage-1 block:

```r
  # Stage 1 - base curve params from singles + control.
  singles <- df[cls %in% c("control", "single"), , drop = FALSE]
  seed <- seed_from_singles(df, response)
  f1 <- fit_model(singles, reference, "reference", response, start = seed,
                  lower = lower, upper = upper, n_starts = n_starts,
                  time_limit = time_limit)
  base_par <- f1$par[base_params]
```

with:

```r
  # Stage 1 - base curve params: supplied (frozen, from the app's three single
  # curves) or fit from the singles + control. A supplied `base` is held fixed
  # downstream exactly as an auto-fit one would be.
  if (is.null(base)) {
    singles <- df[cls %in% c("control", "single"), , drop = FALSE]
    seed <- seed_from_singles(df, response)
    f1 <- fit_model(singles, reference, "reference", response, start = seed,
                    lower = lower, upper = upper, n_starts = n_starts,
                    time_limit = time_limit)
    base_par <- f1$par[base_params]
  } else {
    miss <- setdiff(base_params, names(base))
    if (length(miss))
      stop("fit_ternary_asa: `base` is missing param(s): ",
           paste(miss, collapse = ", "))
    base_par <- base[base_params]
    f1 <- NULL
  }
```

Then update `analyse_ternary()`'s signature and forwarding:

```r
analyse_ternary <- function(df, reference = "CA",
                            response = c("continuous", "binary"),
                            base = NULL, lower = NULL, upper = NULL,
                            n_starts = 1, time_limit = 30) {
  response <- match.arg(response)
  if (!all(c("C1", "C2", "C3") %in% names(df)))
    stop("analyse_ternary requires C1, C2 and C3 columns")
  fit_ternary_asa(df, reference = reference, response = response, base = base,
                  lower = lower, upper = upper, n_starts = n_starts,
                  time_limit = time_limit)
}
```

Also update the roxygen for both: add `@param base Optional named curve params (`max, slope1-3, ec50_1-3`); when supplied the Stage-1 singles fit is skipped and these are held fixed. Defaults to `NULL` (fit from singles).`

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-fit-ternary-asa.R')"`
Expected: PASS (new `base` tests + all pre-existing tests in the file still green — `base = NULL` is the unchanged path).

- [ ] **Step 5: Regenerate docs**

Run: `Rscript -e "devtools::document()"`
Expected: `man/analyse_ternary.Rd` updated with the `base` param; no errors.

- [ ] **Step 6: Commit**

```bash
git add R/fit-ternary-asa.R tests/testthat/test-fit-ternary-asa.R man/analyse_ternary.Rd
git commit -m "feat(ternary): analyse_ternary accepts a frozen base curve set"
```

---

## Task 2: I/O — `upload_schema` + `template_df` for ternary

**Files:**
- Modify: `R/app-io.R` (`upload_schema`, `template_df`)
- Test: `tests/testthat/test-app-io.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-app-io.R`:

```r
test_that("upload_schema knows the ternary stage", {
  expect_equal(upload_schema("ternary", "continuous"), c("C1", "C2", "C3", "Res"))
  expect_equal(upload_schema("ternary", "quantal"),
               c("C1", "C2", "C3", "Affected", "Exposed"))
})

test_that("ternary template has the schema columns and spans every tier", {
  d <- template_df("ternary", "continuous")
  expect_equal(names(d), c("C1", "C2", "C3", "Res"))
  expect_true(all(vapply(d, is.numeric, logical(1))))
  cls <- as.character(classify_rows(d))
  expect_true(all(c("control", "single", "binary", "ternary") %in% cls))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-io.R', filter='ternary')"`
Expected: FAIL — `'arg' should be one of "single", "binary"` (match.arg) / wrong columns.

- [ ] **Step 3: Implement**

In `R/app-io.R`, replace `upload_schema()` body:

```r
upload_schema <- function(stage, response) {
  stage <- match.arg(stage, c("single", "binary", "ternary"))
  response <- match.arg(response, c("continuous", "quantal"))
  if (stage == "single") {
    if (response == "continuous") c("Conc", "Res") else c("Conc", "Affected", "Exposed")
  } else if (stage == "binary") {
    if (response == "continuous") c("C1", "C2", "Res") else c("C1", "C2", "Affected", "Exposed")
  } else {
    if (response == "continuous") c("C1", "C2", "C3", "Res")
    else c("C1", "C2", "C3", "Affected", "Exposed")
  }
}
```

In `template_df()`, add a `ternary` branch. Replace the final `else { … }` (the binary branch) with an `else if (stage == "binary") { … }` keeping its existing body, then append:

```r
  } else {  # ternary: all tiers so the staged fit has data at every stage
    rows <- rbind(
      data.frame(C1 = 0,           C2 = 0,           C3 = 0),            # control
      data.frame(C1 = c(0.1, 0.3, 1), C2 = 0,        C3 = 0),           # chem1 single
      data.frame(C1 = 0,           C2 = c(0.1, 0.3, 1), C3 = 0),        # chem2 single
      data.frame(C1 = 0,           C2 = 0,           C3 = c(0.1, 0.3, 1)), # chem3 single
      data.frame(C1 = c(0.5, 0.5, 0), C2 = c(0.5, 0, 0.5), C3 = c(0, 0.5, 0.5)), # 3 binaries
      data.frame(C1 = c(0.2, 0.4, 0.6), C2 = c(0.2, 0.4, 0.6), C3 = c(0.2, 0.4, 0.6)) # 1:1:1 ternary
    )
    tot <- rows$C1 + rows$C2 + rows$C3
    if (response == "continuous") {
      cbind(rows, Res = round(100 / (1 + tot), 1))          # illustrative decline
    } else {
      cbind(rows, Affected = pmin(round(10 * tot / (1 + tot)), 10), Exposed = 10)
    }
  }
```

(Leave the `single` and `binary` branches unchanged.)

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-io.R', filter='ternary')"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat(io): ternary upload_schema + template_df"
```

---

## Task 3: I/O — `validate_upload` for ternary

**Files:**
- Modify: `R/app-io.R` (`validate_upload`)
- Test: `tests/testthat/test-app-io.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-app-io.R`:

```r
test_that("validate_upload accepts a valid ternary file", {
  good <- data.frame(C1 = c(0, 1, 0, 0, 1), C2 = c(0, 0, 1, 0, 1),
                     C3 = c(0, 0, 0, 1, 1), Res = c(100, 60, 60, 60, 20))
  expect_length(validate_upload(good, "ternary", "continuous"), 0)
})

test_that("validate_upload flags a ternary file with no ternary rows", {
  noternary <- data.frame(C1 = c(0, 1, 0), C2 = c(0, 0, 1),
                          C3 = c(0, 0, 0), Res = c(100, 60, 60))
  expect_match(paste(validate_upload(noternary, "ternary", "continuous"), collapse = " "),
               "ternary rows")
})

test_that("validate_upload catches a negative C3", {
  bad <- data.frame(C1 = c(0, 1), C2 = c(0, 1), C3 = c(-1, 1), Res = c(100, 20))
  expect_match(paste(validate_upload(bad, "ternary", "continuous"), collapse = " "),
               ">= 0|negative|0")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-io.R', filter='ternary')"`
Expected: FAIL — `validate_upload` does not yet accept `"ternary"` / no ternary-rows check.

- [ ] **Step 3: Implement**

In `R/app-io.R`, edit `validate_upload()`. The function already calls `upload_schema(stage, response)` (now ternary-aware). Two changes:

a) Add `C3` to the concentration columns checked. Replace:

```r
  conc_cols <- intersect(c("Conc", "C1", "C2"), present)
```

with:

```r
  conc_cols <- intersect(c("Conc", "C1", "C2", "C3"), present)
```

b) Before `errs` is returned (just above the final `errs`), add the ternary-rows check:

```r
  if (stage == "ternary" && all(c("C1", "C2", "C3") %in% present) &&
      all(vapply(df[c("C1", "C2", "C3")], is.numeric, logical(1)))) {
    nz <- rowSums(as.matrix(df[c("C1", "C2", "C3")]) > 0)
    if (!any(nz == 3, na.rm = TRUE))
      errs <- c(errs, paste0("No ternary rows (all of C1, C2, C3 > 0); the ",
                             "per-ratio A4 step needs at least one ternary mixture."))
  }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-io.R')"`
Expected: PASS (new ternary tests + all existing binary/single validation tests still green).

- [ ] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat(io): ternary validate_upload (C3 + ternary-rows check)"
```

---

## Task 4: I/O — `marginal_df3` + `assemble_curve_params3`

**Files:**
- Modify: `R/app-io.R` (add two functions)
- Test: `tests/testthat/test-app-io.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-app-io.R`:

```r
test_that("marginal_df3 extracts a chemical's single series (other two == 0)", {
  df <- data.frame(C1 = c(0, 1, 2, 0, 0, 3),
                   C2 = c(0, 0, 0, 1, 2, 4),
                   C3 = c(0, 0, 0, 0, 0, 5),
                   Res = c(100, 60, 40, 70, 50, 10))
  m1 <- marginal_df3(df, 1)             # rows where C2 == 0 AND C3 == 0
  expect_equal(m1$C1, c(0, 1, 2))
  expect_equal(m1$Res, c(100, 60, 40))
  expect_false(any(c("C2", "C3") %in% names(m1)))

  m3 <- marginal_df3(df, 3)             # rows where C1 == 0 AND C2 == 0; C3 -> C1
  expect_equal(m3$C1, c(0))             # only the control row qualifies here
  expect_equal(m3$Res, c(100))
})

test_that("assemble_curve_params3 averages max and uses underscore ec50 names", {
  f1 <- list(par = c(max = 870, slope = 4.7, ec50 = 0.13))
  f2 <- list(par = c(max = 872, slope = 13,  ec50 = 5.58))
  f3 <- list(par = c(max = 874, slope = 3.6, ec50 = 0.57))
  p <- assemble_curve_params3(f1, f2, f3)
  expect_equal(names(p), c("max", "slope1", "slope2", "slope3",
                           "ec50_1", "ec50_2", "ec50_3"))
  expect_equal(unname(p[["max"]]), 872)         # mean(870, 872, 874)
  expect_equal(unname(p[["slope2"]]), 13)
  expect_equal(unname(p[["ec50_3"]]), 0.57)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-io.R', filter='marginal_df3|assemble_curve_params3')"`
Expected: FAIL — `could not find function "marginal_df3"`.

- [ ] **Step 3: Implement**

Add to `R/app-io.R` (after `marginal_df` and `assemble_curve_params`):

```r
#' One chemical's single-compound series from a ternary frame
#'
#' Keeps rows where the OTHER TWO chemicals are 0 (so the shared control row is
#' included), drops their columns, and renames this chemical's concentration
#' column to `C1` — the shape a single-chemical fitter expects. The ternary
#' analogue of [marginal_df()].
#' @param df Ternary engine data frame (`C1`, `C2`, `C3`, response columns).
#' @param chem 1, 2 or 3 — which chemical's marginal series to extract.
#' @return A data frame with `C1` and the response columns.
#' @keywords internal
marginal_df3 <- function(df, chem) {
  cols  <- paste0("C", 1:3)
  this  <- paste0("C", chem)
  other <- setdiff(cols, this)
  keep  <- df[[other[1]]] == 0 & df[[other[2]]] == 0
  out   <- df[keep, , drop = FALSE]
  out[other] <- NULL
  names(out)[names(out) == this] <- "C1"
  rownames(out) <- NULL
  out
}

#' Frozen ternary base-parameter vector from three single-chemical fits
#'
#' Builds the named vector [analyse_ternary()] holds fixed as its `base`: a
#' shared `max` (mean of the three per-chemical fits, matching the engine's
#' seeding) plus per-chemical `slope1/2/3` and `ec50_1/2/3`. Names match the
#' ternary registry's base parameters (underscore `ec50_i`, unlike binary's
#' `ec501`).
#' @param fit1,fit2,fit3 Single-fit results (each `par = c(max, slope, ec50)`).
#' @return A named numeric vector: `max`, `slope1-3`, `ec50_1-3`.
#' @keywords internal
assemble_curve_params3 <- function(fit1, fit2, fit3) {
  c(max    = mean(c(fit1$par[["max"]], fit2$par[["max"]], fit3$par[["max"]])),
    slope1 = fit1$par[["slope"]],
    slope2 = fit2$par[["slope"]],
    slope3 = fit3$par[["slope"]],
    ec50_1 = fit1$par[["ec50"]],
    ec50_2 = fit2$par[["ec50"]],
    ec50_3 = fit3$par[["ec50"]])
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-io.R', filter='marginal_df3|assemble_curve_params3')"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat(io): marginal_df3 + assemble_curve_params3 for the ternary tab"
```

---

## Task 5: Bundle the ternary example dataset

**Files:**
- Create: `inst/extdata/ternary_ca_fbsa_cpf_imi_continuous.csv`
- Test: `tests/testthat/test-app-modules.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-app-modules.R`:

```r
test_that("a bundled example dataset ships and validates as ternary continuous", {
  ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled ternary example not installed")
  df <- read_upload(ex)
  expect_length(validate_upload(df, "ternary", "continuous"), 0)
  expect_true(all(c("C1", "C2", "C3", "Res") %in% names(df)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', filter='bundled example dataset ships and validates as ternary')"`
Expected: FAIL/SKIP — file not installed (skip means the assertion never runs; treat skip as "not yet done").

- [ ] **Step 3: Copy the fixture into `inst/extdata`**

Run (PowerShell):

```powershell
Copy-Item "tests/testthat/fixtures/ternary/fbsa_cpf_imi/ternary_fbsa_cpf_imi_continuous.csv" "inst/extdata/ternary_ca_fbsa_cpf_imi_continuous.csv"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', filter='bundled example dataset ships and validates as ternary')"`
Expected: PASS (`devtools::load_all` makes `system.file('extdata', …)` resolve to `inst/extdata`).

- [ ] **Step 5: Commit**

```bash
git add inst/extdata/ternary_ca_fbsa_cpf_imi_continuous.csv tests/testthat/test-app-modules.R
git commit -m "feat(app): bundle ternary FBSA/CPF/IMI example dataset"
```

---

## Task 6: `ternary_ui` + `disable_radio_js`

**Files:**
- Create: `R/app-ternary.R` (UI half + helper)
- Test: `tests/testthat/test-app-modules.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-app-modules.R`:

```r
test_that("ternary_ui builds the sidebar, three curve panels, and the three stages", {
  html <- as.character(ternary_ui("ternary"))
  # sidebar controls + disabled-radio script
  expect_match(html, "ternary-response", fixed = TRUE)
  expect_match(html, "ternary-reference", fixed = TRUE)
  expect_match(html, "coming soon", fixed = TRUE)
  expect_match(html, "prop('disabled', true)", fixed = TRUE)
  # Stage 1: three curve panels
  expect_match(html, "ternary-chem1-autofit", fixed = TRUE)
  expect_match(html, "ternary-chem2-autofit", fixed = TRUE)
  expect_match(html, "ternary-chem3-autofit", fixed = TRUE)
  # Stage 2: fit button + hub table
  expect_match(html, "ternary-fit_asa", fixed = TRUE)
  expect_match(html, "ternary-hub", fixed = TRUE)
  # Stage 3: plots + effect readout
  expect_match(html, "ternary-isoplane", fixed = TRUE)
  expect_match(html, "ternary-sigma_tu", fixed = TRUE)
  expect_match(html, "ternary-effect", fixed = TRUE)
  # no alpha / joint controls (deliberately dropped vs binary)
  expect_false(grepl("ternary-alpha", html, fixed = TRUE))
  expect_false(grepl("ternary-optimize_all", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', filter='ternary_ui builds')"`
Expected: FAIL — `could not find function "ternary_ui"`.

- [ ] **Step 3: Create `R/app-ternary.R` with the helper + UI**

```r
# Ternary Mixture stage. Reuses the binary tab's UX skeleton: a sidebar, three
# single-curve panels (Stage 1), a per-ratio A4 "hub" table (Stage 2), and an
# inspect stage (Stage 3 -- isoplane + Sigma-TU + the selected ratio's effect).
# The science differs from binary: there is no SA/DR/DL model selection. The
# staged Advanced-S/A fit (analyse_ternary) is the method; its punchline is the
# overall-A4 vs per-ratio-A4 contrast (averaging-out). CA + continuous only.

#' Client-side JS to disable one option of a Shiny radio group
#'
#' Shiny renders radio options as `<input name=... value=...>`; jQuery (bundled
#' with Shiny) flips the unsupported one to disabled on load. Used to show the
#' Quantal / IA options greyed-out without adding a JS dependency.
#' @param input_name The namespaced input id (e.g. `ns("response")`).
#' @param value The option value to disable (e.g. `"quantal"`).
#' @keywords internal
disable_radio_js <- function(input_name, value) {
  shiny::tags$script(shiny::HTML(sprintf(
    "$(function(){$('input[name=\"%s\"][value=\"%s\"]').prop('disabled', true);});",
    input_name, value)))
}

#' Ternary Mixture stage UI
#' @param id Module id.
#' @keywords internal
ternary_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous",
                            "Quantal (coming soon)" = "quantal"),
                          selected = "continuous"),
      disable_radio_js(ns("response"), "quantal"),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA) (coming soon)" = "IA"),
                          selected = "CA"),
      disable_radio_js(ns("reference"), "IA"),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example dataset (FBSA + CPF + IMI, continuous) is loaded until you upload your own.")),
      shiny::uiOutput(ns("errors")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced fitting options",
          shiny::helpText("Apply to the staged Advanced-S/A fit."),
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("time_limit"), "time_limit (s/stage)", value = 30, min = 1),
          shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
          shiny::uiOutput(ns("thorough_note"))
        )
      )
    ),

    # Stage 1 -- three single-chemical curve panels.
    bslib::card(
      bslib::card_header("Stage 1 · Single curves"),
      shiny::p("Fit each chemical's dose-response curve (Autofit or Simulate). ",
               "These curves are held fixed as the base for the staged ",
               "Advanced-S/A fit below, which appears once all three are fitted."),
      shiny::div(shiny::h5(shiny::textOutput(ns("chem1_title"))),
                 curve_fit_ui(ns("chem1"))),
      shiny::div(class = "mt-4",
                 shiny::h5(shiny::textOutput(ns("chem2_title"))),
                 curve_fit_ui(ns("chem2"))),
      shiny::div(class = "mt-4",
                 shiny::h5(shiny::textOutput(ns("chem3_title"))),
                 curve_fit_ui(ns("chem3"))),
      shiny::uiOutput(ns("reveal_note"))
    ),

    # Stage 2 -- staged Advanced-S/A fit + the per-ratio A4 hub table.
    bslib::card(
      bslib::card_header("Stage 2 · Advanced S/A"),
      shiny::p("Fit all three single curves in Stage 1, then ",
               shiny::tags$b("Fit Advanced S/A"), " to run the staged fit ",
               "(base → pairwise A1/A2/A3 from the binaries → overall A4 ",
               "→ an individual A4 per ternary ratio)."),
      shiny::div(shiny::actionButton(ns("fit_asa"), "Fit Advanced S/A",
                                     class = "btn-primary")),
      shiny::uiOutput(ns("base_pairwise")),
      shiny::helpText(
        "Each row is one ternary ratio; the ", shiny::tags$b("Overall"),
        " row is the single pooled A4. When per-ratio A4 values differ in sign ",
        "but the overall A4 is near zero, real interaction is being averaged ",
        "away — that contrast is the point of the ternary analysis. ",
        "Click a ratio row to inspect it below."),
      DT::DTOutput(ns("hub"))
    ),

    # Stage 3 -- inspect.
    bslib::card(
      bslib::card_header("Stage 3 · Inspect"),
      bslib::layout_columns(
        bslib::card(bslib::card_header("EC50 isoplane (simplex)"),
                    plotly::plotlyOutput(ns("isoplane"), height = "520px")),
        bslib::card(bslib::card_header("ΣTU vs z"),
                    plotly::plotlyOutput(ns("sigma_tu"), height = "520px"))
      ),
      bslib::card(bslib::card_header("Selected ratio — effect size"),
                  shiny::uiOutput(ns("effect")))
    )
  )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', filter='ternary_ui builds')"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app-ternary.R tests/testthat/test-app-modules.R
git commit -m "feat(app): ternary_ui (sidebar + 3 curves + hub + inspect)"
```

---

## Task 7: `ternary_server`

**Files:**
- Modify: `R/app-ternary.R` (add the server function)
- Test: covered by Task 9 (module smoke test)

- [ ] **Step 1: Append the server to `R/app-ternary.R`**

```r
#' Ternary Mixture stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
ternary_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() "ternary_continuous_template.csv",
      content  = function(file)
        utils::write.csv(template_df("ternary", "continuous"), file, row.names = FALSE)
    )

    output$thorough_note <- shiny::renderUI({
      if (isTRUE(input$thorough))
        shiny::div(class = "text-warning",
                   shiny::tags$small("Multi-start fitting may take several minutes."))
    })

    chem_title <- function(field, n) shiny::renderText({
      nm <- meta[[field]]; if (!is.null(nm) && nzchar(nm)) nm else paste("Chemical", n)
    })
    output$chem1_title <- chem_title("chem1", 1)
    output$chem2_title <- chem_title("chem2", 2)
    output$chem3_title <- chem_title("chem3", 3)

    # Data source: upload, else the bundled FBSA/CPF/IMI example.
    upload_path <- shiny::reactive({
      if (!is.null(input$file)) return(input$file$datapath)
      ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
      if (nzchar(ex)) ex else NULL
    })
    parsed <- shiny::reactive({ shiny::req(upload_path()); read_upload(upload_path()) })
    errs <- shiny::reactive({
      shiny::req(upload_path()); validate_upload(parsed(), "ternary", "continuous")
    })
    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    engine_df <- shiny::reactive({
      shiny::req(upload_path(), length(errs()) == 0)
      to_engine_df(parsed(), "ternary")
    })
    m1 <- shiny::reactive(marginal_df3(engine_df(), 1))
    m2 <- shiny::reactive(marginal_df3(engine_df(), 2))
    m3 <- shiny::reactive(marginal_df3(engine_df(), 3))

    fit1 <- curve_fit_server("chem1", fit_df = m1, meta = meta, chem_field = "chem1")
    fit2 <- curve_fit_server("chem2", fit_df = m2, meta = meta, chem_field = "chem2")
    fit3 <- curve_fit_server("chem3", fit_df = m3, meta = meta, chem_field = "chem3")

    frozen <- shiny::reactive(!is.null(fit1()) && !is.null(fit2()) && !is.null(fit3()))
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2(), fit3())
      assemble_curve_params3(fit1(), fit2(), fit3())
    })

    n_starts_eff <- shiny::reactive(
      if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts)

    output$reveal_note <- shiny::renderUI({
      if (!isTRUE(frozen()))
        shiny::div(class = "text-muted",
                   shiny::tags$small("Fit all three single curves (Autofit or ",
                                     "Simulate) to enable the Advanced-S/A fit."))
    })

    # The staged result; cleared when the curves / reference / response change.
    asa_res <- shiny::reactiveVal(NULL)
    sel_row <- shiny::reactiveVal(integer(0))
    shiny::observeEvent(list(curve_params(), input$reference, input$response), {
      asa_res(NULL); sel_row(integer(0))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$fit_asa, {
      if (!isTRUE(frozen())) {
        shiny::showNotification(
          "Fit all three single curves first (Autofit / Simulate in Stage 1).",
          type = "warning")
        return()
      }
      res <- shiny::withProgress(message = "Fitting Advanced S/A", value = 0.5, {
        tryCatch(
          analyse_ternary(engine_df(), reference = "CA", response = "continuous",
                          base = curve_params(), n_starts = n_starts_eff(),
                          time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
            NULL
          })
      })
      asa_res(res)
      sel_row(1L)                       # default-select the Overall row
    })

    shiny::observeEvent(input$hub_rows_selected, {
      sel_row(input$hub_rows_selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    effect_tbl <- shiny::reactive({
      res <- asa_res(); shiny::req(res)
      ternary_effect_table(res, engine_df())
    })

    # Hub: Overall row + one row per ternary ratio (joined with the effect table).
    hub_df <- shiny::reactive({
      res <- asa_res(); shiny::req(res)
      ind <- res$individual
      eff <- effect_tbl()
      m <- merge(ind, eff[, c("ratio", "pred_CA", "pred_SA", "pred_ASA", "a4_effect")],
                 by = "ratio", all.x = TRUE, sort = FALSE)
      per <- data.frame(
        Ratio    = sprintf("%.2f:%.2f:%.2f", m$C1, m$C2, m$C3),
        n        = m$n,
        A4       = round(m$A4, 4),
        pred_CA  = round(m$pred_CA, 1),
        pred_SA  = round(m$pred_SA, 1),
        pred_ASA = round(m$pred_ASA, 1),
        a4_effect = round(m$a4_effect, 1),
        .key = m$ratio, stringsAsFactors = FALSE, check.names = FALSE)
      overall <- data.frame(
        Ratio = "Overall", n = sum(m$n), A4 = round(res$A4_overall, 4),
        pred_CA = NA_real_, pred_SA = NA_real_, pred_ASA = NA_real_,
        a4_effect = NA_real_, .key = NA_character_,
        stringsAsFactors = FALSE, check.names = FALSE)
      rbind(overall, per)
    })

    selected_ratio <- shiny::reactive({
      h <- hub_df(); r <- sel_row()
      if (!length(r) || r[[1]] < 1 || r[[1]] > nrow(h)) return(NULL)
      k <- h$.key[r[[1]]]
      if (is.na(k)) NULL else k
    })

    output$base_pairwise <- shiny::renderUI({
      res <- asa_res(); shiny::req(res)
      b <- res$base; p <- res$pairwise
      shiny::tags$p(
        shiny::tags$b("Base: "),
        sprintf("max %.1f | slope %.2f/%.2f/%.2f | EC50 %.3f/%.3f/%.3f",
                b[["max"]], b[["slope1"]], b[["slope2"]], b[["slope3"]],
                b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]]),
        shiny::tags$br(),
        shiny::tags$b("Pairwise: "),
        sprintf("A1 %.3f | A2 %.3f | A3 %.3f", p[["A1"]], p[["A2"]], p[["A3"]]))
    })

    output$hub <- DT::renderDT({
      h <- hub_df()
      hidekey <- which(names(h) == ".key") - 1L     # 0-based
      dt <- DT::datatable(
        h, rownames = FALSE,
        selection = list(mode = "single", target = "row",
                         selected = shiny::isolate(sel_row())),
        options = list(dom = "t", ordering = FALSE, scrollX = TRUE,
                       columnDefs = list(list(visible = FALSE, targets = hidekey))))
      DT::formatStyle(dt, "Ratio", valueColumns = "Ratio", target = "row",
                      fontWeight = DT::styleEqual("Overall", "bold"),
                      backgroundColor = DT::styleEqual("Overall", "#d8f0d8"))
    })

    output$isoplane <- plotly::renderPlotly({
      res <- asa_res(); shiny::req(res); plot_isoplane(res, engine_df())
    })
    output$sigma_tu <- plotly::renderPlotly({
      res <- asa_res(); shiny::req(res); plot_sigma_tu(res)
    })

    output$effect <- shiny::renderUI({
      res <- asa_res(); shiny::req(res)
      k <- selected_ratio()
      if (is.null(k))
        return(shiny::tags$p(shiny::tags$small(
          "Select a ratio row above to see its near-EC50 effect size. The ",
          "Overall row pools every ternary mixture.")))
      eff <- effect_tbl(); row <- eff[eff$ratio == k, , drop = FALSE]
      shiny::req(nrow(row) == 1)
      shiny::tags$p(
        shiny::tags$b("At the near-EC50 point "),
        sprintf("(C1=%.3g, C2=%.3g, C3=%.3g): ", row$C1, row$C2, row$C3),
        sprintf("CA %.1f → CA+S/A %.1f → CA+S/A+S/A %.1f  |  A4 effect %.1f",
                row$pred_CA, row$pred_SA, row$pred_ASA, row$a4_effect))
    })

    list(asa_res = asa_res, hub_df = hub_df, effect_tbl = effect_tbl,
         selected_ratio = selected_ratio, frozen = frozen)  # for testability
  })
}
```

- [ ] **Step 2: Sanity-load the package**

Run: `Rscript -e "devtools::load_all('.'); message('ok')"`
Expected: prints `ok` with no parse/collation errors.

- [ ] **Step 3: Commit**

```bash
git add R/app-ternary.R
git commit -m "feat(app): ternary_server (staged ASA fit + hub table + inspect)"
```

---

## Task 8: Wire the tab into the app + add `meta$chem3`

**Files:**
- Modify: `R/app-run.R` (`app_ui`, `app_server`)
- Modify: `R/app-intro.R` (`intro_ui`, `intro_server`)
- Test: `tests/testthat/test-app-run.R`, `tests/testthat/test-app-modules.R`

- [ ] **Step 1: Write the failing tests**

Add to `tests/testthat/test-app-run.R`:

```r
test_that("app_ui includes the Ternary Mixture panel", {
  expect_match(as.character(app_ui()), "ternary-fit_asa", fixed = TRUE)
})
```

Add to `tests/testthat/test-app-modules.R` (inside the existing intro test or as a new one):

```r
test_that("intro_server writes chem3 into the shared meta store", {
  meta <- shiny::reactiveValues()
  shiny::testServer(intro_server, args = list(meta = meta), {
    session$setInputs(chem1 = "CPF", chem2 = "FBSA", chem3 = "IMI")
    expect_equal(meta$chem3, "IMI")
  })
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-run.R'); testthat::test_file('tests/testthat/test-app-modules.R', filter='chem3')"`
Expected: FAIL — no Ternary panel; `meta$chem3` is NULL.

- [ ] **Step 3: Implement the wiring**

In `R/app-run.R`, `app_ui()` — add after the Binary panel:

```r
    bslib::nav_panel("Binary Mixture", binary_ui("binary")),
    bslib::nav_panel("Ternary Mixture", ternary_ui("ternary"))
```

In `R/app-run.R`, `app_server()` — add after the binary server:

```r
  binary_server("binary", meta)
  ternary_server("ternary", meta)
```

In `R/app-intro.R`, `intro_ui()` — add a chem3 field after chem2:

```r
      shiny::textInput(ns("chem2"), "Chemical 2 name"),
      shiny::textInput(ns("chem3"), "Chemical 3 name (ternary)"),
```

In `R/app-intro.R`, `intro_server()` — add inside the `observe`:

```r
      meta$chem2    <- input$chem2
      meta$chem3    <- input$chem3
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-run.R'); testthat::test_file('tests/testthat/test-app-modules.R', filter='chem3')"`
Expected: PASS. Also re-run the existing intro test to confirm no regression:
`Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', filter='intro')"` → PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app-run.R R/app-intro.R tests/testthat/test-app-run.R tests/testthat/test-app-modules.R
git commit -m "feat(app): wire Ternary Mixture tab + meta chem3"
```

---

## Task 9: Ternary module smoke test (testServer)

**Files:**
- Modify: `tests/testthat/test-app-modules.R`

- [ ] **Step 1: Write the test**

Add to `tests/testthat/test-app-modules.R`:

```r
test_that("ternary_server: three curves -> Fit Advanced S/A -> hub + selection", {
  skip_on_cran()
  ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled ternary example not installed")
  meta <- shiny::reactiveValues(chem1 = "CPF", chem2 = "FBSA", chem3 = "IMI")
  shiny::testServer(ternary_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA",
                      thorough = FALSE, n_starts = 1, time_limit = 30)
    expect_length(errs(), 0)                  # bundled example is valid
    expect_gt(nrow(engine_df()), 0)

    # Simulate the three single curves (values near the validated fit).
    session$setInputs(`chem1-val_max` = 872, `chem1-val_slope` = 4.67,
                      `chem1-val_ec50` = 0.1275, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 872, `chem2-val_slope` = 13,
                      `chem2-val_ec50` = 5.58, `chem2-simulate` = 1)
    session$setInputs(`chem3-val_max` = 872, `chem3-val_slope` = 3.63,
                      `chem3-val_ec50` = 0.575, `chem3-simulate` = 1)
    expect_true(isTRUE(frozen()))

    # Run the staged Advanced-S/A fit.
    session$setInputs(fit_asa = 1)
    res <- asa_res()
    expect_false(is.null(res))
    expect_setequal(names(res$pairwise), c("A1", "A2", "A3"))
    expect_gt(nrow(res$individual), 0)

    # Hub: Overall row + one row per ternary ratio; Overall selected by default.
    h <- hub_df()
    expect_equal(h$Ratio[1], "Overall")
    expect_equal(nrow(h), 1L + nrow(res$individual))
    expect_null(selected_ratio())             # row 1 (Overall) -> NULL ratio

    # Clicking a ratio row selects that ratio.
    session$setInputs(hub_rows_selected = 2)
    expect_equal(selected_ratio(), h$.key[2])

    # A curve edit clears the staged result (user re-runs Fit Advanced S/A).
    session$setInputs(`chem1-val_max` = 880, `chem1-simulate` = 2)
    expect_null(asa_res())
  })
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', filter='three curves')"`
Expected: PASS. (The staged fit over the ~419-row fixture with `n_starts = 1` runs in a few seconds.)

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-app-modules.R
git commit -m "test(app): ternary_server smoke test (fit + hub + selection)"
```

---

## Task 10: Final verification

**Files:** none (verification only)

- [ ] **Step 1: Run every test touched by this plan (filtered)**

Run:

```bash
Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-fit-ternary-asa.R'); testthat::test_file('tests/testthat/test-app-io.R'); testthat::test_file('tests/testthat/test-app-modules.R'); testthat::test_file('tests/testthat/test-app-run.R')"
```

Expected: all PASS, 0 failures. (Per project memory, do NOT run the full suite.)

- [ ] **Step 2: Confirm docs are current**

Run: `Rscript -e "devtools::document()"`
Expected: no uncommitted `man/` changes beyond Task 1's. If any appear, commit them:

```bash
git add man/ NAMESPACE
git commit -m "docs: regenerate man pages for ternary tab"
```

- [ ] **Step 3 (optional manual check): launch the app**

Run: `Rscript -e "devtools::load_all('.'); mixdra::run_app()"` and confirm the Ternary Mixture tab loads, the three curve panels fit, Fit Advanced S/A fills the hub table, and clicking a ratio updates the Stage-3 effect readout. (Manual; skip in headless CI.)

---

## Self-Review notes

- **Spec coverage:** §4 engine `base=` → Task 1; §5 I/O (`upload_schema`/`template_df` → Task 2; `validate_upload` → Task 3; `marginal_df3`/`assemble_curve_params3` → Task 4) ; §6 bundled example → Task 5; §7 UI → Task 6, server → Task 7; §8 wiring + `meta$chem3` → Task 8; §10 testing → Tasks 1–9 + Task 10. §3 decisions (disabled radios, hub table, frozen Stage 1) all realized in Tasks 6–7.
- **Out of scope (spec §11):** no IA/quantal path, no joint-refine card, no `alpha` — confirmed absent in Task 6's negative assertions.
- **Type consistency:** base names `max, slope1-3, ec50_1-3` (underscore) used identically in `assemble_curve_params3` (Task 4), the engine `base_params` (Task 1, from the registry), and the test seeds. Hub keys use `res$individual$ratio` joined to `effect_tbl$ratio` (both from `ternary_ratio_key`). Server reactive names (`asa_res`, `hub_df`, `effect_tbl`, `selected_ratio`, `frozen`) match the Task 9 test exactly.
