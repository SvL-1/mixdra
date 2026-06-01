# Staged Fitting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `analyse_mixture()` fit the curve parameters (`max`, `slope*`, `ec50*`) to single-compound data only, fix them, then fit only the interaction parameters (`a`, `b`, …) to the mixture data — and start the solver at `a = 0`, `b = 1`.

**Architecture:** Two-stage orchestration inside `analyse_mixture()`. Stage 1 fits the `reference` model to control + single-compound rows via the existing `fit_model()` (yielding one shared `max` plus per-chemical `slope`/`ec50`). Stage 2 calls `fit_model()` once per deviation with those curve params held `fixed`, so only the interaction params are free. Two small `fit_model()` changes support this: per-name interaction start defaults (`a=0`, `b*=1`), and a fast path when every parameter is fixed.

**Tech Stack:** R package (`mixdra`), `stats::optim`, `testthat` (edition 3), `roxygen2`/`devtools`.

**Spec:** `docs/superpowers/specs/2026-05-30-staged-fitting-design.md`

---

## ✅ STATUS: COMPLETED (2026-05-31, branch `mixdra-engine`, pushed to origin)

All 5 tasks implemented via subagent-driven execution (each with spec + code-quality review); full suite green (205 pass / 0 fail). Checkboxes below were not ticked during execution but every task is done and committed.

**Commits:** `52cbe3b` (a=0/b=1 starts) · `1aeb5ac` (all-fixed fast path) · `63c5370` (`fit_curve_from_singles`) · `a5e97ab`+`7596d8d` (staged `analyse_mixture`) · `7c02629` (re-baseline) · `1f3378d` (bounds-validation fix + extra tests). Doc-fix `6ec0258`.

**Deviations from this plan (as-built):**
- **Task 4 tests:** the planned "chooses reference on noiseless data" assertion was replaced with parameter assertions (`a≈0`) plus a known-interaction recovery test — on noiseless data the reference/SA objectives are ~1e-15 and the LR *selection* is numerically degenerate. Same change applied to `test-analyse.R` and the ternary test in Task 5.
- **Task 5 validation re-baseline:** the real quantal binary result **flipped `reference`→`DR`** (continuous still `DL`, weaker S/A). This was investigated and confirmed correct (see the spec's "Outcome" section), not edited away.
- **Post-review additions (`1f3378d`):** bound-name validation moved *above* the all-fixed fast path (it was being skipped); added a sparse-marginal fallback test and an all-fixed-bounds test.
- **Man pages:** `@keywords internal` functions DO get `.Rd` files (the plan's Step-4 assumption that they wouldn't was wrong); `init_start_par.Rd` / `fit_curve_from_singles.Rd` exist and are correct (unexported in NAMESPACE).

See the spec's **Outcome / as-built** section for the scientific result and the deferred follow-ups.

---

## Conventions for every task

- **Test command** (PowerShell or bash, runs one file):
  ```
  Rscript -e ".libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths())); suppressMessages(devtools::load_all('.')); testthat::test_file('tests/testthat/<FILE>', reporter='summary')"
  ```
  `load_all('.')` makes **internal** (non-exported) functions like `init_start_par`,
  `fit_curve_from_singles`, and `seed_from_singles` callable in tests.
- **Commit convention:** this repo does **not** add a Claude co-author trailer. Use the existing `feat:` / `test:` / `docs:` prefixes. Keep commit bodies short.
- Work on the current branch (`mixdra-engine`). Do not create a PR unless asked.

---

## Task 1: Interaction-parameter start defaults (`a=0`, `b-family=1`)

Extract the start-vector assembly in `fit_model()` into a testable internal helper that defaults `a` to 0 and every `b`-family param (`b`, `b1`, `b2`, `b3`) to 1.

**Files:**
- Modify: `R/fit.R` (currently builds the start vector at lines 44-47)
- Test: `tests/testthat/test-init-start.R` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-init-start.R`:

```r
test_that("init_start_par defaults a=0, b=1, base=0; start overrides", {
  p <- init_start_par(
    params = c("max", "slope1", "slope2", "ec501", "ec502", "a", "b"),
    extra  = c("a", "b"),
    start  = c(max = 800, ec501 = 0.08))
  expect_equal(unname(p[["a"]]), 0)        # a starts at the solver origin 0
  expect_equal(unname(p[["b"]]), 1)        # b starts at 1
  expect_equal(unname(p[["slope1"]]), 0)   # base param not in start -> 0
  expect_equal(unname(p[["max"]]), 800)    # start overrides
  expect_equal(unname(p[["ec501"]]), 0.08) # start overrides
})

test_that("init_start_par sets ternary b1/b2/b3 to 1 and a to 0", {
  p <- init_start_par(
    params = c("max", "slope1", "slope2", "slope3",
               "ec50_1", "ec50_2", "ec50_3", "a", "b1", "b2", "b3"),
    extra  = c("a", "b1", "b2", "b3"),
    start  = numeric(0))
  expect_equal(unname(p[c("b1", "b2", "b3")]), c(1, 1, 1))
  expect_equal(unname(p[["a"]]), 0)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run the test command with `<FILE>` = `test-init-start.R`.
Expected: FAIL — `could not find function "init_start_par"`.

- [ ] **Step 3: Add the helper and use it in `fit_model`**

In `R/fit.R`, **add** this helper above `fit_model` (just before the `#' Fit one mixture model` roxygen block):

```r
#' Assemble the full named starting vector for a model fit
#'
#' Base curve parameters default to 0 (always supplied via `start`); interaction
#' parameters default to their solver origin: `a` at 0 and every `b`-family
#' parameter (`b`, `b1`, `b2`, `b3`) at 1. Values present in `start` override the
#' defaults.
#' @keywords internal
init_start_par <- function(params, extra, start) {
  par <- stats::setNames(numeric(length(params)), params)  # all 0
  bfam <- intersect(extra, c("b", "b1", "b2", "b3"))
  if (length(bfam)) par[bfam] <- 1
  par[names(start)] <- start[names(start)]
  par
}
```

Then in `fit_model`, **replace** the current assembly:

```r
  # Assemble the full starting vector. Zero-initialising means any deviation
  # parameter (a, b, b1, b2, b3) not supplied in `start` defaults to 0.
  par <- stats::setNames(numeric(length(spec$params)), spec$params)
  par[names(start)] <- start[names(start)]
```

with:

```r
  # Assemble the full starting vector: interaction params default to their solver
  # origin (a = 0, b-family = 1); curve params come from `start`.
  par <- init_start_par(spec$params, spec$extra, start)
```

- [ ] **Step 4: Run the test to verify it passes**

Run the test command with `<FILE>` = `test-init-start.R`. Expected: PASS (2 tests).

- [ ] **Step 5: Regression-check the existing fit tests**

Run the test command for `test-fit.R`, then `test-fit-robust.R`, then `test-fit-enrich.R`.
Expected: all PASS (the all-zero default only changed for `b`-family params, which these tests do not pin).

- [ ] **Step 6: Commit**

```bash
git add R/fit.R tests/testthat/test-init-start.R
git commit -m "feat: start interaction solve at a=0, b=1"
```

---

## Task 2: `fit_model()` fast path when every parameter is fixed

The staged `reference` fit holds **all** curve params fixed and has no interaction params, so `free` is empty. `stats::optim` cannot run on a zero-length parameter vector — add a guard that just evaluates the model.

**Files:**
- Modify: `R/fit.R` (insert after `objective_of` is defined, before `theta0 <- par[free]`)
- Test: `tests/testthat/test-fit.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-fit.R`:

```r
test_that("fit_model with all parameters fixed evaluates without optimising", {
  df <- make_binary_df()  # truth: max=800, slope1=4, slope2=1.5, ec501=0.08, ec502=1
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   fixed = names(start))
  expect_equal(fit$df, 0)
  expect_equal(fit$convergence, 0)
  expect_equal(unname(fit$par[["max"]]), 800)     # untouched
  expect_lt(fit$objective, 1e-6)                  # truth params -> ~perfect fit
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run the test command for `test-fit.R`.
Expected: FAIL or ERROR in the new test (zero-length `optim` error / wrong `convergence`).

- [ ] **Step 3: Add the guard**

In `R/fit.R`, immediately **after** the `objective_of` function definition (the block ending `}` around the current line 64) and **before** `theta0 <- par[free]`, insert:

```r
  # When every parameter is fixed (the staged reference fit: curve params held at
  # their single-compound values, no interaction params to free) there is nothing
  # to optimise -- just evaluate the model. `stats::optim` cannot run on a
  # zero-length parameter vector, so short-circuit here.
  if (length(free) == 0) {
    pred <- predict_with(par)
    obs <- if (response == "continuous") df$Res else df$Affected / df$Exposed
    return(list(par = par, objective = objective_of(par), pred = pred,
                residuals = obs - pred, df = 0L, n = nrow(df),
                convergence = 0L, reference = reference, deviation = deviation,
                response = response, conc_cols = conc_cols, n_chem = n_chem,
                kind = "mixture"))
  }
```

- [ ] **Step 4: Run the test to verify it passes**

Run the test command for `test-fit.R`. Expected: all PASS (now includes the new test).

- [ ] **Step 5: Commit**

```bash
git add R/fit.R tests/testthat/test-fit.R
git commit -m "feat: fit_model evaluates directly when all params fixed"
```

---

## Task 3: `fit_curve_from_singles()` — stage-1 curve fit

Fit the `reference` model to control + single-compound rows only; return the named curve-parameter vector. Fall back to the `seed_from_singles()` heuristic if the subset cannot be fit.

**Files:**
- Modify: `R/analyse.R` (add helper after `seed_from_singles`, before `analyse_mixture`)
- Test: `tests/testthat/test-staged.R` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-staged.R`:

```r
# Synthetic binary data with a full single-compound dose series per chemical
# (>= 6 distinct marginal concentrations) so stage 1 can identify each curve,
# plus a few interior mixture points. Noiseless CA-reference surface.
make_staged_binary_df <- function(truth = list(max = 800, slope1 = 4,
                                               slope2 = 1.5, ec501 = 0.08,
                                               ec502 = 1)) {
  c1 <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
  c2 <- c(0.125, 0.25, 0.5, 1, 2, 4)
  singles <- rbind(data.frame(C1 = c(0, c1), C2 = 0),
                   data.frame(C1 = 0,        C2 = c2))
  mix <- expand.grid(C1 = c(0.04, 0.08), C2 = c(0.5, 1))
  d <- unique(rbind(singles, mix))
  d$Res <- ca_bi_vec(d$C1, d$C2, truth$max, truth$slope1, truth$slope2,
                     truth$ec501, truth$ec502)
  d
}

test_that("fit_curve_from_singles recovers curve params from the marginals", {
  d <- make_staged_binary_df()
  base <- fit_curve_from_singles(d, "CA", "continuous", n_starts = 1)
  expect_setequal(names(base), c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_equal(unname(base[["max"]]),   800,  tolerance = 1e-2)
  expect_equal(unname(base[["ec501"]]), 0.08, tolerance = 1e-2)
  expect_equal(unname(base[["ec502"]]), 1,    tolerance = 1e-2)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run the test command for `test-staged.R`.
Expected: FAIL — `could not find function "fit_curve_from_singles"`.

- [ ] **Step 3: Implement the helper**

In `R/analyse.R`, insert **after** the `seed_from_singles` function and **before** the `analyse_mixture` roxygen block:

```r
#' Stage-1 fit: estimate the curve parameters from single-compound data
#'
#' Fits the `reference` model to the control + single-compound rows only (rows
#' where at most one chemical is present). On those rows every mixture model
#' reduces exactly to the independent log-logistic curves, so this identifies the
#' shared `max` and per-chemical `slope`/`ec50` without any interaction
#' parameter. Falls back to the [seed_from_singles()] heuristic if the subset
#' cannot be fit (e.g. a chemical with too sparse a marginal series).
#' @keywords internal
fit_curve_from_singles <- function(df, reference, response,
                                   lower = NULL, upper = NULL,
                                   n_starts = 1, time_limit = 30) {
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  singles <- df[rowSums(df[cols] > 0) <= 1, , drop = FALSE]
  seed <- seed_from_singles(df, response)
  fit <- tryCatch(
    fit_model(singles, reference, "reference", response, start = seed,
              lower = lower, upper = upper, n_starts = n_starts,
              time_limit = time_limit),
    error = function(e) NULL)
  if (is.null(fit) || !all(is.finite(fit$par[names(seed)]))) return(seed)
  fit$par[names(seed)]
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run the test command for `test-staged.R`. Expected: PASS (1 test).

- [ ] **Step 5: Commit**

```bash
git add R/analyse.R tests/testthat/test-staged.R
git commit -m "feat: add stage-1 single-compound curve fit"
```

---

## Task 4: Wire staged fitting into `analyse_mixture()`

Replace the joint fit (seed → 4 full fits) with the two-stage flow: stage-1 curve fit, then 4 deviation fits with the curve params held fixed.

**Files:**
- Modify: `R/analyse.R` (`analyse_mixture`, currently lines ~72-98)
- Test: `tests/testthat/test-staged.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-staged.R`:

```r
test_that("analyse_mixture holds the curve params fixed across all four fits", {
  d <- make_staged_binary_df()
  res <- analyse_mixture(d, "CA", "continuous", n_starts = 1)
  base <- c("max", "slope1", "slope2", "ec501", "ec502")
  ref_base <- res$fits$reference$par[base]
  for (m in c("SA", "DR", "DL"))
    expect_equal(res$fits[[m]]$par[base], ref_base)
  # df now counts only the free interaction params.
  expect_equal(res$fits$reference$df, 0)
  expect_equal(res$fits$SA$df, 1)   # a
  expect_equal(res$fits$DR$df, 2)   # a, b
  expect_equal(res$fits$DL$df, 2)   # a, b
})

test_that("analyse_mixture chooses reference on noiseless reference data", {
  d <- make_staged_binary_df()
  res <- analyse_mixture(d, "CA", "continuous", n_starts = 1)
  expect_equal(res$chosen, "reference")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run the test command for `test-staged.R`.
Expected: FAIL — under the current joint fit the base params differ across fits (and `df` includes the 5 curve params, so `reference$df` is 5 not 0).

- [ ] **Step 3: Rewrite `analyse_mixture`**

In `R/analyse.R`, **replace** the body of `analyse_mixture` from `response <- match.arg(response)` down to the `fits <- lapply(...)` / `names(fits) <- devs` block with:

```r
  response <- match.arg(response)

  # Stage 1: estimate the curve parameters (max, slope*, ec50*) from the
  # single-compound data and hold them fixed. A user-supplied `start` is taken as
  # the fixed curve parameters directly.
  base <- if (is.null(start)) {
    fit_curve_from_singles(df, reference, response, lower, upper,
                           n_starts, time_limit)
  } else {
    start
  }
  base_names <- names(base)

  # Stage 2: with the curve parameters fixed, fit only the interaction
  # parameters (a, b, ...) to the full mixture data, once per deviation.
  devs <- c("reference", "SA", "DR", "DL")
  fits <- lapply(devs, function(d)
    fit_model(df, reference, d, response, start = base, fixed = base_names,
              lower = lower, upper = upper, n_starts = n_starts,
              time_limit = time_limit))
  names(fits) <- devs
```

Leave the `n <- nrow(df)` / `parent_of` / `comparison` / `list(...)` return block unchanged.

- [ ] **Step 4: Update the `analyse_mixture` roxygen**

In the roxygen block above `analyse_mixture`, replace the `@param start` line and add a one-line description of the staged behaviour. Replace:

```r
#' @param start Optional named starting vector; if NULL, seeded from single fits.
```

with:

```r
#' @param start Optional named vector of the curve parameters (`max`, `slope*`,
#'   `ec50*`) to hold fixed. When `NULL` (default) they are estimated from the
#'   single-compound data via [fit_curve_from_singles()].
```

And add, after the existing `@description`/title lines (before `@param df`), this sentence to the title block's details (append a new roxygen paragraph line):

```r
#' Fitting is staged: the curve parameters are fixed from the single-compound
#' data, then only the interaction parameters (`a`, `b`, ...) are fitted to the
#' mixture data. Each fit's `df` therefore counts only its free interaction
#' parameters; likelihood-ratio tests use differences, so they are unaffected.
```

- [ ] **Step 5: Run the tests to verify they pass**

Run the test command for `test-staged.R`. Expected: all PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add R/analyse.R tests/testthat/test-staged.R
git commit -m "feat: stage analyse_mixture into curve-fit then interaction-fit"
```

---

## Task 5: Re-baseline existing analyse_mixture tests + regenerate docs

The tiny synthetic grids in `test-nchem.R` give only one nonzero level per chemical, which the staged stage-1 cannot identify. Give them proper marginals. Then run the full suite and re-baseline the real-data validation expectations.

**Files:**
- Modify: `tests/testthat/test-nchem.R` (ternary end-to-end test)
- Verify (likely no change): `tests/testthat/test-summary.R`, `tests/testthat/test-analyse.R`, `tests/testthat/test-fit-bounds.R`
- Re-baseline if needed: `tests/testthat/test-validation-binary.R`
- Regenerate: `man/analyse_mixture.Rd`, `man/fit_model.Rd`

- [ ] **Step 1: Give the ternary end-to-end test proper marginals**

In `tests/testthat/test-nchem.R`, **replace** the second test (`"analyse_mixture runs end-to-end on a ternary dataset"`, currently lines 10-16) with:

```r
make_staged_ternary_df <- function() {
  # truth: max=800; slopes 4,1.5,1; ec50s 0.08,1,5  (matches ca_tri_vec arg order)
  c1 <- c(0.01, 0.02, 0.04, 0.08, 0.16)
  c2 <- c(0.125, 0.25, 0.5, 1, 2)
  c3 <- c(0.5, 1, 2, 4, 8)
  singles <- rbind(
    data.frame(C1 = c(0, c1), C2 = 0,        C3 = 0),
    data.frame(C1 = 0,        C2 = c(0, c2), C3 = 0),
    data.frame(C1 = 0,        C2 = 0,        C3 = c(0, c3)))
  mix <- expand.grid(C1 = c(0.04, 0.08), C2 = c(0.5, 1), C3 = c(1, 2))
  d <- unique(rbind(singles, mix))
  d$Res <- ca_tri_vec(d$C1, d$C2, d$C3, 800, 4, 1.5, 1, 0.08, 1, 5)
  d
}

test_that("analyse_mixture runs end-to-end on a ternary dataset", {
  g <- make_staged_ternary_df()
  res <- analyse_mixture(g, reference = "CA", response = "continuous")
  expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
  expect_equal(res$chosen, "reference")  # noiseless reference data
})
```

Leave the first test (`"model_spec handles ternary CA selections"`) and the third (`"single-chemical analysis ..."`) unchanged.

- [ ] **Step 2: Run the affected synthetic suites**

Run the test command for each of: `test-nchem.R`, `test-summary.R`, `test-analyse.R`, `test-fit-bounds.R`.
Expected: all PASS.

If `test-summary.R`'s `result_table` test errors (its tiny 2-level grid fails stage 1), replace its data line `g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))` + the `g$Res <- ...` line in **that test only** with `g <- make_staged_binary_df()` — but note `make_staged_binary_df` lives in `test-staged.R`; if needed, inline its body. Only do this if the test actually fails; the structure-only assertions usually still pass.

- [ ] **Step 3: Re-baseline the binary validation test**

Run the test command for `test-validation-binary.R`.

- If it PASSES: do nothing — staged fitting matched the workbook within tolerance. Skip to Step 4.
- If an **objective** expectation fails (e.g. `res$fits$reference$objective`), read the actual staged value from the failure message and update the constant in `tests/testthat/test-validation-binary.R`, keeping the existing `tolerance`. Add a short comment noting it was re-baselined for staged fitting.
- **STOP and report to the user** if a *conclusion* changes — i.e. `res$chosen` is no longer `"DL"` (continuous) / `"reference"` (quantal), or the CA-vs-SA test flips significance. Do **not** edit the assertion to pass; this is a scientific result the user must see. Capture the new `res$chosen`, `res$comparison`, and the four objectives and surface them.

- [ ] **Step 4: Regenerate man pages**

Run:

```
Rscript -e ".libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths())); devtools::document()"
```

Expected: `man/analyse_mixture.Rd` (and possibly `man/fit_model.Rd`) updated; `init_start_par` and `fit_curve_from_singles` produce **no** `.Rd` (they are `@keywords internal` and unexported). Confirm `NAMESPACE` is unchanged (no new exports).

- [ ] **Step 5: Run the entire test suite**

```
Rscript -e ".libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths())); suppressMessages(devtools::load_all('.')); testthat::test_dir('tests/testthat', reporter='summary')"
```

Expected: 0 failures. (`test-validation-binary` may `skip` if its fixture CSVs are absent — that is fine.)

- [ ] **Step 6: Commit**

```bash
git add tests/ man/
git commit -m "test: re-baseline analyse_mixture tests for staged fitting"
```

---

## Self-Review notes (already incorporated)

- **Spec coverage:** stage-1 single fit (Task 3), fix-then-fit-interaction (Task 4), `a=0`/`b=1` starts (Task 1), zero-free reference fit (Task 2), `df` semantics + docs (Task 4), test updates + re-baseline (Task 5). Fallback path (Task 3). Binary + ternary both covered (Tasks 4 & 5).
- **Out of scope (unchanged):** no `fit='staged'|'joint'` flag; ternary DR/DL model selection untouched; no new exports.
- **Type consistency:** helper names `init_start_par(params, extra, start)` and `fit_curve_from_singles(df, reference, response, lower, upper, n_starts, time_limit)` are used identically in their defining and calling tasks. `base_names` is the fixed-parameter set in both stage-1 return and stage-2 `fixed=`.
