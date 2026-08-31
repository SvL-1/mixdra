# mixdra Engine Implementation Plan

**Goal:** Build the headless computational engine of the `mixdra` R package — single-chemical and binary/ternary mixture dose–response fitting (CA/IA references + S/A/DR/DL deviations), maximum-likelihood objectives for continuous and binary data, and likelihood-ratio model comparison — fully validated against the existing Excel workbooks.

**Architecture:** A pure-R package (no Shiny). Scalar model functions are ported from the already-validated `MixTox_shiny_v2/functions/model_functions.R`, then vectorised. A swappable objective (sum-of-squares for continuous, binomial deviance for binary) feeds a staged optimiser built on base `stats::optim` (L-BFGS-B with a Nelder-Mead fallback and multi-start). A comparison layer runs nested likelihood-ratio tests reproducing the workbook's Table-2 output.

**Tech Stack:** R (≥4.1), base `stats` (`optim`), `numDeriv` (Hessian/CIs), `readxl` (fixture extraction), `testthat` (≥3e), `usethis`/`devtools` for scaffolding. (Note: the spec mentioned `bbmle`/`DEoptim`; we use base `optim` instead — `bbmle` is only a thin wrapper over `optim`, and avoiding both keeps dependencies minimal per the spec's lean-dependency goal. They can be reintroduced as fallbacks later if a fit proves stubborn.)

**Reference material (read before starting):**
- Spec: `docs/design/specs/2026-05-30-mixdra-design.md`
- Methodology: Jonker et al. 2005 (Eq. 14 single curve; Eq. 2/3 CA; Eq. 4/5 IA; Eq. 7 S/A; Eq. 8/11 DR; Eq. 12/13 DL; Eq. 16 binary objective).
- Existing validated model math: `MixTox_shiny_v2/functions/model_functions.R`.
- Validation workbook: `MixTox Model_binary_MPs_CPF.xlsm` (sheets listed below).

**Known-answer values (from `MixTox Model_binary_MPs_CPF.xlsm`, cached cell values):**

*Continuous data, CA reference (sheet "CA - continous data model (new)"):*
| Model | max | β1 | β2 | EC50₁ | EC50₂ | a | b | residual SS |
|---|---|---|---|---|---|---|---|---|
| CA | 796.0131 | 6.2505 | 0.3771 | 0.0820 | 100 | – | – | 1,633,769.68 |
| S/A | 797.4987 | 7.5841 | 0.3676 | 0.0758 | 100 | 9.9846 | – | 1,491,922.28 |
| DR | 797.101 | 7.1898 | 0.3723 | 0.0787 | 100 | 629.0985 | -635.663 | 1,483,249.47 |
| DL | 806.9692 | 8.1572 | 0.2105 | 0.0776 | 100 | 529.7646 | 0.7472 | 1,259,961.31 |

- N = 145. Continuous LR statistic: `chi = N * ln(SS_reduced / SS_full)`; df = extra params.
- CA vs S/A: chi=13.17, df=1, p=0.0003. CA vs DR: chi=14.01, df=2, p=0.0009. CA vs DL: chi=37.67, df=2, p≈0.

*Binary data, CA reference (sheet "CA - binary data model (new)"):*
| Model | max | β1 | β2 | EC50₁ | EC50₂ | a | b | residual deviance |
|---|---|---|---|---|---|---|---|---|
| CA | 0.9556 | 0.8805 | 2.4977 | 100 | 0.9168 | – | – | 184.424 |
| S/A | 0.9566 | 0.8248 | 2.4109 | 100 | 0.8902 | 3.7074 | – | 183.5394 |
| DR | 0.9562 | 0.7922 | 2.5060 | 100 | 0.9255 | -5.5163 | 167.2119 | 182.2878 |
| DL | 0.9565 | 0.8218 | 2.4098 | 100 | 0.8926 | 4.3005 | 0.2820 | 183.5356 |

- Binary LR statistic = difference in residual deviance; df = extra params. CA vs S/A: 0.8846, df=1, p=0.347.
- Note: several EC50 values rail to the bound 100 (the less-toxic chemical). The "scientifically equivalent" standard (spec §2) covers this — exact reproduction is **not** required.

**Conventions used throughout this plan:**
- Endpoint direction handled by slope sign: decreasing endpoints (survival, reproduction) use **slope > 0** and the control equals `max`; increasing endpoints use slope < 0. v1 assumes decreasing endpoints (matches the workbooks); sign is data-driven, not a separate code path.
- Mixture data frame columns (internal): `C1`, `C2`, `C3` (concentrations; absent columns are 0), and either `Res` (continuous response) or `Exposed` + `Affected` (binary counts). Observed proportion `pi = Affected / Exposed`.
- Run tests with: `R -q -e 'devtools::test()'` (or `testthat::test_local()`). Each task lists the exact file filter.

---

### Task 0: Package scaffold

**Files:**
- Create: `DESCRIPTION`
- Create: `NAMESPACE`
- Create: `mixdra.Rproj`
- Create: `tests/testthat.R`
- Create: `R/mixdra-package.R`

- [x] **Step 1: Create the package skeleton**

Create `DESCRIPTION`:
```
Package: mixdra
Title: Mixture Dose-Response Analysis (Jonker et al. 2005 Methodology)
Version: 0.0.0.9000
Authors@R:
    person("Mixdra", "Authors", email = "pieter.vklaveren@gmail.com", role = c("aut", "cre"))
Description: Fits single-chemical and binary/ternary mixture dose-response
    models using the concentration addition and independent action reference
    models with synergism/antagonism, dose-ratio-dependent, and
    dose-level-dependent deviation patterns, following Jonker et al. (2005)
    <doi:10.1897/04-431R.1>. Supports continuous and binary (quantal) responses
    and likelihood-ratio model comparison.
License: GPL-3
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.1
Depends: R (>= 4.1)
Imports:
    numDeriv,
    stats
Suggests:
    readxl,
    testthat (>= 3.0.0)
Config/testthat/edition: 3
```

Create `NAMESPACE`:
```
# Generated by roxygen2: do not edit by hand
```

Create `R/mixdra-package.R`:
```r
#' @keywords internal
"_PACKAGE"
```

Create `tests/testthat.R`:
```r
library(testthat)
library(mixdra)

test_check("mixdra")
```

Create `mixdra.Rproj`:
```
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: Default

BuildType: Package
PackageUseDevtools: Yes
PackageInstallArgs: --no-multiarch --with-keep.source
```

- [x] **Step 2: Verify the package loads**

Run: `R -q -e 'devtools::load_all(); cat("LOADED OK\n")'`
Expected: prints `LOADED OK` with no error.

- [x] **Step 3: Commit**

```bash
git add DESCRIPTION NAMESPACE mixdra.Rproj tests/testthat.R R/mixdra-package.R
git commit -m "chore: scaffold mixdra R package"
```

---

### Task 1: Single-chemical log-logistic model and fit

**Files:**
- Create: `R/single.R`
- Test: `tests/testthat/test-single.R`

- [x] **Step 1: Write the failing test**

`tests/testthat/test-single.R`:
```r
test_that("ll3_predict matches the log-logistic formula", {
  # Y = max / (1 + (C/EC50)^slope); at C = EC50, Y = max/2
  expect_equal(ll3_predict(c(0, 1, Inf), max = 100, slope = 2, ec50 = 1),
               c(100, 50, 0))
})

test_that("fit_single recovers known parameters from clean data", {
  conc <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
  truth <- list(max = 800, slope = 2, ec50 = 0.5)
  resp <- ll3_predict(conc, truth$max, truth$slope, truth$ec50)
  fit <- fit_single(conc, resp)
  expect_equal(unname(fit$par["max"]),  800, tolerance = 1e-3)
  expect_equal(unname(fit$par["slope"]), 2,   tolerance = 1e-3)
  expect_equal(unname(fit$par["ec50"]),  0.5, tolerance = 1e-3)
  expect_lt(fit$ssr, 1e-6)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-single.R")'`
Expected: FAIL — `could not find function "ll3_predict"`.

- [x] **Step 3: Write minimal implementation**

`R/single.R`:
```r
#' Three-parameter log-logistic dose-response prediction
#'
#' `Y = max / (1 + (C / ec50)^slope)`. With `slope > 0` the response
#' decreases with dose (control = `max`).
#' @param conc Numeric vector of concentrations (>= 0).
#' @param max,slope,ec50 Curve parameters.
#' @return Numeric vector of predicted responses.
#' @export
ll3_predict <- function(conc, max, slope, ec50) {
  max / (1 + (conc / ec50)^slope)
}

#' Fit a three-parameter log-logistic curve to single-chemical data
#'
#' Minimises the residual sum of squares.
#' @param conc Numeric vector of concentrations.
#' @param resp Numeric vector of responses (same length as `conc`).
#' @return A list with `par` (named: max, slope, ec50), `ssr`, and `convergence`.
#' @export
fit_single <- function(conc, resp) {
  stopifnot(length(conc) == length(resp), length(conc) > 3)
  pos <- conc[conc > 0]
  start <- c(max = max(resp, na.rm = TRUE),
             slope = 1,
             ec50 = stats::median(pos))
  obj <- function(p) {
    pred <- ll3_predict(conc, p[["max"]], p[["slope"]], p[["ec50"]])
    sum((resp - pred)^2)
  }
  lower <- c(max = 1e-8, slope = 1e-3, ec50 = 1e-8)
  upper <- c(max = Inf,  slope = 50,   ec50 = Inf)
  res <- stats::optim(start, obj, method = "L-BFGS-B",
                      lower = lower, upper = upper)
  list(par = res$par, ssr = res$value, convergence = res$convergence)
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-single.R")'`
Expected: PASS (2 tests).

- [x] **Step 5: Commit**

```bash
git add R/single.R tests/testthat/test-single.R
git commit -m "feat: single-chemical log-logistic model and fit"
```

---

### Task 2: Maximum-likelihood objectives (continuous + binary)

**Files:**
- Create: `R/objective.R`
- Test: `tests/testthat/test-objective.R`

- [x] **Step 1: Write the failing test**

`tests/testthat/test-objective.R`:
```r
test_that("obj_ss returns the sum of squared residuals", {
  expect_equal(obj_ss(c(1, 2, 3), c(1, 2, 3)), 0)
  expect_equal(obj_ss(c(1, 2, 3), c(2, 2, 2)), 1 + 0 + 1)
})

test_that("binlik per-row contribution matches the Jonker/VBA formula", {
  # P*log(pi_hat/pi) + (T-P)*log((1-pi_hat)/(1-pi)); pi = P/T
  # T=10, P=8 (pi=0.8), pi_hat=0.9656 -> matches workbook row (Contr.L = -2.0154)
  expect_equal(binlik(exposed = 10, affected = 8, pi_hat = 0.9656),
               8 * log(0.9656 / 0.8) + 2 * log((1 - 0.9656) / (1 - 0.8)),
               tolerance = 1e-9)
  expect_equal(round(binlik(10, 8, 0.9656), 4), -2.0154)
  # Saturated cell (pi = 1) contributes only the positive term
  expect_equal(binlik(10, 10, 0.9656), 10 * log(0.9656 / 1), tolerance = 1e-9)
})

test_that("obj_deviance is -2 * sum of per-row binlik contributions", {
  exposed <- c(10, 10, 10)
  affected <- c(10, 8, 10)
  pi_hat <- c(0.9656, 0.9656, 0.9656)
  expected <- -2 * sum(mapply(binlik, exposed, affected, pi_hat))
  expect_equal(obj_deviance(exposed, affected, pi_hat), expected)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-objective.R")'`
Expected: FAIL — `could not find function "obj_ss"`.

- [x] **Step 3: Write minimal implementation**

`R/objective.R`:
```r
#' Sum-of-squares objective for continuous responses
#' @param obs,pred Numeric vectors of observed and predicted responses.
#' @return The residual sum of squares.
#' @export
obj_ss <- function(obs, pred) {
  sum((obs - pred)^2)
}

#' Per-row binomial log-likelihood contribution (Jonker Eq. 16 / VBA BinLik)
#'
#' Returns `P*log(pi_hat/pi) + (T-P)*log((1-pi_hat)/(1-pi))`, where `pi = P/T`.
#' Saturated terms (pi == 0 or pi == 1) drop the undefined component.
#' @param exposed Number exposed (T).
#' @param affected Number responding (P).
#' @param pi_hat Model-predicted probability in (0, 1).
#' @return A scalar (<= 0) log-likelihood contribution.
#' @export
binlik <- function(exposed, affected, pi_hat) {
  pi <- affected / exposed
  pos <- if (pi > 0) affected * log(pi_hat / pi) else 0
  neg <- if (pi < 1) (exposed - affected) * log((1 - pi_hat) / (1 - pi)) else 0
  pos + neg
}

#' Binomial deviance objective for binary responses
#'
#' `-2 * sum(binlik)`; minimising this maximises the binomial likelihood.
#' @param exposed,affected Numeric vectors of counts.
#' @param pi_hat Numeric vector of predicted probabilities.
#' @return The residual deviance (scalar).
#' @export
obj_deviance <- function(exposed, affected, pi_hat) {
  -2 * sum(mapply(binlik, exposed, affected, pi_hat))
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-objective.R")'`
Expected: PASS (3 tests).

- [x] **Step 5: Commit**

```bash
git add R/objective.R tests/testthat/test-objective.R
git commit -m "feat: continuous and binary ML objectives"
```

---

### Task 3: Port binary mixture model functions

**Files:**
- Create: `R/models-binary.R`
- Test: `tests/testthat/test-models-binary.R`

**Source of truth:** Port the scalar functions verbatim (logic unchanged) from
`MixTox_shiny_v2/functions/model_functions.R`:
- `CA_bi` (lines ~171-227) → `ca_bi`
- `IA_bi` (lines ~1258-1290) → `ia_bi`
- `CA_SA_bi` (lines ~433-488) → `ca_sa_bi`
- `CA_DR_bi` (lines ~870-925) → `ca_dr_bi`
- `CA_DL_bi` (lines ~1117-1171) → `ca_dl_bi`
- `IA_SA_bi` (lines ~1431-1473) → `ia_sa_bi`
- `IA_DR_bi` (lines ~1615-1657) → `ia_dr_bi`
- `IA_DL_bi` (lines ~1817-1862) → `ia_dl_bi`

Rename arguments to a consistent lower-case signature:
`(c1, c2, max, slope1, slope2, ec501, ec502, ...)`. Keep `qnorm`/`pnorm` for the IA
deviation models (they replace Excel's `NormSInv`/`NormSDist`). At the end of the file add
vectorised wrappers so each accepts vector `c1`/`c2`:

```r
# Example wrapper (repeat the pattern for every model function):
ca_bi_vec     <- Vectorize(ca_bi,     vectorize.args = c("c1", "c2"))
ia_bi_vec     <- Vectorize(ia_bi,     vectorize.args = c("c1", "c2"))
ca_sa_bi_vec  <- Vectorize(ca_sa_bi,  vectorize.args = c("c1", "c2"))
ca_dr_bi_vec  <- Vectorize(ca_dr_bi,  vectorize.args = c("c1", "c2"))
ca_dl_bi_vec  <- Vectorize(ca_dl_bi,  vectorize.args = c("c1", "c2"))
ia_sa_bi_vec  <- Vectorize(ia_sa_bi,  vectorize.args = c("c1", "c2"))
ia_dr_bi_vec  <- Vectorize(ia_dr_bi,  vectorize.args = c("c1", "c2"))
ia_dl_bi_vec  <- Vectorize(ia_dl_bi,  vectorize.args = c("c1", "c2"))
```

- [x] **Step 1: Write the failing test**

`tests/testthat/test-models-binary.R`:
```r
test_that("ca_bi reduces to the single-chemical curve when one conc is 0", {
  # C2 = 0 -> Y = max / (1 + (C1/ec501)^slope1)
  expect_equal(ca_bi(c1 = 0.05, c2 = 0, max = 800, slope1 = 6, slope2 = 0.4,
                     ec501 = 0.08, ec502 = 100),
               800 / (1 + (0.05 / 0.08)^6))
  # Both zero, decreasing endpoint (positive slopes) -> control = max
  expect_equal(ca_bi(0, 0, 800, 6, 0.4, 0.08, 100), 800)
})

test_that("ia_bi matches the closed-form independent-action product", {
  f1 <- 1 / (1 + (0.05 / 0.08)^6)
  f2 <- 1 / (1 + (10 / 50)^0.4)
  expect_equal(ia_bi(0.05, 10, 800, 6, 0.4, 0.08, 50), 800 * f1 * f2)
})

test_that("ca_sa_bi reduces to ca_bi when a = 0", {
  expect_equal(ca_sa_bi(0.05, 10, 800, 6, 0.4, 0.08, 50, a = 0),
               ca_bi(0.05, 10, 800, 6, 0.4, 0.08, 50),
               tolerance = 1e-5)
})

test_that("vectorised wrapper maps over concentration vectors", {
  out <- ca_bi_vec(c1 = c(0, 0.05), c2 = c(0, 0),
                   max = 800, slope1 = 6, slope2 = 0.4, ec501 = 0.08, ec502 = 100)
  expect_length(out, 2)
  expect_equal(out[1], 800)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-models-binary.R")'`
Expected: FAIL — `could not find function "ca_bi"`.

- [x] **Step 3: Write the implementation**

Create `R/models-binary.R` by porting the eight scalar functions named above from
`MixTox_shiny_v2/functions/model_functions.R` (logic unchanged, lower-case argument names),
then append the eight `Vectorize` wrappers shown above. Each function must be preceded by a
short roxygen `#'` block and tagged `@keywords internal` (these are engine internals, not
user API).

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-models-binary.R")'`
Expected: PASS (4 tests).

- [x] **Step 5: Commit**

```bash
git add R/models-binary.R tests/testthat/test-models-binary.R
git commit -m "feat: port binary mixture model functions (CA/IA + S/A/DR/DL)"
```

---

### Task 4: Model registry (parameters per model)

**Files:**
- Create: `R/registry.R`
- Test: `tests/testthat/test-registry.R`

Maps a `(reference, deviation, n_chem)` selection to its predictor function, free-parameter
names, and nesting parent. This drives both the fitter (Task 5) and the comparison (Task 6).

- [x] **Step 1: Write the failing test**

`tests/testthat/test-registry.R`:
```r
test_that("binary CA registry entries expose the right parameters", {
  ref <- model_spec("CA", "reference", n_chem = 2)
  expect_equal(ref$params, c("max", "slope1", "slope2", "ec501", "ec502"))
  expect_null(ref$parent)
  expect_identical(ref$fn, ca_bi_vec)

  sa <- model_spec("CA", "SA", n_chem = 2)
  expect_equal(sa$params, c("max", "slope1", "slope2", "ec501", "ec502", "a"))
  expect_equal(sa$parent, "reference")
  expect_equal(sa$extra, "a")

  dr <- model_spec("CA", "DR", n_chem = 2)
  expect_equal(dr$extra, c("a", "b"))
  expect_equal(dr$parent, "SA")

  dl <- model_spec("CA", "DL", n_chem = 2)
  expect_equal(dl$extra, c("a", "b"))
  expect_equal(dl$parent, "SA")
})

test_that("binary IA registry resolves to the IA predictor functions", {
  expect_identical(model_spec("IA", "reference", 2)$fn, ia_bi_vec)
  expect_identical(model_spec("IA", "DR", 2)$fn, ia_dr_bi_vec)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-registry.R")'`
Expected: FAIL — `could not find function "model_spec"`.

- [x] **Step 3: Write minimal implementation**

`R/registry.R`:
```r
#' Look up the specification for a mixture model
#'
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param n_chem Number of chemicals (2 for binary; 3 handled in a later task).
#' @return A list: `fn` (vectorised predictor), `params` (all free parameter
#'   names), `extra` (deviation parameters beyond the reference), `parent`
#'   (the deviation this one nests within, or NULL for the reference).
#' @keywords internal
model_spec <- function(reference, deviation, n_chem) {
  reference <- match.arg(reference, c("CA", "IA"))
  deviation <- match.arg(deviation, c("reference", "SA", "DR", "DL"))
  if (n_chem != 2) stop("model_spec: only n_chem = 2 implemented in this task")

  base_params <- c("max", "slope1", "slope2", "ec501", "ec502")
  extra <- switch(deviation,
                  reference = character(0),
                  SA = "a",
                  DR = c("a", "b"),
                  DL = c("a", "b"))
  parent <- switch(deviation,
                   reference = NULL,
                   SA = "reference",
                   DR = "SA",
                   DL = "SA")
  key <- paste0(tolower(reference), "_",
                switch(deviation, reference = "bi", SA = "sa_bi",
                       DR = "dr_bi", DL = "dl_bi"), "_vec")
  list(fn = get(key, mode = "function"),
       params = c(base_params, extra),
       extra = extra,
       parent = parent)
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-registry.R")'`
Expected: PASS (2 tests).

- [x] **Step 5: Commit**

```bash
git add R/registry.R tests/testthat/test-registry.R
git commit -m "feat: model registry mapping selections to predictors and params"
```

---

### Task 5: Staged mixture fitter

**Files:**
- Create: `R/fit.R`
- Test: `tests/testthat/test-fit.R`

Fits one model (e.g. CA + S/A) to a mixture dataset by minimising the appropriate objective.
Staged seeding: start values come from single-chemical fits (caller-supplied); deviation
parameters start at 0. Supports continuous and binary data and fixing parameters.

- [x] **Step 1: Write the failing test**

`tests/testthat/test-fit.R`:
```r
make_binary_df <- function() {
  # Tiny synthetic 2-chem dataset, decreasing endpoint
  expand <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  truth <- list(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expand$Res <- ca_bi_vec(expand$C1, expand$C2, truth$max, truth$slope1,
                          truth$slope2, truth$ec501, truth$ec502)
  expand
}

test_that("fit_model (continuous CA) recovers parameters and reports SS objective", {
  df <- make_binary_df()
  start <- c(max = 700, slope1 = 2, slope2 = 1, ec501 = 0.1, ec502 = 2)
  fit <- fit_model(df, reference = "CA", deviation = "reference",
                   response = "continuous", start = start)
  expect_equal(unname(fit$par["ec501"]), 0.08, tolerance = 1e-2)
  expect_lt(fit$objective, 1e-2)          # near-perfect fit on noiseless data
  expect_equal(fit$df, 5)                  # 5 free parameters
})

test_that("fit_model respects fixed parameters", {
  df <- make_binary_df()
  start <- c(max = 800, slope1 = 2, slope2 = 1, ec501 = 0.08, ec502 = 1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = start,
                   fixed = c("max", "ec501", "ec502"))
  expect_equal(unname(fit$par["max"]), 800)    # untouched
  expect_equal(unname(fit$par["ec501"]), 0.08) # untouched
  expect_equal(fit$df, 2)                       # only slope1, slope2 free
})

test_that("fit_model (binary) minimises deviance and predicts probabilities", {
  df <- make_binary_df()
  df$Exposed <- 10
  # Build affected counts from a known probability surface (max as proportion)
  p <- ca_bi_vec(df$C1, df$C2, 0.95, 4, 1.5, 0.08, 1)
  df$Affected <- round(p * df$Exposed)
  start <- c(max = 0.9, slope1 = 2, slope2 = 1, ec501 = 0.1, ec502 = 2)
  fit <- fit_model(df, "CA", "reference", response = "binary", start = start)
  expect_true(is.finite(fit$objective))
  expect_true(all(fit$pred >= 0 & fit$pred <= 1))
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-fit.R")'`
Expected: FAIL — `could not find function "fit_model"`.

- [x] **Step 3: Write minimal implementation**

`R/fit.R`:
```r
#' Fit one mixture model to a dataset
#'
#' @param df Data frame with `C1`, `C2` (and `C3` for ternary) plus either
#'   `Res` (continuous) or `Exposed`/`Affected` (binary).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param start Named numeric vector of starting values for the base parameters
#'   (max, slope1, slope2, ec501, ec502). Deviation parameters default to 0.
#' @param fixed Character vector of parameter names to hold fixed at `start`.
#' @param lower_frac,upper_mult Box-constraint multipliers applied to positive
#'   base parameters (deviation params are unconstrained).
#' @return A list with `par` (named fitted parameters), `objective`
#'   (residual SS or deviance), `pred` (fitted values), `residuals`, `df`
#'   (number of free parameters), `n`, and `convergence`.
#' @export
fit_model <- function(df, reference, deviation = "reference",
                      response = c("continuous", "binary"),
                      start, fixed = character(0),
                      lower_frac = 0.1, upper_mult = 10) {
  response <- match.arg(response)
  n_chem <- sum(c("C1", "C2", "C3") %in% names(df))
  spec <- model_spec(reference, deviation, n_chem)

  # Assemble the full starting vector. Zero-initialising means any deviation
  # parameter (a, b, b1, b2, b3) not supplied in `start` defaults to 0.
  par <- stats::setNames(numeric(length(spec$params)), spec$params)
  par[names(start)] <- start[names(start)]

  free <- setdiff(spec$params, fixed)
  conc <- list(c1 = df$C1, c2 = df$C2)

  predict_with <- function(p_full) {
    args <- c(conc, as.list(p_full))
    names(args) <- c("c1", "c2", names(p_full))
    do.call(spec$fn, args)
  }
  objective_of <- function(p_full) {
    pred <- predict_with(p_full)
    if (response == "continuous") {
      obj_ss(df$Res, pred)
    } else {
      obj_deviance(df$Exposed, df$Affected, pmin(pmax(pred, 1e-8), 1 - 1e-8))
    }
  }

  obj_free <- function(theta) {
    p_full <- par
    p_full[free] <- theta
    val <- objective_of(p_full)
    if (!is.finite(val)) 1e12 else val
  }

  theta0 <- par[free]
  base <- setdiff(free, spec$extra)   # bounded curve params; deviation params are unconstrained
  lower <- stats::setNames(rep(-Inf, length(free)), free)
  upper <- stats::setNames(rep(Inf, length(free)), free)
  lower[base] <- pmax(1e-8, theta0[base] * lower_frac)
  upper[base] <- pmax(lower[base] * 1.01, theta0[base] * upper_mult)

  res <- stats::optim(theta0, obj_free, method = "L-BFGS-B",
                      lower = lower[free], upper = upper[free])

  par[free] <- res$par
  pred <- predict_with(par)
  obs <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  list(par = par, objective = res$value, pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = res$convergence)
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-fit.R")'`
Expected: PASS (3 tests).

- [x] **Step 5: Commit**

```bash
git add R/fit.R tests/testthat/test-fit.R
git commit -m "feat: staged mixture fitter (continuous + binary, fixed params)"
```

---

### Task 6: Likelihood-ratio model comparison

**Files:**
- Create: `R/compare.R`
- Test: `tests/testthat/test-compare.R`

Computes the LR statistic and p-value between a model and its nesting parent, and selects the
most parsimonious model. Continuous: `chi = n * ln(SS_parent / SS_child)`. Binary:
`chi = deviance_parent - deviance_child`. df = difference in free parameters.

- [x] **Step 1: Write the failing test**

`tests/testthat/test-compare.R`:
```r
test_that("lr_test (continuous) matches the workbook CA-vs-S/A statistic", {
  # From MixTox Model_binary_MPs_CPF.xlsm, CA continuous: N=145
  out <- lr_test(obj_parent = 1633769.68, obj_child = 1491922.28,
                 df_parent = 5, df_child = 6, n = 145, response = "continuous")
  expect_equal(round(out$chi, 2), 13.17)
  expect_equal(out$df, 1)
  expect_equal(round(out$p, 4), 0.0003)
})

test_that("lr_test (binary) is the deviance difference", {
  # CA vs S/A binary: 184.424 - 183.5394 = 0.8846, df = 1
  out <- lr_test(obj_parent = 184.424, obj_child = 183.5394,
                 df_parent = 5, df_child = 6, n = 145, response = "binary")
  expect_equal(round(out$chi, 4), 0.8846)
  expect_equal(round(out$p, 3), 0.347)
})

test_that("select_parsimonious keeps the simplest model not significantly beaten", {
  # Reference not improved by any deviation -> pick reference
  fits <- list(
    reference = list(objective = 100, df = 5),
    SA        = list(objective = 99.9, df = 6),
    DR        = list(objective = 99.8, df = 7),
    DL        = list(objective = 99.7, df = 7)
  )
  pick <- select_parsimonious(fits, n = 145, response = "continuous", alpha = 0.05)
  expect_equal(pick, "reference")
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-compare.R")'`
Expected: FAIL — `could not find function "lr_test"`.

- [x] **Step 3: Write minimal implementation**

`R/compare.R`:
```r
#' Likelihood-ratio test between a model and its nesting parent
#'
#' @param obj_parent,obj_child Objective values (SS for continuous, deviance for
#'   binary) of the simpler (parent) and more complex (child) models.
#' @param df_parent,df_child Number of free parameters in each model.
#' @param n Number of observations (used for the continuous statistic).
#' @param response "continuous" or "binary".
#' @return A list with `chi`, `df`, and `p` (upper-tail chi-squared p-value).
#' @export
lr_test <- function(obj_parent, obj_child, df_parent, df_child, n, response) {
  ddf <- df_child - df_parent
  chi <- if (response == "continuous") {
    n * log(obj_parent / obj_child)
  } else {
    obj_parent - obj_child
  }
  list(chi = chi, df = ddf, p = stats::pchisq(chi, df = ddf, lower.tail = FALSE))
}

#' Select the most parsimonious model from a set of nested fits
#'
#' Walks reference -> SA -> {DR, DL}; a more complex model is accepted only if
#' it significantly improves on its parent (LR test p < alpha).
#' @param fits Named list of fits (`reference`, `SA`, `DR`, `DL`), each with
#'   `objective` and `df`.
#' @param n Number of observations.
#' @param response "continuous" or "binary".
#' @param alpha Significance threshold.
#' @return The name of the selected model.
#' @export
select_parsimonious <- function(fits, n, response, alpha = 0.05) {
  improves <- function(child, parent) {
    t <- lr_test(fits[[parent]]$objective, fits[[child]]$objective,
                 fits[[parent]]$df, fits[[child]]$df, n, response)
    isTRUE(t$p < alpha)
  }
  best <- "reference"
  if (!is.null(fits$SA) && improves("SA", "reference")) {
    best <- "SA"
    dr_ok <- !is.null(fits$DR) && improves("DR", "SA")
    dl_ok <- !is.null(fits$DL) && improves("DL", "SA")
    if (dr_ok || dl_ok) {
      cand <- c(DR = if (dr_ok) fits$DR$objective else NA,
                DL = if (dl_ok) fits$DL$objective else NA)
      best <- names(which.min(cand))  # smaller objective = better fit
    }
  }
  best
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-compare.R")'`
Expected: PASS (3 tests).

- [x] **Step 5: Commit**

```bash
git add R/compare.R tests/testthat/test-compare.R
git commit -m "feat: likelihood-ratio comparison and parsimonious model selection"
```

---

### Task 7: Full analysis driver (`analyse_mixture`)

**Files:**
- Create: `R/analyse.R`
- Test: `tests/testthat/test-analyse.R`

Ties the pieces together: given a mixture dataset and a reference, fit single curves for the
seeds, fit {reference, SA, DR, DL}, run the comparisons, and assemble the Table-2 result
object (the future B10:P20 block).

- [x] **Step 1: Write the failing test**

`tests/testthat/test-analyse.R`:
```r
test_that("analyse_mixture returns a fit per deviation and a chosen model", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  res <- analyse_mixture(df, reference = "CA", response = "continuous")
  expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
  expect_true(res$chosen %in% c("reference", "SA", "DR", "DL"))
  # On noiseless reference-generated data, the reference should win.
  expect_equal(res$chosen, "reference")
  # Comparison table has one row per non-reference model.
  expect_equal(nrow(res$comparison), 3)
  expect_true(all(c("model", "chi", "df", "p") %in% names(res$comparison)))
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-analyse.R")'`
Expected: FAIL — `could not find function "analyse_mixture"`.

- [x] **Step 3: Write minimal implementation**

`R/analyse.R`:
```r
#' Seed starting values from per-chemical single-curve fits
#' @keywords internal
seed_from_singles <- function(df, response) {
  resp_col <- if (response == "continuous") df$Res else df$Affected / df$Exposed
  fit1 <- fit_single(df$C1[df$C2 == 0], resp_col[df$C2 == 0])
  fit2 <- fit_single(df$C2[df$C1 == 0], resp_col[df$C1 == 0])
  c(max    = mean(c(fit1$par[["max"]], fit2$par[["max"]])),
    slope1 = fit1$par[["slope"]], slope2 = fit2$par[["slope"]],
    ec501  = fit1$par[["ec50"]],  ec502  = fit2$par[["ec50"]])
}

#' Analyse a mixture: fit reference + deviations and compare
#'
#' @param df Mixture data frame (see [fit_model()]).
#' @param reference "CA" or "IA".
#' @param response "continuous" or "binary".
#' @param start Optional named starting vector; if NULL, seeded from single fits.
#' @param alpha Significance threshold for model selection.
#' @return A list: `fits` (named list of model fits), `comparison` (data frame of
#'   LR tests vs each model's parent), `chosen` (selected model name),
#'   `reference`, `response`.
#' @export
analyse_mixture <- function(df, reference, response = c("continuous", "binary"),
                            start = NULL, alpha = 0.05) {
  response <- match.arg(response)
  if (is.null(start)) start <- seed_from_singles(df, response)

  devs <- c("reference", "SA", "DR", "DL")
  fits <- lapply(devs, function(d)
    fit_model(df, reference, d, response, start = start))
  names(fits) <- devs

  n <- nrow(df)
  parent_of <- c(SA = "reference", DR = "SA", DL = "SA")
  comparison <- do.call(rbind, lapply(names(parent_of), function(m) {
    p <- parent_of[[m]]
    t <- lr_test(fits[[p]]$objective, fits[[m]]$objective,
                 fits[[p]]$df, fits[[m]]$df, n, response)
    data.frame(model = m, parent = p, chi = t$chi, df = t$df, p = t$p)
  }))

  list(fits = fits, comparison = comparison,
       chosen = select_parsimonious(fits, n, response, alpha),
       reference = reference, response = response)
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-analyse.R")'`
Expected: PASS (1 test).

- [x] **Step 5: Commit**

```bash
git add R/analyse.R tests/testthat/test-analyse.R
git commit -m "feat: analyse_mixture driver assembling fits and comparison"
```

---

### Task 8: Validation against the binary workbook fixture

**Files:**
- Create: `tests/testthat/fixtures/extract_binary_fixture.R` (one-off generator)
- Create: `tests/testthat/fixtures/binary_mps_cpf_continuous.csv` (generated)
- Test: `tests/testthat/test-validation-binary.R`

Proves the engine reproduces the published analysis to the agreed "scientifically equivalent"
standard: same chosen model and parameters in the right ballpark.

- [x] **Step 1: Generate the data fixture from the workbook**

Create `tests/testthat/fixtures/extract_binary_fixture.R`:
```r
# One-off: extract the continuous CA data columns (A=[CPF], B=[MPs], C=Data)
# from the binary workbook into a CSV checked into the repo.
# Requires the workbook present at the repo root and the `readxl` package.
library(readxl)
wb <- "MixTox Model_binary_MPs_CPF.xlsm"
raw <- read_excel(wb, sheet = "CA - continous data model (new)",
                  range = "A25:C322", col_names = c("C1", "C2", "Res"))
raw <- raw[stats::complete.cases(raw), ]
write.csv(raw, "tests/testthat/fixtures/binary_mps_cpf_continuous.csv",
          row.names = FALSE)
cat("wrote", nrow(raw), "continuous rows\n")

# Also extract the binary (quantal) data: A=[MPs], B=[Imi], C=Survivors, D=Exposed
bin <- read_excel(wb, sheet = "CA - binary data model (new)",
                  range = "A25:D322", col_names = c("C1", "C2", "Affected", "Exposed"))
bin <- bin[stats::complete.cases(bin), ]
write.csv(bin, "tests/testthat/fixtures/binary_mps_cpf_quantal.csv",
          row.names = FALSE)
cat("wrote", nrow(bin), "quantal rows\n")
```

Run: `R -q -e 'source("tests/testthat/fixtures/extract_binary_fixture.R")'`
Expected: prints `wrote 145 continuous rows` and `wrote 145 quantal rows`; both CSVs created.

- [x] **Step 2: Write the validation test**

`tests/testthat/test-validation-binary.R`:
```r
test_that("engine reproduces the binary workbook CA continuous analysis", {
  csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
  skip_if_not(file.exists(csv), "fixture CSV not generated")
  df <- read.csv(csv)

  res <- analyse_mixture(df, reference = "CA", response = "continuous")

  # Reference fit: residual SS within 2% of the workbook (1,633,769.68).
  expect_equal(res$fits$reference$objective, 1633769.68, tolerance = 0.02)
  # EC50 of the toxic chemical (CPF/C1) near the workbook value 0.082.
  expect_equal(unname(res$fits$reference$par["ec501"]), 0.082, tolerance = 0.05)
  # Workbook concludes DL is the most parsimonious (CA vs DL p ~ 0).
  expect_equal(res$chosen, "DL")
  # CA-vs-S/A LR statistic in the right ballpark (workbook 13.17).
  ca_sa <- res$comparison[res$comparison$model == "SA", ]
  expect_gt(ca_sa$chi, 8)
  expect_lt(ca_sa$p, 0.01)
})

test_that("engine reproduces the binary workbook CA quantal analysis", {
  csv <- testthat::test_path("fixtures", "binary_mps_cpf_quantal.csv")
  skip_if_not(file.exists(csv), "quantal fixture CSV not generated")
  df <- read.csv(csv)   # columns C1, C2, Affected, Exposed

  res <- analyse_mixture(df, reference = "CA", response = "binary")

  # Reference residual deviance within 2% of the workbook (184.424).
  expect_equal(res$fits$reference$objective, 184.424, tolerance = 0.02)
  # control proportion (max) near the workbook value 0.9556.
  expect_equal(unname(res$fits$reference$par["max"]), 0.9556, tolerance = 0.02)
  # Workbook: no deviation is significant (all p > 0.3) -> reference chosen.
  expect_equal(res$chosen, "reference")
})
```

- [x] **Step 3: Run the validation tests**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-validation-binary.R")'`
Expected: PASS (2 tests). If a chosen model or objective is off, tune the fitter (multi-start;
widen `upper_mult` for the railing EC50₂; try Nelder-Mead fallback) — see Task 9 — until
scientifically equivalent.

- [x] **Step 4: Commit**

```bash
git add tests/testthat/fixtures/extract_binary_fixture.R tests/testthat/fixtures/binary_mps_cpf_continuous.csv tests/testthat/fixtures/binary_mps_cpf_quantal.csv tests/testthat/test-validation-binary.R
git commit -m "test: validate engine against binary MPs+CPF workbook (continuous + quantal CA)"
```

---

### Task 9: Optimiser robustness — multi-start + fallback

**Files:**
- Modify: `R/fit.R` (extend `fit_model`)
- Modify: `R/analyse.R` (forward `n_starts` from `analyse_mixture` to `fit_model`)
- Test: `tests/testthat/test-fit-robust.R`

If Task 8 already passes cleanly with a single L-BFGS-B start, keep this minimal. Otherwise add
multi-start (perturbed seeds) and a Nelder-Mead fallback to escape local minima and handle the
non-smooth CA bisection surface.

- [x] **Step 1: Write the failing test**

`tests/testthat/test-fit-robust.R`:
```r
test_that("fit_model with n_starts returns the best of several starts", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  # A deliberately poor single start; multi-start should still find a good fit.
  bad <- c(max = 400, slope1 = 0.5, slope2 = 0.5, ec501 = 1, ec502 = 0.1)
  fit <- fit_model(df, "CA", "reference", "continuous", start = bad, n_starts = 8)
  expect_lt(fit$objective, 1)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-fit-robust.R")'`
Expected: FAIL — `unused argument (n_starts = 8)`.

- [x] **Step 3: Extend the implementation**

In `R/fit.R`, add an `n_starts = 1` parameter to `fit_model()`. Wrap the existing optimisation
in a loop: for start index `i`, use `theta0` for `i == 1` and `theta0 * runif(length(theta0), 0.5, 1.5)`
(clamped to `[lower, upper]`) for `i > 1`; run L-BFGS-B, and if it does not converge
(`res$convergence != 0`) retry that start with `method = "Nelder-Mead"`. Keep the fit with the
smallest finite objective. The return structure is unchanged.

Also in `R/analyse.R`: add an `n_starts = 5` argument to `analyse_mixture()` and pass it into
each `fit_model(...)` call so the workbook validations benefit from multi-start. The default of
5 keeps small datasets fast.

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-fit-robust.R")'`
Expected: PASS.

- [x] **Step 5: Re-run the workbook validations with multi-start**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-validation-binary.R")'`
Expected: PASS (2 tests) — confirm Task 8 still holds (now via `analyse_mixture`'s multi-start).

- [x] **Step 6: Commit**

```bash
git add R/fit.R R/analyse.R tests/testthat/test-fit-robust.R
git commit -m "feat: multi-start and Nelder-Mead fallback wired through analyse_mixture"
```

---

### Task 10: Port ternary model functions

**Files:**
- Create: `R/models-ternary.R`
- Test: `tests/testthat/test-models-ternary.R`

**Source of truth:** Port the scalar three-chemical functions from
`MixTox_shiny_v2/functions/model_functions.R` (logic unchanged, lower-case args
`c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, ...`):
- `CA` (lines ~13-165) → `ca_tri`
- `IA` (lines ~1175-1255) → `ia_tri`
- `CA_SA` (lines ~230-429) → `ca_sa_tri`
- `CA_DR` (lines ~676-867) → `ca_dr_tri`
- `CA_DL` (lines ~927-1112) → `ca_dl_tri`
- `IA_SA` (lines ~1293-1429) → `ia_sa_tri`
- `IA_DR` (lines ~1476-1612) → `ia_dr_tri`
- `IA_DL` (lines ~1660-1813) → `ia_dl_tri`

Append `Vectorize(..., vectorize.args = c("c1","c2","c3"))` wrappers named
`ca_tri_vec`, `ia_tri_vec`, `ca_sa_tri_vec`, etc.

- [x] **Step 1: Write the failing test**

`tests/testthat/test-models-ternary.R`:
```r
test_that("ca_tri reduces to a single-chemical curve when two concs are 0", {
  expect_equal(ca_tri(c1 = 0.05, c2 = 0, c3 = 0, max = 800,
                      slope1 = 6, slope2 = 0.4, slope3 = 1,
                      ec50_1 = 0.08, ec50_2 = 100, ec50_3 = 10),
               800 / (1 + (0.05 / 0.08)^6))
})

test_that("ca_tri reduces to ca_bi when the third conc is 0", {
  expect_equal(
    ca_tri(0.05, 10, 0, 800, 6, 0.4, 1, 0.08, 50, 10),
    ca_bi(0.05, 10, 800, 6, 0.4, 0.08, 50),
    tolerance = 1e-4)
})

test_that("ia_tri matches the triple independent-action product", {
  f1 <- 1 / (1 + (0.05 / 0.08)^6)
  f2 <- 1 / (1 + (10 / 50)^0.4)
  f3 <- 1 / (1 + (2 / 10)^1)
  expect_equal(ia_tri(0.05, 10, 2, 800, 6, 0.4, 1, 0.08, 50, 10),
               800 * f1 * f2 * f3)
})

test_that("vectorised ternary wrapper maps over three concentration vectors", {
  out <- ca_tri_vec(c(0, 0.05), c(0, 0), c(0, 0), 800, 6, 0.4, 1, 0.08, 100, 10)
  expect_length(out, 2)
  expect_equal(out[1], 800)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-models-ternary.R")'`
Expected: FAIL — `could not find function "ca_tri"`.

- [x] **Step 3: Write the implementation**

Port the eight ternary functions and add the vectorised wrappers as described above. Tag each
`@keywords internal`.

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-models-ternary.R")'`
Expected: PASS (4 tests).

- [x] **Step 5: Commit**

```bash
git add R/models-ternary.R tests/testthat/test-models-ternary.R
git commit -m "feat: port ternary mixture model functions"
```

---

### Task 11: Generalise registry, fitter, and analysis to 1 and 3 chemicals

**Files:**
- Modify: `R/registry.R`
- Modify: `R/fit.R`
- Modify: `R/analyse.R`
- Test: `tests/testthat/test-nchem.R`

- [x] **Step 1: Write the failing test**

`tests/testthat/test-nchem.R`:
```r
test_that("model_spec handles ternary CA selections", {
  ref <- model_spec("CA", "reference", n_chem = 3)
  expect_equal(ref$params,
               c("max", "slope1", "slope2", "slope3",
                 "ec50_1", "ec50_2", "ec50_3"))
  expect_identical(ref$fn, ca_tri_vec)
  expect_equal(model_spec("CA", "SA", 3)$extra, "a")
})

test_that("analyse_mixture runs end-to-end on a ternary dataset", {
  g <- expand.grid(C1 = c(0, 0.05), C2 = c(0, 0.5), C3 = c(0, 2))
  g$Res <- ca_tri_vec(g$C1, g$C2, g$C3, 800, 4, 1.5, 1, 0.08, 1, 5)
  res <- analyse_mixture(g, reference = "CA", response = "continuous")
  expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
  expect_equal(res$chosen, "reference")  # noiseless reference data
})

test_that("single-chemical analysis returns just the dose-response fit", {
  d <- data.frame(C1 = c(0, 0.01, 0.03, 0.1, 0.3, 1, 3),
                  Res = ll3_predict(c(0, 0.01, 0.03, 0.1, 0.3, 1, 3), 800, 2, 0.5))
  fit <- analyse_single(d)
  expect_equal(unname(fit$par["ec50"]), 0.5, tolerance = 1e-2)
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-nchem.R")'`
Expected: FAIL — ternary branch errors / `analyse_single` not found.

- [x] **Step 3: Extend the implementations**

In `R/registry.R`: extend `model_spec()` for `n_chem == 3`:
- base params `c("max","slope1","slope2","slope3","ec50_1","ec50_2","ec50_3")`;
- key suffix `_tri_vec` (e.g. `ca_dr_tri_vec`);
- `extra` per deviation — **note ternary DR differs from binary**, matching the ported
  function signatures from `model_functions.R`:
  - `SA` → `extra = "a"`, `parent = "reference"` (function arg list: `..., a`);
  - `DR` → `extra = c("a","b1","b2","b3")`, `parent = "SA"` (function arg list: `..., a, b1, b2, b3`);
  - `DL` → `extra = c("a","b")`, `parent = "SA"` (function arg list: `..., a, b`).
  (Verify each against the actual ported argument lists in `R/models-ternary.R` before wiring.)

In `R/fit.R`: build `conc` from whichever of `C1`,`C2`,`C3` are present, naming them
`c1`,`c2`,`c3`, and generalise `predict_with`'s name assignment to
`names(args) <- c(names(conc), names(p_full))`. No change to the `base`/bounds logic is needed
— Task 5 already derives `base <- setdiff(free, spec$extra)`, so the ternary curve parameters
(`slope3`, `ec50_*`) are bounded and the deviation parameters (`a`, `b`, `b1`, `b2`, `b3`)
remain unconstrained automatically.

In `R/analyse.R`: generalise `seed_from_singles()` to fit one single curve per present
chemical (chemical *i* seeded from rows where all other concentrations are 0) and name seeds
`slope{i}`/`ec50_{i}` to match the active registry. Add:
```r
#' Analyse a single chemical (dose-response curve only)
#' @param df Data frame with `C1` and `Res` (or `Exposed`/`Affected`).
#' @return The [fit_single()] result.
#' @export
analyse_single <- function(df) {
  resp <- if ("Res" %in% names(df)) df$Res else df$Affected / df$Exposed
  fit_single(df$C1, resp)
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-nchem.R")'`
Expected: PASS (3 tests).

- [x] **Step 5: Commit**

```bash
git add R/registry.R R/fit.R R/analyse.R tests/testthat/test-nchem.R
git commit -m "feat: generalise engine to 1 and 3 chemicals"
```

---

### Task 12: Confidence intervals and the Table-2 result object

**Files:**
- Create: `R/summary.R`
- Test: `tests/testthat/test-summary.R`

Adds 95% CIs (from the numerically-estimated Hessian of the objective) and a tidy result table
matching the workbook's parameter/comparison block (B10:P20).

- [x] **Step 1: Write the failing test**

`tests/testthat/test-summary.R`:
```r
test_that("param_ci returns finite lower/upper bounds for free parameters", {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1) + rep(c(-2, 2), length.out = 9)
  fit <- fit_model(g, "CA", "reference", "continuous",
                   start = c(max = 800, slope1 = 4, slope2 = 1.5,
                             ec501 = 0.08, ec502 = 1))
  ci <- param_ci(fit, g, "CA", "reference", "continuous")
  expect_true(all(c("parameter", "estimate", "lower", "upper") %in% names(ci)))
  expect_true(all(ci$lower <= ci$estimate & ci$estimate <= ci$upper))
})

test_that("result_table assembles one column per fitted model", {
  g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
  g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
  res <- analyse_mixture(g, "CA", "continuous")
  tab <- result_table(res)
  expect_true(all(c("reference", "SA", "DR", "DL") %in% colnames(tab)))
  expect_true("max" %in% rownames(tab))
  expect_true("objective" %in% rownames(tab))
})
```

- [x] **Step 2: Run test to verify it fails**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-summary.R")'`
Expected: FAIL — `could not find function "param_ci"`.

- [x] **Step 3: Write minimal implementation**

`R/summary.R`:
```r
#' Approximate 95% confidence intervals for fitted parameters
#'
#' Uses the numerically-estimated Hessian of the objective at the optimum.
#' For continuous data the objective is rescaled to the log-likelihood
#' (sigma^2 = SS/n) so the Hessian yields standard errors; for binary data the
#' deviance Hessian is used directly.
#' @param fit A [fit_model()] result.
#' @param df,reference,deviation,response As passed to [fit_model()].
#' @param level Confidence level.
#' @return A data frame: `parameter`, `estimate`, `lower`, `upper`.
#' @export
param_ci <- function(fit, df, reference, deviation, response, level = 0.95) {
  n_chem <- sum(c("C1", "C2", "C3") %in% names(df))
  spec <- model_spec(reference, deviation, n_chem)
  free <- names(fit$par)
  conc <- df[intersect(c("C1", "C2", "C3"), names(df))]
  names(conc) <- paste0("c", seq_along(conc))

  nll <- function(theta) {
    args <- c(as.list(conc), as.list(stats::setNames(theta, free)))
    names(args) <- c(names(conc), free)
    pred <- do.call(spec$fn, args)
    if (response == "continuous") {
      n <- nrow(df); ss <- obj_ss(df$Res, pred)
      0.5 * n * log(ss / n)              # profile Gaussian negative log-likelihood
    } else {
      0.5 * obj_deviance(df$Exposed, df$Affected,
                         pmin(pmax(pred, 1e-8), 1 - 1e-8))
    }
  }
  H <- numDeriv::hessian(nll, fit$par)
  se <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(free)))
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(parameter = free, estimate = unname(fit$par),
             lower = unname(fit$par) - z * se,
             upper = unname(fit$par) + z * se)
}

#' Assemble the Table-2 style result block from an analysis
#' @param res An [analyse_mixture()] result.
#' @return A matrix: parameters + `objective`/`df` in rows, models in columns.
#' @export
result_table <- function(res) {
  models <- names(res$fits)
  all_params <- unique(unlist(lapply(res$fits, function(f) names(f$par))))
  rows <- c(all_params, "objective", "df")
  tab <- matrix(NA_real_, nrow = length(rows), ncol = length(models),
                dimnames = list(rows, models))
  for (m in models) {
    f <- res$fits[[m]]
    tab[names(f$par), m] <- f$par
    tab["objective", m] <- f$objective
    tab["df", m] <- f$df
  }
  tab
}
```

- [x] **Step 4: Run test to verify it passes**

Run: `R -q -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-summary.R")'`
Expected: PASS (2 tests).

- [x] **Step 5: Commit**

```bash
git add R/summary.R tests/testthat/test-summary.R
git commit -m "feat: parameter CIs and Table-2 result assembly"
```

---

### Task 13: Document, generate NAMESPACE, and full check

**Files:**
- Modify: all `R/*.R` (ensure roxygen blocks present)
- Generated: `NAMESPACE`, `man/*.Rd`

- [x] **Step 1: Generate documentation and NAMESPACE**

Run: `R -q -e 'devtools::document()'`
Expected: writes `man/*.Rd` and an updated `NAMESPACE` exporting the public functions
(`ll3_predict`, `fit_single`, `obj_ss`, `binlik`, `obj_deviance`, `fit_model`, `lr_test`,
`select_parsimonious`, `analyse_mixture`, `analyse_single`, `param_ci`, `result_table`).

- [x] **Step 2: Run the full test suite**

Run: `R -q -e 'devtools::test()'`
Expected: all tests across all files PASS, 0 failures.

- [x] **Step 3: R CMD check (engine only)**

Run: `R -q -e 'devtools::check(args = "--no-manual", error_on = "warning")'`
Expected: 0 errors, 0 warnings (notes acceptable).

- [x] **Step 4: Commit**

```bash
git add NAMESPACE man R
git commit -m "docs: roxygen documentation and generated NAMESPACE for the engine"
```

---

## Engine complete

At this point `mixdra` is a working, tested R package whose engine can be driven from a plain
R script:

```r
devtools::load_all()
df  <- read.csv("my_binary_data.csv")          # columns C1, C2, Res
res <- analyse_mixture(df, reference = "CA", response = "continuous")
res$chosen
result_table(res)
```

**What was "Plan 2" has been split into focused, sequential plans** (each its own spec +
implementation plan), all built on this engine:

- **Plan 2 — plotting layer (single + binary).** The four plots (dose-response, observed-vs-
  predicted, 3-D `plotly` surface, 2-D isobole/contour) as exported package functions, all
  `plotly`, ported from Skylar's validated binary script. Design:
  `docs/design/specs/2026-05-30-mixdra-binary-plots-design.md`.
- **Plan 3 — ternary plotting.** EC50 isoplane surface + ΣTU / z-value plots (needs an
  EC50-isobole solver over the simplex + ΣTU computation; not a simple slice).
- **Plan 4 — Shiny app.** `bslib` app + stage modules + `inst/app/run_app()`, hosting the plots.
- **Plan 5 — file I/O.** Downloadable Excel/CSV templates and upload/export.
- **Plan 6 — docs.** README + methodology vignette.
