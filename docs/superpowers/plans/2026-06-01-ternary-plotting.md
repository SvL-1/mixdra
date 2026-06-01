# Ternary Plotting (T2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute the EC50 isoplane and ΣTU/z curves from a fitted ternary Advanced-S/A model and render them as interactive plotly plots, validated against the precomputed `Isoplains_Mixtox.xlsx` reference.

**Architecture:** Three layers (Approach B from the spec). A new pure, exported computation module `R/ternary-isoplane.R` (`ec50_isoplane`, `sigma_tu_curve`, `ec50_markers`) consumes an `analyse_ternary()` result and produces data frames of scientific outputs. Pure builders in `R/plot-data-ternary.R` reshape those for plotting. Thin plotly renderers in `R/plot.R` (`plot_isoplane`, `plot_sigma_tu`) consume the builders. Computation is validated numerically against the Excel fixture; renderers are smoke-tested.

**Tech Stack:** R package (`devtools`/`testthat`), `plotly` (Suggests, guarded), `readxl` (fixture extraction). Spec: `docs/superpowers/specs/2026-06-01-ternary-plotting-design.md`.

**The math (settled):** At the EC50 surface the implicit CA equation collapses to a closed form. For a TU-direction `z = (z1,z2,z3)` with `Σz = 1`:

```
F4(z)    = exp(A1*z1*z2 + A2*z1*z3 + A3*z2*z3 + A4*z1*z2*z3)
sigma_tu = F4(z)                 (ΣTU on the EC50 isoplane)
C_i      = EC50_i * z_i * F4(z)  (the isoplane point for direction z)
```

Model `"SA"` uses `A4 = 0`; `"ASA"` uses an A4. All A = 0 ⇒ `F4 ≡ 1` (pure-CA reference: ΣTU ≡ 1, `C_i = EC50_i * z_i`).

**The `analyse_ternary()` result (`res`) shape** (from `R/fit-ternary-asa.R`, already implemented):
- `res$base` — named numeric: `max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3`.
- `res$pairwise` — named numeric: `A1, A2, A3`.
- `res$A4_overall` — scalar.
- `res$individual` — data frame: `ratio, C1, C2, C3, A4, n`, where `C1/C2/C3` are **concentration proportions** (`C_i / ΣC`, from `ternary_ratio_key`), NOT TU fractions.

## Running R and tests here (read first)

- Required packages live in the **user** library; non-interactive R does NOT auto-load them. Prefix every R invocation with `.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))`.
- Run a single test file headlessly (writes a marker to a file, then Read it — the terminal mangles streamed stdout):
  ```r
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  devtools::load_all(".")
  testthat::test_file("tests/testthat/test-ternary-isoplane.R")
  ```
  Or in RStudio: open `mixdra.Rproj`, Ctrl+Shift+L (load_all), then run the file.
- `plotly` and `readxl` are installed in the user library.
- `R CMD check` / `rcmdcheck` CANNOT run here (no Rtools). Do not gate tasks on it; it is a release-time step.
- Commit after every task. Do NOT add a `Co-Authored-By: Claude` trailer (project rule).

## File structure

```
R/ternary-isoplane.R                               NEW  exported: ec50_isoplane, sigma_tu_curve, ec50_markers; internal: simplex_grid
R/plot-data-ternary.R                              NEW  pure builders: isoplane_plot_data, sigma_tu_plot_data
R/plot.R                                           EDIT add plot_isoplane, plot_sigma_tu (reuses require_plotly)
tests/testthat/fixtures/extract_isoplane_fixture.R NEW  tracked extraction tool
tests/testthat/fixtures/isoplane_points.csv        NEW  (generated)
tests/testthat/fixtures/isoplane_zvalues.csv       NEW  (generated)
tests/testthat/fixtures/isoplane_ec50.csv          NEW  (generated)
tests/testthat/test-ternary-isoplane.R             NEW  unit tests for the API
tests/testthat/test-plot-data-ternary.R            NEW  unit tests for builders
tests/testthat/test-plot-ternary.R                 NEW  renderer smoke tests
tests/testthat/test-validation-isoplane.R          NEW  validation vs Excel fixture
NAMESPACE, man/*.Rd                                 EDIT via roxygen2 document()
```

---

> **AS-BUILT amendment (2026-06-01).** Tasks 1 and 8 changed during execution.
> The Excel reference `Isoplains_Mixtox.xlsx` was found to compute its isoplane by
> a different, interior-approximate method (~9% interior divergence), while our
> `ec50_isoplane` is provably exact for the fitted model. Per the project owner's
> decision, **the Excel validation was dropped**: Task 1's extracted Excel
> fixtures (`extract_isoplane_fixture.R`, `isoplane_*.csv`) were removed, and
> Task 8 now validates by **self-consistency** — every isoplane / marker / curve
> point, fed back through the `ca_asa_tri` predictor, must return `max/2`
> (tolerance `1e-4`). The Task 1 and Task 8 steps below are kept as the original
> script for the historical record; see the spec §4 as-built amendment for the
> final shape. The z-path convention work in Task 1 (still valid) is now done
> directly from the verified `TU-zValues` layout rather than at runtime.

## Task 1: Extract isoplane fixtures from `Isoplains_Mixtox.xlsx`

**Files:**
- Create: `tests/testthat/fixtures/extract_isoplane_fixture.R`
- Create (generated): `tests/testthat/fixtures/isoplane_points.csv`, `isoplane_zvalues.csv`, `isoplane_ec50.csv`

Mirrors `tests/testthat/fixtures/extract_ternary_fixture.R`. The reference workbook is `recieved/R scripts & excel/R scripts & excel/Isoplains_Mixtox.xlsx`. Verified sheet layouts (inspected 2026-06-01):
- `FBSA_CPF_IMI`: cols `CPF, FBSA, IMI, Mod`; `Mod ∈ {"CA_SA","CA_SA+SA","Exp"}`; 760 rows.
- `TU-zValues`: cols `z, TU, Chem, Mod, Exp`; FBSA/CPF/IMI mixture is `Exp == "Mixtox1"`; `Mod ∈ {"CA_SA","CA_SA+SA"}`; `Chem ∈ {"CPF","FBSA","IMI"}` for this mixture.
- `TER_EC50`: cols `z, TU, Chem, Exp, Ratio`; FBSA/CPF/IMI is `Exp == "Mixtox1"`.

- [ ] **Step 1: Write the extraction script**

Create `tests/testthat/fixtures/extract_isoplane_fixture.R`:

```r
# One-off: extract the precomputed FBSA/CPF/IMI isoplane reference data from
# Skylar's Isoplains_Mixtox.xlsx into CSV fixtures, for validating the
# model-driven ec50_isoplane / sigma_tu_curve / ec50_markers functions.
# Run from the repo root:
#   R -q -e 'source("tests/testthat/fixtures/extract_isoplane_fixture.R")'
# Sheet layouts verified 2026-06-01 (see plan Task 1).
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required to regenerate this fixture.")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
wb <- "recieved/R scripts & excel/R scripts & excel/Isoplains_Mixtox.xlsx"

# Isoplane point cloud: one row per (point, model). Concentrations map
# CPF->C1, FBSA->C2, IMI->C3 (matches the C1/C2/C3 convention of the fixture
# used by analyse_ternary, where CPF=C1, FBSA=C2, IMI=C3).
pts <- as.data.frame(readxl::read_excel(wb, sheet = "FBSA_CPF_IMI"))
pts <- pts[, c("CPF", "FBSA", "IMI", "Mod")]
names(pts) <- c("C1", "C2", "C3", "model")
pts <- pts[stats::complete.cases(pts), , drop = FALSE]
write.csv(pts, "tests/testthat/fixtures/isoplane_points.csv", row.names = FALSE)
cat("wrote", nrow(pts), "isoplane points (",
    paste(unique(pts$model), collapse = ", "), ")\n")

# z-value curves for this mixture.
zv <- as.data.frame(readxl::read_excel(wb, sheet = "TU-zValues"))
zv <- zv[zv$Exp == "Mixtox1", c("z", "TU", "Chem", "Mod"), drop = FALSE]
zv <- zv[stats::complete.cases(zv), , drop = FALSE]
write.csv(zv, "tests/testthat/fixtures/isoplane_zvalues.csv", row.names = FALSE)
cat("wrote", nrow(zv), "z-value rows\n")

# EC50 marker points for this mixture.
ec <- as.data.frame(readxl::read_excel(wb, sheet = "TER_EC50"))
ec <- ec[ec$Exp == "Mixtox1", c("z", "TU", "Chem", "Ratio"), drop = FALSE]
ec <- ec[stats::complete.cases(ec), , drop = FALSE]
write.csv(ec, "tests/testthat/fixtures/isoplane_ec50.csv", row.names = FALSE)
cat("wrote", nrow(ec), "EC50 marker rows\n")
```

- [ ] **Step 2: Run the extraction script**

Run from the repo root:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); source("tests/testthat/fixtures/extract_isoplane_fixture.R")'
```
Expected: three "wrote N ..." lines; `isoplane_points.csv` reports models `CA_SA, CA_SA+SA, Exp`; the three CSVs exist under `tests/testthat/fixtures/`.

- [ ] **Step 3: Sanity-check the fixtures**

Run:
```bash
Rscript -e 'p<-read.csv("tests/testthat/fixtures/isoplane_points.csv"); print(table(p$model)); print(head(p)); z<-read.csv("tests/testthat/fixtures/isoplane_zvalues.csv"); print(table(z$Chem,z$Mod))'
```
Expected: `model` has all three of `CA_SA`/`CA_SA+SA`/`Exp`; `isoplane_zvalues.csv` has CPF/FBSA/IMI × CA_SA/CA_SA+SA rows. Confirm `C1` (CPF) includes 0 values (binary edges of the simplex are present).

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/fixtures/extract_isoplane_fixture.R tests/testthat/fixtures/isoplane_points.csv tests/testthat/fixtures/isoplane_zvalues.csv tests/testthat/fixtures/isoplane_ec50.csv
git commit -m "test: extract FBSA/CPF/IMI isoplane reference fixtures from Isoplains_Mixtox.xlsx"
```

---

## Task 2: `simplex_grid` + `ec50_isoplane`

**Files:**
- Create: `R/ternary-isoplane.R`
- Test: `tests/testthat/test-ternary-isoplane.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-ternary-isoplane.R`:

```r
# A minimal analyse_ternary-shaped result for unit testing the pure isoplane API
# (no fitting needed — the functions only read these fields).
mock_res <- function(A1 = 0, A2 = 0, A3 = 0, A4 = 0,
                     ec50 = c(2, 5, 0.5), max = 100) {
  list(
    reference = "CA", response = "continuous",
    base = c(max = max, slope1 = 3, slope2 = 3, slope3 = 3,
             ec50_1 = ec50[1], ec50_2 = ec50[2], ec50_3 = ec50[3]),
    pairwise = c(A1 = A1, A2 = A2, A3 = A3),
    A4_overall = A4,
    individual = data.frame(
      ratio = c("a", "b"), C1 = c(1/3, 0.6), C2 = c(1/3, 0.2),
      C3 = c(1/3, 0.2), A4 = c(0.5, -0.3), n = c(5L, 5L),
      stringsAsFactors = FALSE))
}

test_that("simplex_grid has (n+1)(n+2)/2 rows summing to 1", {
  g <- simplex_grid(4)
  expect_equal(nrow(g), (4 + 1) * (4 + 2) / 2)   # 15
  expect_equal(rowSums(g[c("z1", "z2", "z3")]), rep(1, nrow(g)), tolerance = 1e-12)
  expect_true(all(g >= 0))
})

test_that("ec50_isoplane with no interaction is the pure-CA reference", {
  res <- mock_res()                       # all A = 0
  d <- ec50_isoplane(res, "SA", n = 10)
  expect_equal(d$sigma_tu, rep(1, nrow(d)), tolerance = 1e-12)
  # C_i = EC50_i * z_i when F4 = 1
  expect_equal(d$C1, res$base[["ec50_1"]] * d$z1, tolerance = 1e-12)
  expect_equal(d$C2, res$base[["ec50_2"]] * d$z2, tolerance = 1e-12)
  expect_equal(d$C3, res$base[["ec50_3"]] * d$z3, tolerance = 1e-12)
})

test_that("ec50_isoplane vertices give single-chemical EC50 points", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 2)
  d <- ec50_isoplane(res, "ASA", n = 6)
  v1 <- d[d$z1 == 1, ]                     # vertex z = (1,0,0)
  expect_equal(nrow(v1), 1)
  expect_equal(v1$sigma_tu, 1, tolerance = 1e-12)        # F4 = exp(0) = 1
  expect_equal(v1$C1, res$base[["ec50_1"]], tolerance = 1e-12)
  expect_equal(c(v1$C2, v1$C3), c(0, 0), tolerance = 1e-12)
})

test_that("ec50_isoplane SA ignores A4, ASA uses A4_overall", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 5)
  sa  <- ec50_isoplane(res, "SA", n = 8)
  asa <- ec50_isoplane(res, "ASA", n = 8)
  # interior point with all z > 0 differs between SA and ASA; edges (a z == 0) match
  interior <- sa$z1 > 0 & sa$z2 > 0 & sa$z3 > 0
  expect_false(isTRUE(all.equal(sa$sigma_tu[interior], asa$sigma_tu[interior])))
  expect_equal(sa$sigma_tu[!interior], asa$sigma_tu[!interior], tolerance = 1e-12)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-ternary-isoplane.R")'
```
Expected: FAIL — `could not find function "simplex_grid"` / `"ec50_isoplane"`.

- [ ] **Step 3: Write the implementation**

Create `R/ternary-isoplane.R`:

```r
# Model-driven EC50 isoplane and Sigma-TU / z computation for the ternary
# Advanced-S/A fit. Pure (no plotly); these are scientific outputs in their own
# right, consumed by the plot builders in R/plot-data-ternary.R and reusable by
# the app. At the EC50 surface the implicit CA equation collapses to a closed
# form (see plan/spec): on the isoplane, sigma_tu = F4(z) and C_i = EC50_i*z_i*F4.

#' Triangular grid over the 2-simplex
#'
#' All `(i, j, k)` with `i + j + k = n` and `i, j, k >= 0`, scaled to
#' proportions summing to 1. Vertices (e.g. `(1,0,0)`) and edges (one coordinate
#' 0) are included.
#' @param n Grid resolution (points per simplex edge).
#' @return A data frame with columns `z1`, `z2`, `z3`.
#' @keywords internal
simplex_grid <- function(n) {
  rows <- list()
  for (i in 0:n) for (j in 0:(n - i)) rows[[length(rows) + 1L]] <- c(i, j, n - i - j)
  m <- do.call(rbind, rows) / n
  data.frame(z1 = m[, 1], z2 = m[, 2], z3 = m[, 3])
}

#' F4 deviation factor for a ternary Advanced-S/A model
#' @keywords internal
.f4 <- function(z1, z2, z3, A1, A2, A3, A4) {
  exp(A1 * z1 * z2 + A2 * z1 * z3 + A3 * z2 * z3 + A4 * z1 * z2 * z3)
}

#' EC50 isoplane points for a ternary Advanced-S/A fit
#'
#' Returns the EC50 surface in concentration space over a triangular grid of
#' TU-directions `z` on the 2-simplex. On the EC50 surface `sigma_tu = F4(z)` and
#' `C_i = EC50_i * z_i * F4(z)`.
#' @param res An [analyse_ternary()] result.
#' @param model `"SA"` (pairwise only, A4 = 0) or `"ASA"` (adds the overall
#'   three-way `res$A4_overall`).
#' @param n Grid resolution per simplex edge (default 30).
#' @return A data frame: `z1, z2, z3, C1, C2, C3, sigma_tu, model`.
#' @export
ec50_isoplane <- function(res, model = c("SA", "ASA"), n = 30) {
  model <- match.arg(model)
  b <- res$base
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]
  A4 <- if (model == "ASA") res$A4_overall else 0
  g <- simplex_grid(n)
  F4 <- .f4(g$z1, g$z2, g$z3, A1, A2, A3, A4)
  data.frame(z1 = g$z1, z2 = g$z2, z3 = g$z3,
             C1 = b[["ec50_1"]] * g$z1 * F4,
             C2 = b[["ec50_2"]] * g$z2 * F4,
             C3 = b[["ec50_3"]] * g$z3 * F4,
             sigma_tu = F4, model = model, stringsAsFactors = FALSE)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-ternary-isoplane.R")'
```
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add R/ternary-isoplane.R tests/testthat/test-ternary-isoplane.R
git commit -m "feat: ec50_isoplane + simplex_grid for ternary Advanced-S/A"
```

---

## Task 3: `sigma_tu_curve`

**Files:**
- Modify: `R/ternary-isoplane.R`
- Test: `tests/testthat/test-ternary-isoplane.R`

z-path convention (verified against `TU-zValues` at `z_x = 0`, 2026-06-01): for chemical *x*, vary `z_x` from 0→1 while the **other two are held at equal z** (`(1 - z_x)/2` each). At `z_x = 1` it is the pure-x vertex (`F4 = 1`); at `z_x = 0` it is the equal-split binary of the other two.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ternary-isoplane.R`:

```r
test_that("sigma_tu_curve: no interaction gives sigma_tu == 1 everywhere", {
  res <- mock_res()                       # all A = 0
  d <- sigma_tu_curve(res, "SA", n = 11)
  expect_setequal(unique(d$chem), c("C1", "C2", "C3"))
  expect_equal(d$sigma_tu, rep(1, nrow(d)), tolerance = 1e-12)
})

test_that("sigma_tu_curve: each chemical's vertex (z = 1) has sigma_tu == 1", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 2)
  d <- sigma_tu_curve(res, "ASA", n = 11)
  tip <- d[d$z == 1, ]
  expect_equal(nrow(tip), 3)              # one per chemical
  expect_equal(tip$sigma_tu, rep(1, 3), tolerance = 1e-12)
})

test_that("sigma_tu_curve: C1 at z = 0 is the equal-split FBSA/IMI binary point", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 2)
  d <- sigma_tu_curve(res, "SA", n = 11)
  c1_0 <- d[d$chem == "C1" & d$z == 0, "sigma_tu"]
  # other two equal at z = 0.5 each => F4 = exp(A3 * 0.5 * 0.5)
  expect_equal(c1_0, exp(res$pairwise[["A3"]] * 0.25), tolerance = 1e-12)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-ternary-isoplane.R")'
```
Expected: FAIL — `could not find function "sigma_tu_curve"`.

- [ ] **Step 3: Write the implementation**

Append to `R/ternary-isoplane.R`:

```r
#' Sigma-TU vs z curves for a ternary Advanced-S/A fit
#'
#' For each chemical *x*, varies its TU-fraction `z_x` from 0 to 1 while holding
#' the other two chemicals at equal z (`(1 - z_x)/2` each), and reports `ΣTU =
#' F4(z)` along that path. This is the curve underlying the z-plot: deviation of
#' `sigma_tu` from 1 is the interaction (`< 1` synergy, `> 1` antagonism).
#' @param res An [analyse_ternary()] result.
#' @param model `"SA"` (pairwise only) or `"ASA"` (adds `res$A4_overall`).
#' @param n Number of z points per chemical (default 21, i.e. steps of 0.05).
#' @return A data frame: `chem` ("C1"/"C2"/"C3"), `z`, `sigma_tu`, `model`.
#' @export
sigma_tu_curve <- function(res, model = c("SA", "ASA"), n = 21) {
  model <- match.arg(model)
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]
  A4 <- if (model == "ASA") res$A4_overall else 0
  zx <- seq(0, 1, length.out = n)
  other <- (1 - zx) / 2
  chems <- c("C1", "C2", "C3")
  do.call(rbind, lapply(seq_along(chems), function(ix) {
    z <- matrix(other, nrow = n, ncol = 3)
    z[, ix] <- zx
    data.frame(chem = chems[ix], z = zx,
               sigma_tu = .f4(z[, 1], z[, 2], z[, 3], A1, A2, A3, A4),
               model = model, stringsAsFactors = FALSE)
  }))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-ternary-isoplane.R")'
```
Expected: PASS (7 tests).

- [ ] **Step 5: Commit**

```bash
git add R/ternary-isoplane.R tests/testthat/test-ternary-isoplane.R
git commit -m "feat: sigma_tu_curve (Sigma-TU vs z) for ternary Advanced-S/A"
```

---

## Task 4: `ec50_markers`

**Files:**
- Modify: `R/ternary-isoplane.R`
- Test: `tests/testthat/test-ternary-isoplane.R`

`res$individual$C1/C2/C3` are concentration proportions `p_i`, so convert to TU-fractions `z_i = (p_i/EC50_i) / Σ(p_j/EC50_j)` before applying the isoplane formula with that ratio's individual `A4`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ternary-isoplane.R`:

```r
test_that("ec50_markers returns one row per individual ratio", {
  res <- mock_res(A1 = 0.7, A2 = -0.3, A3 = -1.8)
  m <- ec50_markers(res, df = NULL)
  expect_equal(nrow(m), nrow(res$individual))
  expect_equal(m$ratio, res$individual$ratio)
  expect_true(all(is.finite(m$sigma_tu)))
  expect_true(all(c("z1", "z2", "z3", "C1", "C2", "C3", "sigma_tu") %in% names(m)))
})

test_that("ec50_markers TU-fraction conversion: equal proportions are NOT equal z", {
  # ratio "a" has equal concentration proportions (1/3 each) but unequal EC50s,
  # so its z (TU fractions) must be unequal and weighted by 1/EC50_i.
  res <- mock_res(A1 = 0, A2 = 0, A3 = 0, ec50 = c(2, 5, 0.5))
  m <- ec50_markers(res, df = NULL)
  ra <- m[m$ratio == "a", ]
  ec <- c(2, 5, 0.5); z_expected <- (1 / ec) / sum(1 / ec)
  expect_equal(c(ra$z1, ra$z2, ra$z3), z_expected, tolerance = 1e-12)
  expect_equal(ra$sigma_tu, 1, tolerance = 1e-12)   # all A = 0 here
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-ternary-isoplane.R")'
```
Expected: FAIL — `could not find function "ec50_markers"`.

- [ ] **Step 3: Write the implementation**

Append to `R/ternary-isoplane.R`:

```r
#' EC50 marker points for each tested ternary ratio
#'
#' For every ratio in `res$individual`, converts its concentration proportions to
#' TU-fractions `z` and evaluates the EC50 isoplane point using that ratio's own
#' individual `A4`. Plotting these against the overall isoplane/curve is the
#' visual form of the per-ratio-vs-overall A4 "averaging-out" comparison.
#' @param res An [analyse_ternary()] result.
#' @param df Unused (kept for signature symmetry with the renderers); the marker
#'   geometry comes entirely from `res$individual` and `res$base`.
#' @return A data frame: `ratio, z1, z2, z3, C1, C2, C3, sigma_tu`.
#' @export
ec50_markers <- function(res, df = NULL) {
  b <- res$base
  ec <- c(b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]])
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]
  ind <- res$individual
  if (nrow(ind) == 0)
    return(data.frame(ratio = character(0), z1 = numeric(0), z2 = numeric(0),
                      z3 = numeric(0), C1 = numeric(0), C2 = numeric(0),
                      C3 = numeric(0), sigma_tu = numeric(0),
                      stringsAsFactors = FALSE))
  do.call(rbind, lapply(seq_len(nrow(ind)), function(i) {
    p <- c(ind$C1[i], ind$C2[i], ind$C3[i])    # concentration proportions
    tu <- p / ec
    z <- tu / sum(tu)
    F4 <- .f4(z[1], z[2], z[3], A1, A2, A3, ind$A4[i])
    data.frame(ratio = ind$ratio[i], z1 = z[1], z2 = z[2], z3 = z[3],
               C1 = ec[1] * z[1] * F4, C2 = ec[2] * z[2] * F4,
               C3 = ec[3] * z[3] * F4, sigma_tu = F4, stringsAsFactors = FALSE)
  }))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-ternary-isoplane.R")'
```
Expected: PASS (9 tests).

- [ ] **Step 5: Commit**

```bash
git add R/ternary-isoplane.R tests/testthat/test-ternary-isoplane.R
git commit -m "feat: ec50_markers (per-ratio EC50 points) for ternary Advanced-S/A"
```

---

## Task 5: Roxygen documentation + NAMESPACE exports

**Files:**
- Modify: `NAMESPACE` (generated), `man/*.Rd` (generated)

- [ ] **Step 1: Regenerate docs and NAMESPACE**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::document(".")'
```
Expected: writes `man/ec50_isoplane.Rd`, `man/sigma_tu_curve.Rd`, `man/ec50_markers.Rd`; adds `export(ec50_isoplane)`, `export(sigma_tu_curve)`, `export(ec50_markers)` to `NAMESPACE`. `simplex_grid` and `.f4` are NOT exported (`@keywords internal` / no `@export`).

- [ ] **Step 2: Verify exports**

Run:
```bash
Rscript -e 'cat(grep("isoplane|sigma_tu|ec50_markers", readLines("NAMESPACE"), value = TRUE), sep = "\n")'
```
Expected: three `export(...)` lines for the new functions; no `export(simplex_grid)`.

- [ ] **Step 3: Commit**

```bash
git add NAMESPACE man/
git commit -m "docs: man pages + NAMESPACE for ternary isoplane API"
```

---

## Task 6: Plot builders (`R/plot-data-ternary.R`)

**Files:**
- Create: `R/plot-data-ternary.R`
- Test: `tests/testthat/test-plot-data-ternary.R`

Builders assemble the API output into a single long-format frame per plot, so the renderers stay trivial. `isoplane_plot_data` stacks the SA cloud, ASA cloud, and EC50 markers with a `series` label. `sigma_tu_plot_data` stacks the per-chemical SA and ASA curves with a `series` label (the ΣTU = 1 reference line is added by the renderer).

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-plot-data-ternary.R`:

```r
mock_res <- function(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 1,
                     ec50 = c(2, 5, 0.5), max = 100) {
  list(reference = "CA", response = "continuous",
       base = c(max = max, slope1 = 3, slope2 = 3, slope3 = 3,
                ec50_1 = ec50[1], ec50_2 = ec50[2], ec50_3 = ec50[3]),
       pairwise = c(A1 = A1, A2 = A2, A3 = A3), A4_overall = A4,
       individual = data.frame(ratio = c("a", "b"), C1 = c(1/3, 0.6),
                               C2 = c(1/3, 0.2), C3 = c(1/3, 0.2),
                               A4 = c(0.5, -0.3), n = c(5L, 5L),
                               stringsAsFactors = FALSE))
}

test_that("isoplane_plot_data stacks SA, ASA and markers with a series label", {
  res <- mock_res()
  d <- isoplane_plot_data(res, df = NULL, n = 6)
  expect_true(all(c("C1", "C2", "C3", "series") %in% names(d)))
  expect_setequal(unique(d$series), c("CA+S/A", "CA+S/A+S/A", "EC50"))
  expect_equal(sum(d$series == "EC50"), nrow(res$individual))
})

test_that("sigma_tu_plot_data stacks SA and ASA per-chemical curves", {
  res <- mock_res()
  d <- sigma_tu_plot_data(res, n = 11)
  expect_true(all(c("chem", "z", "sigma_tu", "series") %in% names(d)))
  expect_setequal(unique(d$series), c("CA+S/A", "CA+S/A+S/A"))
  expect_setequal(unique(d$chem), c("C1", "C2", "C3"))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-plot-data-ternary.R")'
```
Expected: FAIL — `could not find function "isoplane_plot_data"`.

- [ ] **Step 3: Write the implementation**

Create `R/plot-data-ternary.R`:

```r
# Pure builders for the ternary plotting layer. No plotly here; they reshape the
# R/ternary-isoplane.R outputs into single long-format frames with a `series`
# label so the renderers in R/plot.R stay thin.

#' Combine isoplane clouds + EC50 markers into one plotting frame
#'
#' Stacks the CA+S/A isoplane, the CA+S/A+S/A isoplane, and the per-ratio EC50
#' markers, tagging each with a `series` label.
#' @param res An [analyse_ternary()] result.
#' @param df Passed through to [ec50_markers()] (unused there).
#' @param n Isoplane grid resolution.
#' @return A data frame: `C1, C2, C3, series`.
#' @keywords internal
isoplane_plot_data <- function(res, df = NULL, n = 30) {
  sa  <- ec50_isoplane(res, "SA", n)
  asa <- ec50_isoplane(res, "ASA", n)
  mk  <- ec50_markers(res, df)
  rbind(
    data.frame(C1 = sa$C1,  C2 = sa$C2,  C3 = sa$C3,  series = "CA+S/A",
               stringsAsFactors = FALSE),
    data.frame(C1 = asa$C1, C2 = asa$C2, C3 = asa$C3, series = "CA+S/A+S/A",
               stringsAsFactors = FALSE),
    data.frame(C1 = mk$C1,  C2 = mk$C2,  C3 = mk$C3,  series = "EC50",
               stringsAsFactors = FALSE))
}

#' Combine SA and ASA Sigma-TU curves into one plotting frame
#' @param res An [analyse_ternary()] result.
#' @param n Number of z points per chemical.
#' @return A data frame: `chem, z, sigma_tu, series`.
#' @keywords internal
sigma_tu_plot_data <- function(res, n = 21) {
  sa  <- sigma_tu_curve(res, "SA", n)
  asa <- sigma_tu_curve(res, "ASA", n)
  rbind(
    data.frame(chem = sa$chem,  z = sa$z,  sigma_tu = sa$sigma_tu,
               series = "CA+S/A", stringsAsFactors = FALSE),
    data.frame(chem = asa$chem, z = asa$z, sigma_tu = asa$sigma_tu,
               series = "CA+S/A+S/A", stringsAsFactors = FALSE))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-plot-data-ternary.R")'
```
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add R/plot-data-ternary.R tests/testthat/test-plot-data-ternary.R
git commit -m "feat: ternary plot-data builders (isoplane + sigma-TU)"
```

---

## Task 7: Renderers `plot_isoplane` + `plot_sigma_tu`

**Files:**
- Modify: `R/plot.R`
- Test: `tests/testthat/test-plot-ternary.R`

Reuses the existing `require_plotly()` helper at the top of `R/plot.R`. Colors mirror Skylar's Isoplane figures (CA+S/A orange, CA+S/A+S/A blue, EC50 markers red).

- [ ] **Step 1: Write the failing smoke tests**

Create `tests/testthat/test-plot-ternary.R`:

```r
mock_res <- function(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 1,
                     ec50 = c(2, 5, 0.5), max = 100) {
  list(reference = "CA", response = "continuous",
       base = c(max = max, slope1 = 3, slope2 = 3, slope3 = 3,
                ec50_1 = ec50[1], ec50_2 = ec50[2], ec50_3 = ec50[3]),
       pairwise = c(A1 = A1, A2 = A2, A3 = A3), A4_overall = A4,
       individual = data.frame(ratio = c("a", "b"), C1 = c(1/3, 0.6),
                               C2 = c(1/3, 0.2), C3 = c(1/3, 0.2),
                               A4 = c(0.5, -0.3), n = c(5L, 5L),
                               stringsAsFactors = FALSE))
}

test_that("plot_isoplane returns a plotly object with three traces", {
  skip_if_not_installed("plotly")
  p <- plot_isoplane(mock_res(), df = NULL, n = 8)
  expect_s3_class(p, "plotly")
  expect_equal(length(plotly::plotly_build(p)$x$data), 3)  # SA, ASA, markers
})

test_that("plot_sigma_tu returns a plotly object with the reference line", {
  skip_if_not_installed("plotly")
  p <- plot_sigma_tu(mock_res(), n = 11)
  expect_s3_class(p, "plotly")
  # 3 chems x 2 models = 6 line traces + 1 reference line = 7
  expect_equal(length(plotly::plotly_build(p)$x$data), 7)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-plot-ternary.R")'
```
Expected: FAIL — `could not find function "plot_isoplane"`.

- [ ] **Step 3: Write the implementation**

Append to `R/plot.R`:

```r
#' Plot the EC50 isoplane of a ternary mixture
#'
#' 3-D scatter of the CA+S/A and CA+S/A+S/A EC50 isoplane point-clouds with the
#' per-ratio EC50 markers overlaid, as an interactive plotly object. Mirrors the
#' MixTox isoplane figures, sourced from the fitted model.
#' @param res An [analyse_ternary()] result.
#' @param df The data frame the fit was built from (passed to [ec50_markers()]).
#' @param n Isoplane grid resolution per simplex edge (default 30).
#' @return A plotly object.
#' @export
plot_isoplane <- function(res, df = NULL, n = 30) {
  require_plotly()
  d <- isoplane_plot_data(res, df, n = n)
  cols <- c("CA+S/A" = "orange", "CA+S/A+S/A" = "blue", "EC50" = "red")
  sizes <- c("CA+S/A" = 3, "CA+S/A+S/A" = 3, "EC50" = 7)
  p <- plotly::plot_ly()
  for (s in c("CA+S/A", "CA+S/A+S/A", "EC50")) {
    seg <- d[d$series == s, ]
    p <- plotly::add_trace(p, x = seg$C1, y = seg$C2, z = seg$C3,
                           type = "scatter3d", mode = "markers", name = s,
                           marker = list(size = sizes[[s]], color = cols[[s]]))
  }
  plotly::layout(p, scene = list(xaxis = list(title = "C1"),
                                 yaxis = list(title = "C2"),
                                 zaxis = list(title = "C3")))
}

#' Plot Sigma-TU vs z for a ternary mixture
#'
#' Per-chemical ΣTU-vs-z curves under CA+S/A (solid) and CA+S/A+S/A (dashed),
#' with the additivity reference line at ΣTU = 1, as an interactive plotly
#' object. Deviation from 1 is the interaction (< 1 synergy, > 1 antagonism).
#' @param res An [analyse_ternary()] result.
#' @param n Number of z points per chemical (default 21).
#' @return A plotly object.
#' @export
plot_sigma_tu <- function(res, n = 21) {
  require_plotly()
  d <- sigma_tu_plot_data(res, n = n)
  p <- plotly::plot_ly()
  for (ch in c("C1", "C2", "C3")) {
    for (s in c("CA+S/A", "CA+S/A+S/A")) {
      seg <- d[d$chem == ch & d$series == s, ]
      p <- plotly::add_lines(p, x = seg$z, y = seg$sigma_tu,
                             name = paste(ch, s),
                             line = list(dash = if (s == "CA+S/A") "solid" else "dash"))
    }
  }
  p <- plotly::add_lines(p, x = c(0, 1), y = c(1, 1), name = "ΣTU = 1",
                         line = list(color = "black", width = 1))
  plotly::layout(p, xaxis = list(title = "z (chemical TU fraction)"),
                 yaxis = list(title = "ΣTU"))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); testthat::test_file("tests/testthat/test-plot-ternary.R")'
```
Expected: PASS (2 tests).

- [ ] **Step 5: Regenerate docs and commit**

Run:
```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::document(".")'
git add R/plot.R man/ NAMESPACE tests/testthat/test-plot-ternary.R
git commit -m "feat: plot_isoplane + plot_sigma_tu plotly renderers"
```

---

## Task 8: Validation against the Excel fixture

**Files:**
- Test: `tests/testthat/test-validation-isoplane.R`

Validates the computation numerically against the precomputed Excel data (Task 1 fixtures). Robust strategy: rather than matching grid parameterizations, verify the **identity** holds at the Excel's own isoplane points — for each `CA_SA` reference point `(C1,C2,C3)`, recompute `z` and `ΣTU` from our fitted `EC50_i` and assert `ΣTU ≈ F4_SA(z)` with our fitted `A1/A2/A3`. The fit itself was already validated in T1, so a pass here confirms the isoplane math + conventions. Also asserts the z-curve reproduces `TU-zValues` and markers land near `TER_EC50`.

SLOW: calls `analyse_ternary` on the 419-row fixture (~2 min). Run in the background and Read the marker, per the project's R-here gotchas.

- [ ] **Step 1: Write the validation test**

Create `tests/testthat/test-validation-isoplane.R`:

```r
# Validation gate: does the model-driven isoplane/Sigma-TU computation reproduce
# Skylar's precomputed Isoplains_Mixtox.xlsx reference for the FBSA/CPF/IMI
# mixture? The ternary fit is already validated in test-validation-ternary.R, so
# this isolates the isoplane math + z-path conventions.
#
# CPF=C1, FBSA=C2, IMI=C3 (matches extract_isoplane_fixture.R and the ternary
# fixture). SLOW: analyse_ternary on 419 rows ~ 2 min.

test_that("ec50_isoplane / sigma_tu_curve reproduce Isoplains_Mixtox.xlsx", {
  fx <- testthat::test_path("fixtures", "ternary_fbsa_cpf_imi_continuous.csv")
  pts_f <- testthat::test_path("fixtures", "isoplane_points.csv")
  zv_f  <- testthat::test_path("fixtures", "isoplane_zvalues.csv")
  skip_if_not(file.exists(fx) && file.exists(pts_f) && file.exists(zv_f),
              "isoplane fixtures missing")

  df <- read.csv(fx)
  set.seed(42)
  res <- analyse_ternary(df, "CA", "continuous", n_starts = 10, time_limit = 120,
                         lower = c(ec50_2 = 5.579), upper = c(ec50_2 = 5.581))
  ec <- c(res$base[["ec50_1"]], res$base[["ec50_2"]], res$base[["ec50_3"]])
  A1 <- res$pairwise[["A1"]]; A2 <- res$pairwise[["A2"]]; A3 <- res$pairwise[["A3"]]

  ## --- Isoplane identity at the Excel's CA_SA points -------------------------
  pts <- read.csv(pts_f)
  sa <- pts[pts$model == "CA_SA", ]
  # recompute z and sigma_tu from our EC50s at each reference point
  tu <- cbind(sa$C1 / ec[1], sa$C2 / ec[2], sa$C3 / ec[3])
  s  <- rowSums(tu)
  z  <- tu / s
  f4 <- exp(A1 * z[, 1] * z[, 2] + A2 * z[, 1] * z[, 3] + A3 * z[, 2] * z[, 3])
  # On the EC50 isoplane sigma_tu (= s) must equal F4(z). Median relative error
  # tolerance absorbs the workbook's rounding and small A-parameter differences.
  rel <- abs(s - f4) / pmax(f4, 1e-6)
  expect_lt(stats::median(rel), 0.02)

  ## --- Sigma-TU curve vs TU-zValues (CA_SA) ----------------------------------
  zv <- read.csv(zv_f)
  zv_sa <- zv[zv$Mod == "CA_SA", ]
  chem_map <- c(CPF = "C1", FBSA = "C2", IMI = "C3")
  cur <- sigma_tu_curve(res, "SA", n = 21)
  errs <- vapply(seq_len(nrow(zv_sa)), function(i) {
    ch <- chem_map[[ zv_sa$Chem[i] ]]
    j  <- which.min(abs(cur$z[cur$chem == ch] - zv_sa$z[i]))
    ours <- cur$sigma_tu[cur$chem == ch][j]
    abs(ours - zv_sa$TU[i]) / max(zv_sa$TU[i], 1e-6)
  }, numeric(1))
  expect_lt(stats::median(errs), 0.05)
})
```

- [ ] **Step 2: Run the validation test in the background**

The run is long. Use a driver that writes a marker, run detached, then Read the marker file:

```r
# diag_run_iso_val.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))
out <- "diag_iso_val_out.txt"
writeLines("VAL: starting", out)
r <- as.data.frame(testthat::test_file(
  "tests/testthat/test-validation-isoplane.R", reporter = "summary"))
cat(sprintf("VERDICT: fail=%d warn=%d\n",
            sum(r$failed, na.rm = TRUE), sum(r$warning, na.rm = TRUE)),
    file = out, append = TRUE)
```
Run:
```bash
Rscript diag_run_iso_val.R > diag_iso_val_console.txt 2>&1
```
Expected (in `diag_iso_val_out.txt`): `VERDICT: fail=0 warn=0`.

If the median-relative-error assertions fail, the likely cause is the z-path convention or a CPF/FBSA/IMI ↔ C1/C2/C3 mismatch — inspect a few rows of `isoplane_points.csv` against `ec50_isoplane(res, "SA")` before adjusting tolerances. Do NOT loosen tolerances to mask a real convention bug.

- [ ] **Step 3: Clean up the driver scratch**

```bash
rm -f diag_run_iso_val.R diag_iso_val_out.txt diag_iso_val_console.txt
```
(These match the `diag_*` gitignore pattern, so they never enter `git status`.)

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-validation-isoplane.R
git commit -m "test: validate ternary isoplane/sigma-TU against Isoplains_Mixtox.xlsx"
```

---

## Task 9: Full-suite check + sync

**Files:** none (verification + git).

- [ ] **Step 1: Run the full suite in the background**

```r
# diag_run_suite.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(devtools))
out <- "diag_suite_out.txt"; writeLines("VAL: starting", out)
r <- as.data.frame(devtools::test(reporter = "summary"))
cat(sprintf("VERDICT: pass=%d fail=%d warn=%d\n", sum(r$passed, na.rm = TRUE),
            sum(r$failed, na.rm = TRUE), sum(r$warning, na.rm = TRUE)),
    file = out, append = TRUE)
```
Run:
```bash
Rscript diag_run_suite.R > diag_suite_console.txt 2>&1
```
Expected: `VERDICT: fail=0 ...` (pass count = prior 252 + the new isoplane/builder/renderer tests; pre-existing ~49 Nelder-Mead warnings from ternary fitting are non-fatal).

- [ ] **Step 2: Clean up scratch**

```bash
rm -f diag_run_suite.R diag_suite_out.txt diag_suite_console.txt
```

- [ ] **Step 3: Fast-forward main and push (only when asked)**

Per the project's merge cadence, fast-forward `main` to the branch tip in place (no checkout) and push:
```bash
git branch -f main mixdra-engine
git push origin main
```
(`R CMD check` remains a release-time step on an Rtools machine — not part of this plan.)

---

## Notes

- New exported API: `ec50_isoplane()`, `sigma_tu_curve()`, `ec50_markers()`, `plot_isoplane()`, `plot_sigma_tu()`. Internal: `simplex_grid()`, `.f4()`, `isoplane_plot_data()`, `sigma_tu_plot_data()`.
- All units are pure/plotly-free except the two renderers (guarded on `plotly`), matching the Plan-2 plotting layer.
- App wiring (T3), IA/quantal, and arbitrary effect levels remain out of scope (see spec §1).
```
