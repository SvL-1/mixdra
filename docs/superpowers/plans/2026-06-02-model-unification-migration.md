# Model-Unification Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the unified `mix_response()` the single production predictor for the 16 SA/DR/DL/reference mixture models, deleting ~1,592 lines of duplicated maths while preserving behaviour and the equivalence proof.

**Architecture:** `registry.R::model_spec()` is the only seam every consumer uses. We change *what* `$fn` resolves to: for `reference`/`SA`/`DR`/`DL` it returns an adapter closure that repacks the historical by-name call into `mix_response()`'s vector signature; for `ASA` it stays `get("ca_asa_tri_vec")`. The 16 verbatim functions move to a testthat helper so the equivalence test keeps cross-checking against them, then are deleted from the package.

**Tech Stack:** R package (`mixdra`), `testthat` edition 3, `devtools`/`roxygen2`. Tests run filtered via `devtools::test(filter=...)` (the full suite is deliberately avoided — it is slow).

**Spec:** `docs/superpowers/specs/2026-06-02-model-unification-migration-design.md`

> **✅ STATUS: COMPLETED (2026-06-02).** Implemented on branch `model-unification-migration`
> (7 commits `20be962`..`033cc9d`), merged fast-forward into `main`. Net **−1,539 lines** in `R/`.
> All 10 relevant test groups green on merged `main` (764 assertions, 0 failures), incl. the
> 576-assertion equivalence contract. One ~5e-8 tolerance re-baseline (`test-predict-mixture.R`).
> ASA left separate as planned. Each task passed spec + code-quality review plus a final
> holistic review ("ready to merge").

---

## Prerequisites (read before starting)

- **Working-tree state.** The current branch `binary-explain-fit` has uncommitted edits to `R/analyse.R`, `R/app-run.R`, `R/registry.R`, `tests/testthat/test-nchem.R` from in-progress work. **Commit (or stash) those first** so this migration's commits are isolated and reviewable. The `registry.R` edit already in the tree removed the ternary DR/DL dispatch entries — this plan builds on top of that state.
- **Test command idiom.** This repo runs filtered tests, e.g. `Rscript -e 'devtools::test(filter="mixture-predict")'` (the `filter` matches the test file name minus the `test-` prefix and `.R` suffix). If `devtools::load_all()` hits a libPaths/user-lib issue, set `R_LIBS_USER` to the project's library as the existing test setup does. Confirm the baseline is green before Task 1: `Rscript -e 'devtools::test(filter="mixture-predict")'` → all PASS.
- **Do not touch** `R/models-ternary-asa.R` or the `ASA` path anywhere. ASA is out of scope.

---

## File Structure

| File | Responsibility | Change |
| --- | --- | --- |
| `R/mix-response.R` | Unified predictor **+ the production adapter** (`make_adapter`, `.adapter_bmap`) | Modify: append adapter helpers |
| `R/registry.R` | `model_spec()` dispatch | Modify: return adapter for the 16, keep ASA via `get()` |
| `R/models-binary.R` | 8 verbatim binary fns + Vectorize wrappers | **Delete** (relocated to helper) |
| `R/models-ternary.R` | 8 verbatim ternary fns + Vectorize wrappers | **Delete** (relocated to helper) |
| `R/models-ternary-asa.R` | `ca_asa_tri` (ASA) | Unchanged |
| `tests/testthat/helper-legacy-models.R` | The 16 functions, as the retained equivalence oracle | **Create** |
| `tests/testthat/test-adapter.R` | Unit test for `make_adapter` glue (param unpack + b-mapping) | **Create** |
| `tests/testthat/test-model-spec-dispatch.R` | `model_spec(...)$fn` reproduces the legacy predictors | **Create** |
| `tests/testthat/test-mixture-predict.R` | Equivalence proof | Unchanged (now resolves the oracle from the helper) |
| `man/ca_*.Rd`, `man/ia_*.Rd` (16 files) | Auto-generated docs for the deleted internal fns | **Delete** via `devtools::document()` |

---

## Task 1: Adapter closure that wraps `mix_response()`

Add the production adapter to `R/mix-response.R` (co-located with the predictor it wraps, keeping `registry.R` minimal). Unit-test it directly against the still-present legacy `*_vec` functions.

**Files:**
- Modify: `R/mix-response.R` (append after the `mix_response` definition)
- Test: `tests/testthat/test-adapter.R` (create)

- [x] **Step 1: Write the failing test**

Create `tests/testthat/test-adapter.R`:

```r
# make_adapter() must reproduce the legacy *_vec predictors exactly (same maths,
# same by-name call contract). The legacy functions are still in R/ at this point.

test_that("make_adapter reproduces binary CA reference/SA over a grid", {
  m <- 800; sl <- c(6, 0.4); e <- c(0.08, 50); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))

  fn_ref <- make_adapter("CA", "reference", 2)
  got <- fn_ref(c1 = g$c1, c2 = g$c2, max = m,
                slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2])
  want <- ca_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2])
  expect_equal(unname(got), unname(want), tolerance = 1e-4)

  fn_sa <- make_adapter("CA", "SA", 2)
  got <- fn_sa(c1 = g$c1, c2 = g$c2, max = m,
               slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2], a = a)
  want <- ca_sa_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a)
  expect_equal(unname(got), unname(want), tolerance = 1e-4)
})

test_that("make_adapter maps binary DR b -> c(b, 0) and DL b scalar", {
  m <- 800; sl <- c(6, 0.4); e <- c(0.08, 50); a <- 0.5; b <- 0.3
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))

  fn_dr <- make_adapter("IA", "DR", 2)
  got <- fn_dr(c1 = g$c1, c2 = g$c2, max = m,
               slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2], a = a, b = b)
  want <- ia_dr_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, b)
  expect_equal(unname(got), unname(want), tolerance = 1e-8)

  fn_dl <- make_adapter("IA", "DL", 2)
  got <- fn_dl(c1 = g$c1, c2 = g$c2, max = m,
               slope1 = sl[1], slope2 = sl[2], ec501 = e[1], ec502 = e[2], a = a, b = 0.2)
  want <- ia_dl_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, 0.2)
  expect_equal(unname(got), unname(want), tolerance = 1e-8)
})

test_that("make_adapter reproduces ternary reference/SA", {
  m <- 800; sl <- c(6, 0.4, 2); e <- c(0.08, 50, 10); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60), c3 = c(0, 5, 30))

  fn <- make_adapter("CA", "SA", 3)
  got <- fn(c1 = g$c1, c2 = g$c2, c3 = g$c3, max = m,
            slope1 = sl[1], slope2 = sl[2], slope3 = sl[3],
            ec50_1 = e[1], ec50_2 = e[2], ec50_3 = e[3], a = a)
  want <- ca_sa_tri_vec(g$c1, g$c2, g$c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3], a)
  expect_equal(unname(got), unname(want), tolerance = 1e-4)
})
```

- [x] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'devtools::test(filter="adapter")'`
Expected: FAIL — `could not find function "make_adapter"`.

- [x] **Step 3: Implement `make_adapter` + `.adapter_bmap`**

Append to `R/mix-response.R`:

```r
# --- production dispatch adapter ----------------------------------------------
# model_spec() returns make_adapter(...) as `$fn` for the reference/SA/DR/DL
# family. It preserves the historical by-name call contract used by
# `.mixture_eval()` / `param_ci()` -- do.call(fn, c(conc_list, as.list(par))) --
# and repacks each measurement row into mix_response()'s vector signature.

# Map the flat `b` parameter onto mix_response()'s expected shape.
#   DR binary : scalar b -> c(b, 0)   (b weights chemical 1's dose fraction)
#   DR ternary: b1,b2,b3 -> c(b1,b2,b3)   (not reached via model_spec today, but
#               kept faithful to mix_response()'s contract / the equivalence test)
#   DL        : scalar b passes through
#   reference / SA: b unused -> 0
.adapter_bmap <- function(deviation, n_chem, A) {
  if (deviation == "DR") {
    if (n_chem == 2) c(A$b, 0) else c(A$b1, A$b2, A$b3)
  } else if (deviation == "DL") {
    A$b
  } else {
    0
  }
}

#' Build a vectorised, by-name predictor over `mix_response()`
#'
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param n_chem 2 or 3.
#' @return A function accepting named concentration vectors (`c1`, `c2`[, `c3`])
#'   and named scalar parameters (`max`, `slope1..`, `ec50..`, `a`, `b`/`b1..b3`),
#'   returning one prediction per row.
#' @keywords internal
make_adapter <- function(reference, deviation, n_chem) {
  slope_names <- if (n_chem == 2) c("slope1", "slope2")
                 else            c("slope1", "slope2", "slope3")
  ec50_names  <- if (n_chem == 2) c("ec501", "ec502")
                 else            c("ec50_1", "ec50_2", "ec50_3")
  cc_names    <- paste0("c", seq_len(n_chem))

  function(...) {
    A      <- list(...)
    concs  <- lapply(cc_names, function(k) A[[k]])
    slopes <- vapply(slope_names, function(k) A[[k]], numeric(1))
    ec50s  <- vapply(ec50_names,  function(k) A[[k]], numeric(1))
    a      <- if (is.null(A$a)) 0 else A$a
    b      <- .adapter_bmap(deviation, n_chem, A)
    n_pts  <- length(concs[[1]])
    vapply(seq_len(n_pts), function(i) {
      mix_response(vapply(concs, function(v) v[[i]], numeric(1)),
                   A$max, slopes, ec50s, reference, deviation, a = a, b = b)
    }, numeric(1))
  }
}
```

- [x] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'devtools::test(filter="adapter")'`
Expected: PASS (3 tests).

- [x] **Step 5: Commit**

```bash
git add R/mix-response.R tests/testthat/test-adapter.R
git commit -m "feat: make_adapter() wraps mix_response() for production dispatch"
```

---

## Task 2: Wire `model_spec()` to the adapter

Switch `model_spec()` so the 16 models resolve through `make_adapter`; ASA keeps `get()`. Add a dispatch test proving `model_spec(...)$fn` reproduces the legacy predictors, and confirm the full equivalence suite still passes.

**Files:**
- Modify: `R/registry.R` (the `dev_key` / `key` / `list(fn = ...)` block)
- Test: `tests/testthat/test-model-spec-dispatch.R` (create)

- [x] **Step 1: Write the failing test**

Create `tests/testthat/test-model-spec-dispatch.R`:

```r
# model_spec(...)$fn is the single production seam. It must reproduce the legacy
# *_vec predictors for every dispatched (reference, deviation, n_chem) combo and
# leave ASA resolving to ca_asa_tri_vec.

test_that("binary dispatch reproduces legacy predictors", {
  m <- 800; sl <- c(6, 0.4); e <- c(0.08, 50); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))
  call_fn <- function(fn, extra = list())
    do.call(fn, c(list(c1 = g$c1, c2 = g$c2, max = m,
                       slope1 = sl[1], slope2 = sl[2],
                       ec501 = e[1], ec502 = e[2]), extra))

  expect_equal(unname(call_fn(model_spec("CA", "reference", 2)$fn)),
               unname(ca_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2])),
               tolerance = 1e-4)
  expect_equal(unname(call_fn(model_spec("CA", "DR", 2)$fn, list(a = a, b = 0.3))),
               unname(ca_dr_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, 0.3)),
               tolerance = 1e-4)
  expect_equal(unname(call_fn(model_spec("IA", "DL", 2)$fn, list(a = a, b = 0.2))),
               unname(ia_dl_bi_vec(g$c1, g$c2, m, sl[1], sl[2], e[1], e[2], a, 0.2)),
               tolerance = 1e-8)
})

test_that("ternary dispatch reproduces legacy reference/SA", {
  m <- 800; sl <- c(6, 0.4, 2); e <- c(0.08, 50, 10); a <- 0.5
  g <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60), c3 = c(0, 5, 30))
  fn <- model_spec("IA", "SA", 3)$fn
  got <- fn(c1 = g$c1, c2 = g$c2, c3 = g$c3, max = m,
            slope1 = sl[1], slope2 = sl[2], slope3 = sl[3],
            ec50_1 = e[1], ec50_2 = e[2], ec50_3 = e[3], a = a)
  expect_equal(unname(got),
               unname(ia_sa_tri_vec(g$c1, g$c2, g$c3, m, sl[1], sl[2], sl[3],
                                    e[1], e[2], e[3], a)),
               tolerance = 1e-8)
})

test_that("ASA still resolves to the dedicated ternary predictor", {
  fn <- model_spec("CA", "ASA", 3)$fn
  expect_identical(fn, ca_asa_tri_vec)
})

test_that("dispatched params are unchanged", {
  expect_equal(model_spec("CA", "DR", 2)$params,
               c("max", "slope1", "slope2", "ec501", "ec502", "a", "b"))
  expect_equal(model_spec("CA", "ASA", 3)$extra, c("A1", "A2", "A3", "A4"))
})
```

- [x] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'devtools::test(filter="model-spec-dispatch")'`
Expected: FAIL — the binary/ternary `expect_identical`/`expect_equal` for the adapter path fail because `$fn` is still the legacy `get()` result (the ASA and params tests pass, the reproduce tests fail only if numerics differ — if they already pass, proceed; the meaningful assertion is the post-change green run plus `expect_identical(fn, ca_asa_tri_vec)` which must stay true).

> Note: the "reproduce" assertions may already pass against the legacy `$fn` (it *is* the legacy fn). That is fine — this test's job is to lock the contract so Step 3's swap cannot regress it. Treat a green Step 2 as an acceptable starting state; the gate is Step 4.

- [x] **Step 3: Swap the dispatch in `model_spec()`**

In `R/registry.R`, replace the `dev_key` / `key` / `list(...)` tail of `model_spec()` (currently building a `_vec` key and calling `get(key, mode = "function")` for every deviation) with ASA-only `get()` plus the adapter for everything else:

```r
  parent <- switch(deviation,
                   reference = NULL,
                   SA = "reference",
                   DR = "SA",
                   DL = "SA",
                   ASA = "SA")

  fn <- if (deviation == "ASA") {
    get(paste0(tolower(reference), "_asa_", suffix, "_vec"), mode = "function")
  } else {
    make_adapter(reference, deviation, n_chem)
  }

  list(fn = fn,
       params = c(base_params, extra),
       extra = extra,
       parent = parent)
```

Delete the now-unused `dev_key <- switch(...)` and `key <- paste0(...)` lines.

- [x] **Step 4: Run dispatch + equivalence tests**

Run: `Rscript -e 'devtools::test(filter="model-spec-dispatch")'`
Expected: PASS (4 tests), including `expect_identical(fn, ca_asa_tri_vec)`.

Run: `Rscript -e 'devtools::test(filter="mixture-predict")'`
Expected: PASS (the equivalence test is unaffected — it calls `mix_response()` directly).

Run: `Rscript -e 'devtools::test(filter="registry")'`
Expected: PASS (existing registry tests; if any asserted the old `$fn` identity against a `*_vec`, update it to the new contract — adapter for the 16, `ca_asa_tri_vec` for ASA).

- [x] **Step 5: Commit**

```bash
git add R/registry.R tests/testthat/test-model-spec-dispatch.R
git commit -m "feat: model_spec() dispatches the 16 models through mix_response()"
```

---

## Task 3: Relocate the 16 functions to a test helper, delete the package copies

Production no longer references the legacy functions — only the equivalence and `*-models-*` tests do. Move them verbatim into a testthat helper (auto-sourced before tests) and delete `R/models-binary.R` / `R/models-ternary.R`.

**Files:**
- Create: `tests/testthat/helper-legacy-models.R`
- Delete: `R/models-binary.R`, `R/models-ternary.R`

- [x] **Step 1: Create the helper from the existing sources**

```bash
{
  echo "# Legacy verbatim mixture predictors (CA/IA x reference/SA/DR/DL x binary/ternary)."
  echo "# RELOCATED from R/models-binary.R + R/models-ternary.R. These are NO LONGER"
  echo "# part of the package: production resolves mix_response() via model_spec()."
  echo "# They are retained here solely as the equivalence oracle for"
  echo "# test-mixture-predict.R (and the subjects of test-models-binary/ternary.R)."
  echo "# testthat auto-sources helper-*.R before running tests."
  echo ""
  cat R/models-binary.R
  echo ""
  cat R/models-ternary.R
} > tests/testthat/helper-legacy-models.R
```

(The roxygen `#'` lines come along harmlessly as comments; roxygen only scans `R/`, so they generate no `.Rd`.)

- [x] **Step 2: Delete the package copies**

```bash
git rm R/models-binary.R R/models-ternary.R
```

- [x] **Step 3: Run the oracle-dependent tests**

Run: `Rscript -e 'devtools::test(filter="mixture-predict")'`
Expected: PASS — `mix_response()` (package) vs the legacy functions (helper) still 576/576.

Run: `Rscript -e 'devtools::test(filter="models-binary")'`
Run: `Rscript -e 'devtools::test(filter="models-ternary")'`
Expected: PASS — these tests call the legacy functions, now resolved from the helper.

Run: `Rscript -e 'devtools::test(filter="adapter")'`
Run: `Rscript -e 'devtools::test(filter="model-spec-dispatch")'`
Expected: PASS — the adapter/dispatch tests also reference the legacy `*_vec` as oracle (now from the helper).

- [x] **Step 4: Verify no package code still references the legacy names**

Run: `Rscript -e 'devtools::load_all(); cat("loaded OK\n")'`
Expected: `loaded OK` with no "object not found" errors.

Run (must return nothing): `grep -rnE "\b(ca|ia)_(bi|tri|sa_bi|dr_bi|dl_bi|sa_tri|dr_tri|dl_tri)(_vec)?\b" R/`
Expected: empty output (only `ca_asa_tri`/`ca_asa_tri_vec` remain in `R/`, which this pattern excludes).

- [x] **Step 5: Commit**

```bash
git add tests/testthat/helper-legacy-models.R
git commit -m "refactor: relocate the 16 legacy predictors to a test-only oracle helper"
```

---

## Task 4: Regenerate documentation

The 16 deleted functions had `@keywords internal` `.Rd` pages. Regenerate so `man/` and `NAMESPACE` reflect reality.

**Files:**
- Delete (auto): `man/ca_bi.Rd`, `man/ia_bi.Rd`, `man/ca_sa_bi.Rd`, `man/ca_dr_bi.Rd`, `man/ca_dl_bi.Rd`, `man/ia_sa_bi.Rd`, `man/ia_dr_bi.Rd`, `man/ia_dl_bi.Rd`, `man/ca_tri.Rd`, `man/ia_tri.Rd`, `man/ca_sa_tri.Rd`, `man/ca_dr_tri.Rd`, `man/ca_dl_tri.Rd`, `man/ia_sa_tri.Rd`, `man/ia_dr_tri.Rd`, `man/ia_dl_tri.Rd`
- Possibly create: `man/make_adapter.Rd` (from the new roxygen block)

- [x] **Step 1: Regenerate**

Run: `Rscript -e 'devtools::document()'`
Expected: roxygen removes the 16 internal `.Rd` files and writes `man/make_adapter.Rd`. `NAMESPACE` is unchanged (none of these were exported — confirm the diff shows no `NAMESPACE` change).

- [x] **Step 2: Confirm the man/ delta is exactly the 16 removals + make_adapter**

Run: `git status --short man/`
Expected: 16 `D man/<legacy>.Rd` and one `?? man/make_adapter.Rd` (or `A` after add). No other `.Rd` touched.

- [x] **Step 3: Stage and commit**

```bash
git add -A man/ NAMESPACE
git commit -m "docs: regenerate man/ after removing the 16 legacy model functions"
```

---

## Task 5: Re-baseline numerics and broad verification

`mix_response()`'s CA bisection (200 iters / `1e-10`) is tighter than the legacy 50-iter / `1e-6` loop, so optim-derived fits and reported SSR/R²/F/P shift in the last digits. Run the broader filtered groups, confirm any failure is a sub-`1e-3` relative shift (not a logic break), and re-baseline the affected golden values.

**Files:**
- Modify (only where a numeric expectation shifts): likely `tests/testthat/test-fit.R`, `test-staged.R`, `test-validation-binary.R`, `test-validation-ternary.R`, `test-recovery.R`, `test-summary.R`, `test-fit-bounds.R`, `test-predict-mixture.R`

- [x] **Step 1: Run each broad group and record failures**

Run, one at a time:
```
Rscript -e 'devtools::test(filter="fit")'
Rscript -e 'devtools::test(filter="staged")'
Rscript -e 'devtools::test(filter="validation-binary")'
Rscript -e 'devtools::test(filter="validation-ternary")'
Rscript -e 'devtools::test(filter="recovery")'
Rscript -e 'devtools::test(filter="summary")'
Rscript -e 'devtools::test(filter="predict-mixture")'
Rscript -e 'devtools::test(filter="fit-bounds")'
```
Expected: mostly PASS. For any FAIL, note the file, the expectation, the expected-vs-actual numbers.

- [x] **Step 2: Triage each failure with the decision rule**

For every failing expectation:
- Compute the relative difference `abs(actual - expected) / max(abs(expected), 1e-8)`.
- **If `< 1e-3`** (a last-digit convergence shift): this is the expected re-baseline. Update the literal expected value in the test to `actual`, and/or widen an over-tight `tolerance =` to a documented value (e.g. `1e-3`). Add a one-line comment: `# re-baselined for mix_response() 200-iter bisection (was legacy 50-iter)`.
- **If `>= 1e-3` or a structural/logic error** (NaN, wrong length, wrong sign): STOP — this is not a re-baseline. Investigate via the systematic-debugging skill; do not paper over it by loosening tolerance.

- [x] **Step 3: Re-run the edited groups to green**

Re-run each `devtools::test(filter=...)` you edited.
Expected: PASS.

- [x] **Step 4: Smoke-check dispatch across every reachable combo**

Run:
```bash
Rscript -e '
devtools::load_all()
combos <- list(c("CA","reference"), c("CA","SA"), c("CA","DR"), c("CA","DL"),
               c("IA","reference"), c("IA","SA"), c("IA","DR"), c("IA","DL"))
for (cb in combos) {
  fn <- model_spec(cb[1], cb[2], 2)$fn
  y <- fn(c1=c(0,0.1,0.3), c2=c(0,5,20), max=800, slope1=3, slope2=2,
          ec501=0.1, ec502=10, a=0.4, b=0.2)
  stopifnot(all(is.finite(y)), length(y)==3)
}
for (cb in list(c("CA","reference"), c("CA","SA"), c("IA","reference"), c("IA","SA"))) {
  fn <- model_spec(cb[1], cb[2], 3)$fn
  y <- fn(c1=c(0,0.1), c2=c(0,5), c3=c(0,2), max=800, slope1=3, slope2=2, slope3=1.5,
          ec50_1=0.1, ec50_2=10, ec50_3=4, a=0.4)
  stopifnot(all(is.finite(y)), length(y)==2)
}
fn <- model_spec("CA","ASA",3)$fn
y <- fn(c1=c(0,0.1), c2=c(0,5), c3=c(0,2), max=800, slope1=3, slope2=2, slope3=1.5,
        ec50_1=0.1, ec50_2=10, ec50_3=4, A1=0.2, A2=0.1, A3=0.1, A4=0.05)
stopifnot(all(is.finite(y)), length(y)==2)
cat("all dispatch combos finite OK\n")
'
```
Expected: `all dispatch combos finite OK`.

- [x] **Step 5: Commit**

```bash
git add tests/
git commit -m "test: re-baseline numerics for the unified mix_response() bisection"
```

---

## Task 6: Final sweep and memory update

- [x] **Step 1: Run the app-layer and remaining model tests to catch indirect breakage**

```
Rscript -e 'devtools::test(filter="app-modules")'
Rscript -e 'devtools::test(filter="plot-data")'
Rscript -e 'devtools::test(filter="simulate")'
Rscript -e 'devtools::test(filter="compare")'
Rscript -e 'devtools::test(filter="models-ternary-asa")'
```
Expected: PASS. (These exercise `model_spec` indirectly via fitting/plotting/simulation.)

- [x] **Step 2: Confirm the line-count win**

Run: `git diff --stat main -- R/`
Expected: `R/models-binary.R` (−408) and `R/models-ternary.R` (−1184) deleted; `R/mix-response.R` and `R/registry.R` net small additions. Net ≈ −1,500 lines in `R/`.

- [x] **Step 3: Update the project memory**

Edit `C:\Users\jelle\.claude\projects\D--sam\memory\mixdra-model-unification.md` to record that the migration is complete: production dispatches through `mix_response()` via a `make_adapter()` closure in `model_spec()`; the 16 legacy functions now live only in `tests/testthat/helper-legacy-models.R` as the equivalence oracle; ASA remains separate. Keep the `MEMORY.md` pointer line in sync.

- [x] **Step 4: Final commit**

```bash
git add C:/Users/jelle/.claude/projects/D--sam/memory/
git commit -m "docs: record completed mix_response() production migration in memory"
```

---

## Self-Review notes (for the executor)

- **Equivalence preserved:** `test-mixture-predict.R` is untouched and keeps comparing `mix_response()` (package) against the legacy functions (helper). If it goes red at any task, the adapter/bisection — not the test — is wrong.
- **The one fiddly transform** is `.adapter_bmap` (binary DR `b -> c(b, 0)`). Task 1's DR test and Task 2's DR dispatch test both pin it.
- **Do not loosen `mixture-predict` tolerances** to make Task 5 pass — those tolerances are the equivalence contract. Only golden-value fit tests get re-baselined.
- **ASA is never touched.** `model_spec("CA","ASA",3)$fn` must stay `expect_identical` to `ca_asa_tri_vec` (Task 2, Step 1).
