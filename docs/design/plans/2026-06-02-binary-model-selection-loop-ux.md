# Binary Model-Selection Joint-Fit Loop (UX Spine) Implementation Plan

> **STATUS (2026-06-03): executed on branch `binary-joint-selection-loop`, with a post-implementation reversal.** Tasks 1–5 were implemented joint as written, then the loop was switched back to the **staged** fit (Sam's required method — joint masks interactions; see the spec's AMENDMENT and `mixdra-staged-fitting`). The joint engine fns (`joint_fit_one`, `analyse_mixture_joint`) were removed; `compare_fits` + `selection_chain_order` remain. Where tasks below say "joint", the shipped loop calls `fit_model(..., start = curve_params(), fixed = names(curve_params()))` per model instead. The UX, live stepper, and LR comparison are unchanged.

**Goal:** Make the binary tab's automatic interaction-model choice the spine of the UI: a single "compare all models" action runs a warm-started, fully-joint fit per model (reference → SA → DR/DL), fills a live results table row by row, and highlights the parsimonious winner — retiring the separate staged "Find best" and manual "Optimize all" steps.

**Architecture:** A small engine layer adds a per-model joint fit (`joint_fit_one`) and a full-chain convenience (`analyse_mixture_joint`), reusing the existing `fit_model` optimiser and the `lr_test`/`select_parsimonious` comparison (extracted into a regime-agnostic `compare_fits` helper). The Shiny module drives the chain step-by-step with a timer (`invalidateLater`) so each model's row paints as it finishes, then runs `compare_fits` once all four are stored. The UI collapses from 4 stages to 3 (seed curves → compare-all hero → inspect a model).

**Tech Stack:** R package (`mixdra`), `stats::optim` (L-BFGS-B + Nelder-Mead), Shiny + bslib + DT + plotly, testthat (`devtools::test`).

**Spec:** `docs/design/specs/2026-06-02-binary-model-selection-loop-ux-design.md`

---

## File Structure

- **`R/compare.R`** — add `compare_fits()` (builds the comparison frame + chosen model from a set of fits; regime-agnostic). Existing `lr_test`/`select_parsimonious` unchanged.
- **`R/analyse.R`** — add `joint_chain_order()`, `joint_fit_one()` (one warm-started joint model fit), `analyse_mixture_joint()` (full chain). Refactor `analyse_mixture()` to call `compare_fits()`.
- **`R/app-binary.R`** — restructure `binary_ui()` to 3 stages; replace the `find_best`/`optimize_all` server observers with the live stepped loop; demote manual Autofit/Simulate to an "explore" accordion.
- **`tests/testthat/test-compare.R`** — cover `compare_fits()`.
- **`tests/testthat/test-joint-loop.R`** (new) — cover `joint_fit_one()` monotonicity + `analyse_mixture_joint()` recovery/selection.
- **`tests/testthat/test-app-modules.R`** — rewrite the `binary_server` workspace test (loop instead of find_best/optimize_all) and the `binary_ui` HTML assertions for the 3-stage layout.

---

## Task 1: Extract `compare_fits()` and refactor `analyse_mixture()`

**Files:**
- Modify: `R/compare.R` (add `compare_fits()` after `select_parsimonious`, ~line 49)
- Modify: `R/analyse.R:136-148` (replace the inline comparison block)
- Test: `tests/testthat/test-compare.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-compare.R`:

```r
test_that("compare_fits builds the LR comparison frame and chooses a model", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2), C2 = c(0, 0.5, 5))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  seed <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  fits <- lapply(c("reference", "SA", "DR", "DL"), function(d)
    fit_model(df, "CA", d, "continuous", start = seed, fixed = names(seed)))
  names(fits) <- c("reference", "SA", "DR", "DL")

  cmp <- compare_fits(fits, n = nrow(df), response = "continuous", alpha = 0.05)
  expect_true(all(c("model", "parent", "chi", "df", "p") %in% names(cmp$comparison)))
  expect_equal(nrow(cmp$comparison), 3)          # SA, DR, DL each vs their parent
  expect_true(cmp$chosen %in% c("reference", "SA", "DR", "DL"))
})

test_that("compare_fits handles a ternary (reference + SA only) fit set", {
  fits <- list(
    reference = list(objective = 10, df = 7),
    SA        = list(objective = 9,  df = 8))
  cmp <- compare_fits(fits, n = 27, response = "continuous", alpha = 0.05)
  expect_equal(nrow(cmp$comparison), 1)          # only SA vs reference
  expect_equal(cmp$comparison$model, "SA")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter='compare')"`
Expected: FAIL — `could not find function "compare_fits"`.

- [ ] **Step 3: Add `compare_fits()` to `R/compare.R`**

Append after `select_parsimonious()` (after line 49):

```r
#' Build the nested-model comparison frame and choose a model
#'
#' Regime-agnostic: works whether the fits are staged (curves fixed, only a/b
#' free) or fully joint (all parameters free). It only reads each fit's
#' `objective` and `df`, so the likelihood-ratio differences are identical in
#' form either way. The parent chain is `SA` vs `reference`, `DR`/`DL` vs `SA`;
#' models absent from `fits` are dropped (e.g. ternary reference + SA only).
#' @param fits Named list of fits (subset of `reference`, `SA`, `DR`, `DL`),
#'   each with `objective` and `df`.
#' @param n Number of observations.
#' @param response "continuous" or "binary".
#' @param alpha Significance threshold for [select_parsimonious()].
#' @return A list: `comparison` (data frame of LR tests, or `NULL` if no
#'   non-reference model is present) and `chosen` (selected model name).
#' @export
compare_fits <- function(fits, n, response, alpha = 0.05) {
  parent_of <- c(SA = "reference", DR = "SA", DL = "SA")
  parent_of <- parent_of[intersect(names(parent_of), names(fits))]
  comparison <- if (length(parent_of)) {
    do.call(rbind, lapply(names(parent_of), function(m) {
      p <- parent_of[[m]]
      t <- lr_test(fits[[p]]$objective, fits[[m]]$objective,
                   fits[[p]]$df, fits[[m]]$df, n, response)
      data.frame(model = m, parent = p, chi = t$chi, df = t$df, p = t$p)
    }))
  } else {
    NULL
  }
  list(comparison = comparison,
       chosen = select_parsimonious(fits, n, response, alpha))
}
```

- [ ] **Step 4: Refactor `analyse_mixture()` to use it**

In `R/analyse.R`, replace lines 136-148 (from `n <- nrow(df)` through the closing `list(...)`) with:

```r
  cmp <- compare_fits(fits, nrow(df), response, alpha)
  list(fits = fits, comparison = cmp$comparison, chosen = cmp$chosen,
       reference = reference, response = response)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter='compare')"`
Run: `Rscript -e "devtools::test(filter='analyse')"`
Expected: PASS for both (the existing `analyse` test still sees the same `$comparison`/`$chosen` shape).

- [ ] **Step 6: Commit**

```bash
git add R/compare.R R/analyse.R tests/testthat/test-compare.R
git commit -m "refactor(engine): extract regime-agnostic compare_fits() helper"
```

---

## Task 2: `joint_fit_one()` — one warm-started joint model fit

**Files:**
- Modify: `R/analyse.R` (add `joint_chain_order()` + `joint_fit_one()` after `fit_curve_from_singles()`, ~line 75)
- Test: `tests/testthat/test-joint-loop.R` (new)

**Background:** A joint fit of one model is just `fit_model()` with the curve params left *free* and seeded from a starting vector. Warm-starting each child from its parent's joint fit — and seeding the child's NEW interaction parameter at its neutral value (`0`, which reduces SA→reference, and DR/DL→SA) — guarantees the child can do no worse than the parent, so `child$objective <= parent$objective` and the LR statistic stays >= 0.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-joint-loop.R`:

```r
# Joint model-selection loop: every model is fully optimised (curves free),
# warm-started down the nesting chain. See
# docs/design/specs/2026-06-02-binary-model-selection-loop-ux-design.md

test_that("joint_chain_order is reference->SA->DR->DL for binary, reference->SA for ternary", {
  expect_equal(joint_chain_order(2), c("reference", "SA", "DR", "DL"))
  expect_equal(joint_chain_order(3), c("reference", "SA"))
})

test_that("joint_fit_one frees the curves and (warm-started) never worsens the parent", {
  df <- (function() {
    g <- expand.grid(C1 = c(0, 0.05, 0.2, 0.8), C2 = c(0, 0.5, 5, 20))
    g$Res <- ca_bi_vec(g$C1, g$C2, 800, 4, 1.5, 0.08, 1)
    g
  })()
  seed <- c(max = 700, slope1 = 3, slope2 = 2, ec501 = 0.1, ec502 = 1.5)  # off-true seed

  ref <- joint_fit_one(df, "CA", "reference", "continuous", seed)
  expect_equal(ref$df, 5L)                       # all five curve params were FREE
  expect_false(isTRUE(ref$simulated))

  sa <- joint_fit_one(df, "CA", "SA", "continuous", seed, parent_fit = ref)
  expect_equal(sa$df, 6L)                         # + a
  expect_true("a" %in% names(sa$par))
  expect_lte(sa$objective, ref$objective + 1e-6)  # warm-start monotonicity

  dr <- joint_fit_one(df, "CA", "DR", "continuous", seed, parent_fit = sa)
  expect_equal(dr$df, 7L)                         # + b
  expect_lte(dr$objective, sa$objective + 1e-6)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter='joint-loop')"`
Expected: FAIL — `could not find function "joint_chain_order"`.

- [ ] **Step 3: Implement in `R/analyse.R`**

Insert after `fit_curve_from_singles()` (after line 75):

```r
#' Model order for the joint selection chain
#'
#' Binary mixtures walk reference -> SA -> {DR, DL}; ternary mixtures support
#' only reference -> SA here (Advanced S/A is fitted separately).
#' @keywords internal
joint_chain_order <- function(n_chem) {
  if (n_chem == 2) c("reference", "SA", "DR", "DL") else c("reference", "SA")
}

#' Fit one mixture model jointly (curves + interaction), warm-started
#'
#' Unlike the staged [fit_curve_from_singles()] + fixed-curve fits used by
#' [analyse_mixture()], this leaves EVERY parameter free and optimises them
#' together (the Excel-faithful joint fit). It is seeded from `parent_fit` when
#' given (the parent in the nesting chain), else from `seed_curves`. The child's
#' new interaction parameter (the one its parent lacks) is seeded at its neutral
#' value `0` -- `a = 0` reduces SA to the reference, `b = 0` reduces DR/DL to SA
#' -- so the child starts at a point reproducing the parent's fit and can only
#' improve on it (keeping the nested LR comparison sound).
#' @param df Mixture data frame (see [fit_model()]).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param seed_curves Named numeric of curve params (max, slope*, ec50*) used as
#'   the start when there is no parent (the reference fit).
#' @param parent_fit The parent model's joint fit, or `NULL` for the reference.
#' @param n_starts,time_limit Forwarded to [fit_model()].
#' @return An enriched mixture fit (as [fit_model()]).
#' @keywords internal
joint_fit_one <- function(df, reference, deviation, response,
                          seed_curves, parent_fit = NULL,
                          n_starts = 1, time_limit = 30) {
  n_chem <- length(intersect(c("C1", "C2", "C3"), names(df)))
  spec   <- model_spec(reference, deviation, n_chem)
  start  <- if (is.null(parent_fit)) seed_curves else parent_fit$par
  # Neutral-seed any interaction parameter the child adds beyond its start
  # (a or b): 0 makes the new term vanish, reproducing the parent.
  new_extra <- setdiff(spec$extra, names(start))
  if (length(new_extra))
    start <- c(start, stats::setNames(rep(0, length(new_extra)), new_extra))
  fit_model(df, reference, deviation, response, start = start,
            n_starts = n_starts, time_limit = time_limit)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::test(filter='joint-loop')"`
Expected: PASS (all three `test_that` blocks).

- [ ] **Step 5: Commit**

```bash
git add R/analyse.R tests/testthat/test-joint-loop.R
git commit -m "feat(engine): joint_fit_one() warm-started per-model joint fit"
```

---

## Task 3: `analyse_mixture_joint()` — full chain + comparison

**Files:**
- Modify: `R/analyse.R` (add `analyse_mixture_joint()` after `joint_fit_one()`)
- Test: `tests/testthat/test-joint-loop.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-joint-loop.R`:

```r
test_that("analyse_mixture_joint fits all four jointly and recovers a reference truth", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  df  <- simulate_mixture(par, "CA", "reference", "continuous")
  res <- analyse_mixture_joint(df, reference = "CA", response = "continuous",
                               n_starts = 1)
  expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
  # joint fits free the curves, so df are the FULL counts (5/6/7), not 0/1/2
  expect_equal(res$fits$reference$df, 5L)
  expect_equal(res$fits$SA$df, 6L)
  expect_equal(res$fits$DR$df, 7L)
  # curves recovered jointly from the full data
  expect_equal(unname(res$fits$reference$par[["max"]]), 800, tolerance = 1e-2)
  expect_equal(unname(res$fits$reference$par[["ec501"]]), 0.08, tolerance = 1e-2)
  # noiseless reference truth -> no invented interaction
  expect_equal(unname(res$fits$SA$par[["a"]]), 0, tolerance = 1e-2)
  expect_equal(nrow(res$comparison), 3)
  expect_true(res$chosen %in% c("reference", "SA", "DR", "DL"))
})

test_that("analyse_mixture_joint selects DR on DR-generated data", {
  skip_on_cran()
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1,
           a = 2.5, b = 1.5)
  df  <- simulate_mixture(par, "CA", "DR", "continuous")
  set.seed(1)
  res <- analyse_mixture_joint(df, reference = "CA", response = "continuous",
                               n_starts = 10)
  expect_equal(unname(res$fits$DR$par[["a"]]), 2.5, tolerance = 0.1)
  expect_equal(unname(res$fits$DR$par[["b"]]), 1.5, tolerance = 0.1)
  expect_equal(res$chosen, "DR")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter='joint-loop')"`
Expected: FAIL — `could not find function "analyse_mixture_joint"`.

- [ ] **Step 3: Implement in `R/analyse.R`**

Insert after `joint_fit_one()`:

```r
#' Analyse a mixture by fully-joint fits (model selection via the joint loop)
#'
#' The joint counterpart of [analyse_mixture()]. Seeds the curve parameters from
#' the single-compound data (or a user-supplied `start`), then fits each model in
#' the nesting chain JOINTLY (all parameters free) via [joint_fit_one()],
#' warm-starting each from its parent. Selection reuses [compare_fits()]. This is
#' the engine entry point behind the binary tab's "Fit & compare all models"
#' action; the Shiny layer drives the same chain step-by-step for live updates.
#' @inheritParams analyse_mixture
#' @return A list: `fits`, `comparison`, `chosen`, `reference`, `response`.
#' @export
analyse_mixture_joint <- function(df, reference,
                                  response = c("continuous", "binary"),
                                  start = NULL, alpha = 0.05,
                                  n_starts = 1, time_limit = 30) {
  response <- match.arg(response)
  n_chem <- length(intersect(c("C1", "C2", "C3"), names(df)))
  seed <- if (is.null(start)) {
    fit_curve_from_singles(df, reference, response,
                           n_starts = n_starts, time_limit = time_limit)
  } else {
    start
  }
  fits <- list()
  for (dev in joint_chain_order(n_chem)) {
    parent_key <- model_spec(reference, dev, n_chem)$parent
    parent <- if (is.null(parent_key)) NULL else fits[[parent_key]]
    fits[[dev]] <- joint_fit_one(df, reference, dev, response, seed,
                                 parent_fit = parent,
                                 n_starts = n_starts, time_limit = time_limit)
  }
  cmp <- compare_fits(fits, nrow(df), response, alpha)
  list(fits = fits, comparison = cmp$comparison, chosen = cmp$chosen,
       reference = reference, response = response)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::test(filter='joint-loop')"`
Expected: PASS. (The second block runs only under `devtools::test()`, which sets `NOT_CRAN`.)

- [ ] **Step 5: Update NAMESPACE/docs**

Run: `Rscript -e "devtools::document()"`
Expected: `NAMESPACE` gains `export(compare_fits)` and `export(analyse_mixture_joint)`; `.Rd` files written.

- [ ] **Step 6: Commit**

```bash
git add R/analyse.R NAMESPACE man/ tests/testthat/test-joint-loop.R
git commit -m "feat(engine): analyse_mixture_joint() full joint selection chain"
```

---

## Task 4: Binary module — live stepped loop server logic

**Files:**
- Modify: `R/app-binary.R` (server `binary_server`, lines ~239-591)
- Test: `tests/testthat/test-app-modules.R:40-115` (rewrite the workspace test)

This task changes only the **server**. The UI still has the old Stage 2/3/4 cards (Task 5 restructures them); we add the `compare_all` input id in Task 5. To let the test run now, the rewritten test sets `compare_all` via `session$setInputs` (testServer does not require the input to exist in the UI).

- [ ] **Step 1: Rewrite the failing test**

Replace the whole `test_that("binary_server workspace: ...")` block (`tests/testthat/test-app-modules.R:40-115`) with:

```r
test_that("binary_server: compare-all loop fills the store live and selects a model", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30, model = "reference",
                      file = list(datapath = csv, name = "binary.csv"))
    expect_length(errs(), 0)

    # supply both single curves deterministically via Simulate, distinct max
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)
    expect_true(isTRUE(frozen()))
    expect_null(current_fit())

    # Kick the joint compare-all loop. The first model fits synchronously when
    # the queue is armed; the rest advance on timer ticks.
    session$setInputs(compare_all = 1)
    for (i in 1:8) session$elapse(5)             # drain the stepper queue

    # all four joint fits stored, with FULL (joint) df counts
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    expect_equal(fits_store()[["reference"]]$df, 5L)
    expect_equal(fits_store()[["DR"]]$df, 7L)
    # warm-start monotonicity along the chain
    expect_lte(fits_store()[["SA"]]$objective,
               fits_store()[["reference"]]$objective + 1e-6)

    # comparison + chosen produced once the chain completed
    expect_setequal(names(last_compare()$fits), c("reference", "SA", "DR", "DL"))
    expect_equal(nrow(last_compare()$comparison), 3)
    chosen <- last_compare()$chosen
    expect_true(chosen %in% c("reference", "SA", "DR", "DL"))

    # picker follows the chosen model (testServer doesn't echo updateSelectInput)
    session$setInputs(model = chosen)
    expect_equal(current_fit()$deviation, chosen)

    # "explore by hand" still works: Autofit the displayed model (staged a/b)
    session$setInputs(model = "SA", autofit = 1)
    expect_equal(current_fit()$deviation, "SA")
    expect_true("a" %in% names(current_fit()$par))

    # editing a single curve keeps the workspace and re-evaluates stored fits
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_true(isTRUE(frozen()))
    expect_gt(length(fits_store()), 0)
    expect_null(last_compare())                  # comparison invalidated by curve change
  })
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: FAIL — `compare_all` does nothing; `last_compare()` stays `NULL` / store stays empty.

- [ ] **Step 3: Remove the obsolete observers and outputs**

In `R/app-binary.R`, delete:
- the `output$has_fit` block (lines 326-328),
- the `optimize_pre`/`optimize_post` reactiveVals (lines 318-319),
- the entire `observeEvent(input$find_best, ...)` block (lines 444-467),
- the entire `observeEvent(input$optimize_all, ...)` block (lines 481-507),
- the `output$optimize_readout` block (lines 516-524).

In the `observeEvent(input$autofit, ...)` and `observeEvent(input$simulate, ...)` blocks, remove the now-undefined `optimize_pre(NULL); optimize_post(NULL)` lines (lines 413 and 440). In the `observeEvent(input$model, ...)` block remove the same line (476).

- [ ] **Step 4: Add the loop reactives + stepper**

In `R/app-binary.R`, immediately after the `fits_store`/`last_compare` declarations (after line 316), add:

```r
    # Number of chemicals present (binary tab -> 2; future-proofs the loop order).
    n_chem <- shiny::reactive(
      length(intersect(c("C1", "C2", "C3"), names(engine_df()))))

    # Live joint compare-all loop. `loop_queue` holds the models still to fit;
    # `stepper_on` arms the timer-driven stepper. We fit ONE model per tick and
    # let Shiny flush (paint the new table row) between ticks via invalidateLater.
    loop_queue <- shiny::reactiveVal(NULL)
    stepper_on <- shiny::reactiveVal(FALSE)
```

After the `observeEvent(list(fit1(), fit2()), ...)` re-evaluation block (after line 385), add the kick + stepper:

```r
    # Kick: clear the store and arm the joint chain (reference -> SA -> DR -> DL).
    shiny::observeEvent(input$compare_all, {
      shiny::req(frozen(), curve_params())
      fits_store(list())
      last_compare(NULL)
      loop_queue(joint_chain_order(n_chem()))
      stepper_on(TRUE)
    })

    # Stepper: re-runs on each timer tick while armed. Reads the queue with
    # isolate() so only the timer (not its own writes) re-triggers it, which is
    # what forces a client paint between models.
    shiny::observe({
      if (!isTRUE(stepper_on())) return()
      shiny::invalidateLater(0)                 # schedule the next tick
      q <- shiny::isolate(loop_queue())
      if (is.null(q) || length(q) == 0) { stepper_on(FALSE); return() }
      shiny::isolate({
        dev    <- q[[1]]
        s      <- fits_store()
        pkey   <- model_spec(input$reference, dev, n_chem())$parent
        parent <- if (is.null(pkey)) NULL else s[[pkey]]
        fit <- tryCatch(
          joint_fit_one(engine_df(), input$reference, dev, engine_response(),
                        seed_curves = curve_params(), parent_fit = parent,
                        n_starts = n_starts_eff(), time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(
              paste0("Fit failed (", dev, "): ", conditionMessage(e)), type = "error")
            NULL
          })
        if (!is.null(fit)) { s[[dev]] <- fit; fits_store(s) }
        rest <- q[-1]
        loop_queue(rest)
        if (length(rest) == 0) {
          stepper_on(FALSE)
          cmp <- compare_fits(fits_store(), nrow(engine_df()),
                              engine_response(), input$alpha)
          last_compare(list(fits = fits_store(), comparison = cmp$comparison,
                            chosen = cmp$chosen, reference = input$reference,
                            response = engine_response()))
          shiny::updateSelectInput(session, "model", selected = cmp$chosen)
          ch <- fits_store()[[cmp$chosen]]
          shiny::updateNumericInput(session, "val_a",
            value = if ("a" %in% names(ch$par)) round(ch$par[["a"]], 4) else NA)
          shiny::updateNumericInput(session, "val_b",
            value = if ("b" %in% names(ch$par)) round(ch$par[["b"]], 4) else NA)
        }
      })
    })
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: PASS for the rewritten `binary_server` loop test. (The `binary_ui` HTML tests in the same file will still FAIL — they assert the old layout; Task 5 fixes them.)

- [ ] **Step 6: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): live stepped joint compare-all loop (server)"
```

---

## Task 5: Binary module — 3-stage UI restructure

**Files:**
- Modify: `R/app-binary.R` (`binary_ui`, lines 77-233)
- Test: `tests/testthat/test-app-modules.R` (the `binary_ui` HTML-assertion blocks)

- [ ] **Step 1: Rewrite the failing UI tests**

Replace the `test_that("binary_ui: auto-reveal ...")` block (lines 137-157) and the `test_that("binary_ui: Optimize-all is a plain button ...")` block (lines 176-184) with:

```r
test_that("binary_ui: 3-stage layout with compare-all hero and explore accordion", {
  html <- as.character(binary_ui("binary"))
  # no manual Freeze step
  expect_false(grepl("Freeze curves", html, fixed = TRUE))
  # Stage 2 hero: the compare-all action + the promoted results table + alpha
  expect_match(html, "binary-compare_all", fixed = TRUE)
  expect_match(html, "Fit &amp; compare all models", fixed = TRUE)
  expect_match(html, "binary-results", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "significance threshold for the model comparison", fixed = TRUE)
  # Stage 3 inspect: model picker, per-model help, plots, CIs
  expect_match(html, "binary-model", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
  expect_match(html, "binary-surface", fixed = TRUE)
  expect_match(html, "binary-cis", fixed = TRUE)
  # "explore by hand" demoted: manual a/b value + Autofit/Simulate still present
  expect_match(html, "binary-val_a", fixed = TRUE)
  expect_match(html, "Autofit (a, b)", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  # the standalone Optimize-all stage is gone
  expect_false(grepl("binary-optimize_all", html, fixed = TRUE))
  expect_false(grepl("Optimize all params", html, fixed = TRUE))
  # the old buried staged Find-best button is gone (replaced by compare_all)
  expect_false(grepl("binary-find_best", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: FAIL — `binary-compare_all` not found in the old UI.

- [ ] **Step 3: Replace Stages 2-4 in `binary_ui`**

In `R/app-binary.R`, replace everything from the Stage 2 panel (line 127, the comment `# Stage 2 -- per-model interaction workspace ...`) through the end of the Stage 4 panel (line 231, the `)` closing the last `conditionalPanel`) with the three blocks below. Keep Stage 1 (lines 109-125) and the sidebar (lines 79-107) unchanged. Reword the Stage 1 intro `<p>` (lines 113-116) to say the curves *seed* the joint fits:

```r
      shiny::p("Fit each chemical's dose-response curve (Autofit or Simulate). ",
               "These seed the joint fits below — they are a starting point, ",
               "not frozen. The comparison workspace appears once both are fitted."),
```

Then the new Stages 2-3:

```r
    # Stage 2 -- the hero: fit & compare ALL interaction models jointly, live.
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Compare interaction models"),
        shiny::p("Fits every interaction model jointly (curves + interaction ",
                 "together), warm-started reference → S/A → DR/DL, and ",
                 "compares them. Rows fill in as each model finishes; the best ",
                 "model is highlighted."),
        bslib::layout_columns(
          col_widths = c(7, 5),
          shiny::div(class = "mt-2",
                     shiny::actionButton(ns("compare_all"),
                                         "Fit & compare all models",
                                         class = "btn-primary")),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05,
                              min = 0, max = 1, step = 0.01)),
        shiny::helpText(
          "alpha (α) is the significance threshold for the model comparison: ",
          "a more complex model is kept only if it improves the fit at p < α ",
          "(default 0.05). The table shows every fitted model and highlights the ",
          "selected (best) one."),
        DT::DTOutput(ns("results"))
      )
    ),

    # Stage 3 -- inspect the chosen (or any) model: diagnostics, CIs, and an
    # optional by-hand explore panel.
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 3 · Inspect a model"),
        shiny::selectInput(
          ns("model"), "Model",
          choices = c("No interaction (reference)" = "reference",
                      "Similar action (S/A)"       = "SA",
                      "Dose-ratio (DR)"            = "DR",
                      "Dose-level (DL)"            = "DL")),
        shiny::uiOutput(ns("interaction_help")),
        shiny::uiOutput(ns("objective")),
        bslib::layout_columns(
          bslib::card(bslib::card_header("3-D response surface"),
                      plotly::plotlyOutput(ns("surface"), height = "520px")),
          bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                      plotly::plotlyOutput(ns("isobole"), height = "520px"))
        ),
        bslib::card(bslib::card_header("Observed vs predicted"),
                    plotly::plotlyOutput(ns("op"))),
        bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                    DT::DTOutput(ns("cis"))),

        # Explore by hand: fit/enter the displayed model's a/b at the current
        # curves (staged) -- optional, off the main compare-all path.
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel(
            "Explore by hand",
            shiny::helpText("Fit or enter this model's interaction params at the ",
                            "current curves, without re-running the comparison."),
            shiny::conditionalPanel(
              condition = "input.model != 'reference'", ns = ns,
              bslib::card(
                bslib::card_header("Parameters"),
                shiny::fluidRow(
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
                  shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Value")))
                ),
                interaction_param_row(ns, "a", "a",
                  "Overall strength & direction (a > 0 antagonism, a < 0 synergism)."),
                shiny::conditionalPanel(
                  condition = "input.model == 'DR' || input.model == 'DL'", ns = ns,
                  interaction_param_row(ns, "b", "b",
                    "How the interaction shifts with the mixture ratio / dose level."))
              )
            ),
            shiny::div(
              shiny::actionButton(ns("autofit"), "Autofit (a, b)", class = "btn-primary"),
              shiny::actionButton(ns("simulate"), "Simulate")
            )
          )
        )
      )
    )
```

- [ ] **Step 4: Run the UI test to verify it passes**

Run: `Rscript -e "devtools::test(filter='app-modules')"`
Expected: PASS for the new `binary_ui` test and the `binary_server` loop test.

- [ ] **Step 5: Sanity-check the app boots**

Run: `Rscript -e "devtools::load_all(); ui <- mixdra:::binary_ui('binary'); stopifnot(inherits(ui, c('shiny.tag','shiny.tag.list','bslib_fragment'))); cat('UI OK\n')"`
Expected: prints `UI OK` with no error.

- [ ] **Step 6: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): 3-stage UX (seed -> compare-all hero -> inspect)"
```

---

## Task 6: Docs, memory, and full filtered verification

**Files:**
- Modify: `R/app-binary.R:1-3` (top-of-file flow comment)
- Modify: memory files under `C:\Users\jelle\.claude\projects\D--sam\memory\`

- [ ] **Step 1: Update the module header comment**

Replace `R/app-binary.R` lines 1-3 with:

```r
# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows), seed each chemical's curve, then "Fit & compare all models"
# fits CA/IA reference + SA/DR/DL JOINTLY (curves + interaction together),
# warm-started down the nesting chain, and selects the parsimonious winner.
```

- [ ] **Step 2: Update the staged-fitting memory**

Edit `C:\Users\jelle\.claude\projects\D--sam\memory\mixdra-staged-fitting.md`: add a line noting the binary app's model **selection** loop is now fully joint (warm-started chain via `analyse_mixture_joint`/`joint_fit_one`); `analyse_mixture` (staged) is retained as an engine/test entry point but is no longer what the binary tab uses to choose a model. Link `[[excel-joint-fitting]]`.

- [ ] **Step 3: Update the Excel-joint-fitting memory**

Edit `C:\Users\jelle\.claude\projects\D--sam\memory\excel-joint-fitting.md`: note the standalone "Optimize all params (joint)" card was **removed**; joint fitting now happens inside the compare-all loop for every model (so model selection matches Excel's joint regime). Link `[[mixdra-staged-fitting]]`.

- [ ] **Step 4: Run the full set of affected test files**

Run: `Rscript -e "devtools::test(filter='compare')"`
Run: `Rscript -e "devtools::test(filter='analyse')"`
Run: `Rscript -e "devtools::test(filter='joint-loop')"`
Run: `Rscript -e "devtools::test(filter='app-modules')"`
Run: `Rscript -e "devtools::test(filter='recovery')"`
Expected: PASS for each (run individually — the user prefers filtered runs over the full suite).

- [ ] **Step 5: Commit**

```bash
git add R/app-binary.R
git commit -m "docs(binary): update flow comment for the joint compare-all loop"
```

---

## Self-Review Notes

- **Spec coverage:** joint per-model fits seeded from singles (Tasks 2-3); warm-started nesting chain for monotonic LR soundness (Task 2); live row-by-row table via `invalidateLater` stepper (Task 4); 3-stage UX with compare-all hero + demoted explore (Task 5); `Optimize all` retired (Tasks 4-5); LR statistic kept as recommended (no `compare.R` test-statistic change). All covered.
- **Decisions honored:** kept `lr_test` as-is (spec recommendation); used the stepped-`invalidateLater` mechanism (spec recommendation), not async.
- **Type consistency:** `joint_fit_one(df, reference, deviation, response, seed_curves, parent_fit, n_starts, time_limit)`, `joint_chain_order(n_chem)`, `compare_fits(fits, n, response, alpha) -> list(comparison, chosen)`, `analyse_mixture_joint(...) -> list(fits, comparison, chosen, reference, response)` are used identically across engine and app tasks. `last_compare` keeps its existing shape (`fits`, `comparison`, `chosen`, `reference`, `response`) so `output$results`/`output$isobole` need no change.
- **Reused, not duplicated:** `fit_model`, `lr_test`, `select_parsimonious`, `model_spec$parent`/`$extra`, `result_table`, `param_ci` all reused unchanged.
```
