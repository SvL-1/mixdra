# Binary Interaction-Fit Workspace Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the Binary tab's opaque "Freeze → fit interactions" step into a per-model workspace: Freeze only locks the curves, then the user picks a deviation model and **Autofits** its `a`/`b`, **Simulates** it with entered `a`/`b`, or runs **Find best model** for the four-way comparison.

**Architecture:** One small engine helper `eval_mixture()` (mixture twin of `eval_single()`, delegating to `fit_model()`'s all-fixed evaluation path) powers Simulate. The Binary tab's Stage 2 becomes a workspace mirroring the curve panel: a model picker, an `a`/`b` grid (`interaction_param_row()` — curve-panel layout with inert "—" bounds), and Autofit/Simulate/Find-best actions. The server keeps a per-model fit store keyed by model name so switching the picker (including the programmatic switch Find-best performs) never re-triggers a fit — avoiding an update feedback loop.

**Tech Stack:** R, Shiny modules (`moduleServer`/`NS`), bslib, plotly, DT, testthat (`testServer` + unit tests). Run tests **filtered** (`testthat::test_file(...)`), never the full suite. testServer tests are `skip_on_cran`, so run them with `NOT_CRAN=true`.

**Reference spec:** `docs/superpowers/specs/2026-06-02-binary-interaction-workspace-design.md`

---

## Important context for the implementer

- **Branch:** `binary-explain-fit` (already has `interaction_help()`, the relocated
  controls, and the inline explanation). This plan continues on the same branch.
- **This supersedes part of the prior work:** the auto-fit-on-freeze
  (`res_r <- eventReactive(frozen(), …)`) is removed, and the Stage-1 intro
  sentence + Fit-options group move into Stage 2. The `interaction_help()` helper
  and unit tests (`tests/testthat/test-interaction-help.R`) are unchanged and stay green.
- **Engine facts you will rely on (already verified):**
  - `model_spec(reference, deviation, n_chem)` returns `$params` (all params),
    `$extra` (interaction params), `$fn`. For binary: `reference` extra `none`;
    `SA` → `a`; `DR`/`DL` → `a`, `b`.
  - `fit_model(df, reference, deviation, response, start, fixed, n_starts, time_limit)`
    fits one model; when every param is in `fixed` (so `free` is empty) it
    **evaluates only** and returns `list(par, objective, pred, residuals, df = 0,
    n, convergence = 0, reference, deviation, response, conc_cols, n_chem,
    kind = "mixture")`. Bounds may name **only** curve params — never `a`/`b`.
  - `analyse_mixture(df, reference, response, start = curve_params, alpha,
    n_starts, time_limit)` returns `list(fits, comparison, chosen, reference,
    response)`; `comparison` has columns `model/parent/chi/df/p` (3 rows for
    binary).
  - Plots consume a single fit: `plot_surface(fit, df)`, `plot_obs_pred(fit, df)`,
    `plot_isobole(fit, df, reference_fit = …)` (accepts `reference_fit = NULL`).
    `result_table(res)` takes a full `analyse_mixture` result; `param_ci(fit, df,
    reference, deviation, response)` returns `parameter/estimate/lower/upper`.
  - `assemble_curve_params(fit1, fit2)`, `marginal_df()`, `curve_fit_server()` are
    unchanged and still produce `curve_params()` (`max`, `slope1`, `slope2`,
    `ec501`, `ec502`).
- **Commit scope:** stage only the files named in each task. Do **not** stage the
  unrelated working-tree files (`single_cpf_continuous.csv`, the
  `R/app-run.R` scroll fix, or any edits to `analyse.R`/`registry.R`/`test-nchem.R`).
- **No `Co-Authored-By` trailer. Do not push.**

---

## File Structure

- **`R/fit.R`** (modify) — add `eval_mixture()` next to `fit_model()`.
- **`R/app-curve-fit.R`** (modify) — add `interaction_param_row()` next to `param_row()` (sibling grid-row helper).
- **`R/app-binary.R`** (modify) — Freeze-only button; Stage 2 workspace; server rewire (fit store, Autofit/Simulate/Find-best, rewired Stage 3).
- **`tests/testthat/test-eval-mixture.R`** (create) — unit tests for `eval_mixture()`.
- **`tests/testthat/test-app-modules.R`** (modify) — rewrite the `binary_server` testServer flow + the `binary_ui` smoke test; add an `interaction_param_row()` smoke test.

---

## Task 1: `eval_mixture()` engine helper

The Simulate primitive. Assemble the full parameter vector from the frozen curve params + entered interaction params and evaluate the model without optimising, returning the enriched `fit_model()` shape (plus `simulated = TRUE`).

**Files:**
- Create: `tests/testthat/test-eval-mixture.R`
- Modify: `R/fit.R` (add `eval_mixture()` after `fit_model()`)

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-eval-mixture.R`:

```r
# eval_mixture(): forward-evaluate a mixture model at fixed curve + interaction
# params, returning the fit_model() shape. Uses the noiseless CA surfaces from
# the engine's own vectorised model functions as ground truth.

cp <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)

make_grid <- function() {
  d <- expand.grid(C1 = c(0, 0.04, 0.08, 0.16), C2 = c(0, 0.5, 1, 2))
  d
}

test_that("eval_mixture reproduces the CA reference surface (objective ~ 0)", {
  d <- make_grid()
  d$Res <- ca_bi_vec(d$C1, d$C2, cp[["max"]], cp[["slope1"]], cp[["slope2"]],
                     cp[["ec501"]], cp[["ec502"]])
  fit <- eval_mixture(d, "CA", "reference", "continuous", cp)
  expect_equal(fit$objective, 0, tolerance = 1e-8)
  expect_equal(unname(fit$pred), d$Res, tolerance = 1e-8)
  expect_true(fit$simulated)
  expect_equal(fit$kind, "mixture")
  expect_equal(fit$deviation, "reference")
  expect_setequal(names(fit$par), c("max", "slope1", "slope2", "ec501", "ec502"))
})

test_that("eval_mixture reproduces a known CA + S/A surface from entered a", {
  d <- make_grid()
  d$Res <- ca_sa_bi_vec(d$C1, d$C2, cp[["max"]], cp[["slope1"]], cp[["slope2"]],
                        cp[["ec501"]], cp[["ec502"]], a = 5)
  fit <- eval_mixture(d, "CA", "SA", "continuous", cp, interaction = c(a = 5))
  expect_equal(fit$objective, 0, tolerance = 1e-6)
  expect_equal(fit$par[["a"]], 5)
  expect_true("a" %in% names(fit$par))
  expect_equal(length(fit$pred), nrow(d))
})

test_that("eval_mixture drops interaction params the model does not use", {
  d <- make_grid()
  d$Res <- ca_sa_bi_vec(d$C1, d$C2, cp[["max"]], cp[["slope1"]], cp[["slope2"]],
                        cp[["ec501"]], cp[["ec502"]], a = 5)
  # pass a stray `b`; S/A has no b, so it must be ignored (not error)
  fit <- eval_mixture(d, "CA", "SA", "continuous", cp, interaction = c(a = 5, b = 99))
  expect_false("b" %in% names(fit$par))
  expect_equal(fit$par[["a"]], 5)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-eval-mixture.R')"`
Expected: FAIL — `could not find function "eval_mixture"`.

- [ ] **Step 3: Implement `eval_mixture()`**

In `R/fit.R`, after the end of `fit_model()` (after its closing `}`), add:

```r
#' Forward-evaluate a mixture model at fixed parameters (Simulate)
#'
#' Mixture counterpart of [eval_single()]. Assembles the full parameter vector
#' from the frozen curve parameters plus caller-supplied interaction parameters
#' and evaluates the model without optimising, returning the enriched shape of
#' [fit_model()] (so the plotting layer consumes it unchanged) with an added
#' `simulated = TRUE` flag. Delegates to [fit_model()]'s all-fixed path.
#'
#' @param df Mixture data frame (`C1`/`C2`[/`C3`] + response columns).
#' @param reference "CA" or "IA".
#' @param deviation "reference", "SA", "DR", or "DL".
#' @param response "continuous" or "binary".
#' @param curve_params Named numeric of curve params (max, slope1, slope2,
#'   ec501, ec502).
#' @param interaction Named numeric of interaction params (e.g. `c(a = .., b = ..)`).
#'   Entries the chosen model does not use are dropped.
#' @return An enriched mixture fit list (as [fit_model()]) with `simulated = TRUE`.
#' @keywords internal
eval_mixture <- function(df, reference, deviation, response,
                         curve_params, interaction = numeric(0)) {
  n_chem <- length(intersect(c("C1", "C2", "C3"), names(df)))
  spec <- model_spec(reference, deviation, n_chem)
  use  <- interaction[intersect(names(interaction), spec$extra)]
  par  <- c(curve_params, use)
  fit  <- fit_model(df, reference, deviation, response,
                    start = par, fixed = spec$params)   # all fixed -> evaluate only
  fit$simulated <- TRUE
  fit
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-eval-mixture.R')"`
Expected: PASS — all assertions green, no warnings.

- [ ] **Step 5: Commit**

```bash
git add R/fit.R tests/testthat/test-eval-mixture.R
git commit -m "feat(engine): eval_mixture() forward-evaluation helper for Simulate"
```

---

## Task 2: `interaction_param_row()` grid-row helper

A presentational helper for the `a`/`b` grid: the same five-column layout as the curve panel's `param_row()`, but with inert "—" placeholders in the Lower/Upper columns (the engine leaves `a`/`b` unconstrained) and only an editable Value input.

**Files:**
- Modify: `R/app-curve-fit.R` (add `interaction_param_row()` after `param_row()`)
- Modify: `tests/testthat/test-app-modules.R` (add a smoke test)

- [ ] **Step 1: Write the failing smoke test**

In `tests/testthat/test-app-modules.R`, add at the end of the file:

```r
test_that("interaction_param_row renders a value input with inert bound cells", {
  ns <- shiny::NS("binary")
  html <- as.character(interaction_param_row(ns, "a", "a", "overall strength"))
  expect_match(html, "binary-val_a", fixed = TRUE)   # editable Value input present
  expect_match(html, "—", fixed = TRUE)          # em-dash placeholder present
  expect_false(grepl("binary-lo_a", html, fixed = TRUE))  # no Lower input
  expect_false(grepl("binary-hi_a", html, fixed = TRUE))  # no Upper input
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', desc='interaction_param_row renders a value input with inert bound cells')"`
Expected: FAIL — `could not find function "interaction_param_row"`.

- [ ] **Step 3: Implement `interaction_param_row()`**

In `R/app-curve-fit.R`, immediately after the `param_row()` function (after its closing `}`), add:

```r
#' One interaction-parameter row: label, meaning, inert bound cells, value input
#'
#' Same five-column layout as [param_row()] for visual consistency with the
#' curve grid, but `a`/`b` are unconstrained by the engine, so the Lower/Upper
#' columns are inert "—" placeholders and only the Value input is editable.
#' @param ns Module namespace function.
#' @param param Parameter key (`a`/`b`); drives the `val_<param>` input id.
#' @param label Display label.
#' @param meaning One-line explanation.
#' @keywords internal
interaction_param_row <- function(ns, param, label, meaning) {
  shiny::fluidRow(
    shiny::column(2, shiny::tags$b(label)),
    shiny::column(4, shiny::tags$small(meaning)),
    shiny::column(2, shiny::tags$small("—")),
    shiny::column(2, shiny::tags$small("—")),
    shiny::column(2, shiny::numericInput(ns(paste0("val_", param)), NULL, value = NA))
  )
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', desc='interaction_param_row renders a value input with inert bound cells')"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app-curve-fit.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): interaction_param_row() a/b grid-row helper"
```

---

## Task 3: Rebuild `binary_ui` — Freeze-only button + Stage 2 workspace

Rename the Freeze button (no longer fits), remove the Stage-1 intro sentence and Fit-options group added previously, and build the Stage 2 workspace: model picker, `interaction_help` block, `a`/`b` grid, Autofit/Simulate buttons + objective readout, the relocated Fit-options group, and the Find-best group with the comparison table. Stage 3 keeps its outputs but the CI card stays a `DTOutput`.

**Files:**
- Modify: `R/app-binary.R` (`binary_ui` only — server is Task 4)
- Modify: `tests/testthat/test-app-modules.R` (replace the `binary_ui` smoke test)

- [ ] **Step 1: Replace the `binary_ui` smoke test (failing)**

In `tests/testthat/test-app-modules.R`, replace the existing test that starts
`test_that("binary_ui builds and exposes the relocated fit options + intro", {`
(through its closing `})`) with:

```r
test_that("binary_ui: Freeze-only button + Stage 2 workspace", {
  html <- as.character(binary_ui("binary"))
  # Freeze no longer fits
  expect_match(html, "Freeze curves", fixed = TRUE)
  expect_false(grepl("Freeze curves → fit interactions", html, fixed = TRUE))
  # Stage 2 workspace: model picker, a/b value input, the three actions
  expect_match(html, "binary-model", fixed = TRUE)
  expect_match(html, "binary-val_a", fixed = TRUE)
  expect_match(html, "Autofit (a, b)", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "Find best model", fixed = TRUE)
  # relocated fit options + the per-model explanation slot
  expect_match(html, "binary-n_starts", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', desc='binary_ui: Freeze-only button + Stage 2 workspace')"`
Expected: FAIL — `binary-model`, `Autofit (a, b)`, `Find best model` not present; the old "Freeze curves → fit interactions" string still present.

- [ ] **Step 3: Restructure Stage 1 (Freeze-only button)**

In `R/app-binary.R`, in the Stage 1 `bslib::card(...)`, replace the current block
(the intro `shiny::p(...)`, the `Fit options` group, the `thorough` checkbox, the
`thorough_note`, and the Freeze button) — i.e. replace this:

```r
      shiny::p(class = "text-muted",
               "Freezing locks both curves above and fits only the interaction terms ",
               "(a, b) to the mixture rows, then compares four models ",
               "(no interaction → S/A → dose-ratio → dose-level) ",
               "and flags the most parsimonious."),
      shiny::tags$b("Fit options"),
      bslib::layout_columns(
        shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
        shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
        shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1)
      ),
      shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
      shiny::uiOutput(ns("thorough_note")),
      shiny::actionButton(ns("freeze"), "Freeze curves → fit interactions",
                          class = "btn-primary"),
      shiny::uiOutput(ns("freeze_note"))
    ),
```

with this:

```r
      shiny::actionButton(ns("freeze"), "Freeze curves", class = "btn-primary"),
      shiny::uiOutput(ns("freeze_note"))
    ),
```

- [ ] **Step 4: Rebuild the Stage 2 card as the workspace**

In `R/app-binary.R`, replace the entire Stage 2 `conditionalPanel(...)` block —
i.e. replace this:

```r
    # Stage 2 -- interaction-model comparison + model picker (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Interaction models"),
        shiny::selectInput(ns("model"), "Model to display", choices = NULL),
        shiny::helpText("The most parsimonious model is pre-selected; ",
                        "override above to inspect another."),
        shiny::uiOutput(ns("interaction_help")),
        DT::DTOutput(ns("comparison"))
      )
    ),
```

with this:

```r
    # Stage 2 -- per-model interaction workspace (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Interaction model"),
        shiny::selectInput(
          ns("model"), "Model",
          choices = c("No interaction (reference)" = "reference",
                      "Similar action (S/A)"       = "SA",
                      "Dose-ratio (DR)"            = "DR",
                      "Dose-level (DL)"            = "DL")),
        shiny::uiOutput(ns("interaction_help")),

        # a/b parameter grid (hidden for the reference model; b shown for DR/DL).
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
        ),
        shiny::uiOutput(ns("objective")),

        shiny::tags$b("Fit options"),
        bslib::layout_columns(
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1)
        ),
        shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
        shiny::uiOutput(ns("thorough_note")),

        shiny::hr(),
        shiny::p(shiny::tags$b("Find best model"),
                 shiny::tags$small(" — fit all four models and pick the most parsimonious.")),
        bslib::layout_columns(
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::actionButton(ns("find_best"), "Find best model")
        ),
        DT::DTOutput(ns("comparison"))
      )
    ),
```

- [ ] **Step 5: Run the smoke test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', desc='binary_ui: Freeze-only button + Stage 2 workspace')"`
Expected: PASS.

Note: the `binary_server` testServer test in the same file will now FAIL (it
still references the removed `res_r`/`shown_fit`). That is expected — Task 4
rewrites it. Do not run the whole file green yet.

- [ ] **Step 6: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): Freeze-only button + Stage 2 interaction workspace UI"
```

---

## Task 4: Rewire `binary_server` — fit store, Autofit/Simulate/Find-best, Stage 3

Remove the auto-fit on freeze. Add a per-model fit store and a derived `current_fit`, plus `last_compare`. Wire the three action observers, the picker→grid sync, and rewire Stage 3 to read `current_fit()`/`last_compare()`. Drive `interaction_help` from `input$model`.

**Files:**
- Modify: `R/app-binary.R` (`binary_server`)
- Modify: `tests/testthat/test-app-modules.R` (rewrite the `binary_server` testServer test)

- [ ] **Step 1: Rewrite the `binary_server` testServer test (failing)**

In `tests/testthat/test-app-modules.R`, replace the whole test that starts
`test_that("binary_server stages the fit: split, freeze gate, average max, override, invalidate", {`
(through its closing `})`) with:

```r
test_that("binary_server workspace: freeze gate, autofit/simulate/find-best, invalidate", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
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
    expect_equal(unname(curve_params()[["max"]]), 650)

    # freeze locks curves but does NOT fit
    session$setInputs(freeze = 1)
    expect_true(isTRUE(frozen()))
    expect_null(current_fit())

    # Autofit the selected (S/A) model -> a fit for SA is stored
    session$setInputs(model = "SA")
    session$setInputs(autofit = 1)
    expect_equal(current_fit()$deviation, "SA")
    expect_true("a" %in% names(current_fit()$par))

    # Simulate S/A with an entered a -> a simulated fit replaces the stored one
    session$setInputs(val_a = 3, simulate = 1)
    expect_true(isTRUE(current_fit()$simulated))
    expect_equal(current_fit()$par[["a"]], 3)

    # Find best model -> comparison of all four + a chosen model
    session$setInputs(find_best = 1)
    expect_setequal(names(last_compare()$fits), c("reference", "SA", "DR", "DL"))
    expect_equal(nrow(last_compare()$comparison), 3)
    expect_true(last_compare()$chosen %in% c("reference", "SA", "DR", "DL"))

    # editing a single curve after freezing invalidates everything
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_false(isTRUE(frozen()))
    expect_equal(length(fits_store()), 0)
    expect_null(last_compare())
  })
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "Sys.setenv(NOT_CRAN='true'); devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R', desc='binary_server workspace: freeze gate, autofit/simulate/find-best, invalidate')"`
Expected: FAIL — `current_fit`/`fits_store`/`last_compare` not found (still the old `res_r` server).

- [ ] **Step 3: Replace the freeze/fit/state block in `binary_server`**

In `R/app-binary.R`, replace the block that begins at the freeze checkpoint
comment and runs through the `res_r` definition — i.e. replace this:

```r
    # Checkpoint state. `frozen` gates Stages 2-3 (exposed to the UI as an output).
    frozen <- shiny::reactiveVal(FALSE)
    output$frozen <- shiny::reactive(isTRUE(frozen()))
    shiny::outputOptions(output, "frozen", suspendWhenHidden = FALSE)

    # Frozen curve-parameter vector (shared max = average of the two fits).
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2())
      assemble_curve_params(fit1(), fit2())
    })

    output$freeze_note <- shiny::renderUI({
      if (is.null(fit1()) || is.null(fit2()))
        shiny::div(class = "text-muted",
                   shiny::tags$small(
                     "Fit both single curves (Autofit or Simulate) before freezing."))
    })

    # Freeze requires both curves; marking frozen triggers the interaction fit.
    shiny::observeEvent(input$freeze, {
      if (is.null(fit1()) || is.null(fit2())) {
        shiny::showNotification("Fit both single curves before freezing.", type = "warning")
        return()
      }
      frozen(TRUE)
    })

    # A changed curve or header model invalidates the freeze (result can't go stale).
    shiny::observeEvent(
      list(fit1(), fit2(), input$reference, input$response),
      { if (isTRUE(frozen())) frozen(FALSE) },
      ignoreInit = TRUE)

    # Stage 2/3: fit reference + deviations with the frozen curve params fixed.
    res_r <- shiny::eventReactive(frozen(), {
      shiny::req(isTRUE(frozen()), fit1(), fit2())
      df <- engine_df()
      engine_response <- if (input$response == "quantal") "binary" else "continuous"
      n_starts <- if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts
      start <- curve_params()
      shiny::withProgress(message = "Fitting interaction models...", value = 0.5, {
        tryCatch(
          analyse_mixture(df, reference = input$reference, response = engine_response,
                          start = start, alpha = input$alpha,
                          n_starts = n_starts, time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
            NULL
          })
      })
    })

    # Refresh the model picker after each fit, defaulting to the chosen model.
    shiny::observeEvent(res_r(), {
      shiny::updateSelectInput(session, "model",
                               choices = names(res_r()$fits), selected = res_r()$chosen)
    })

    # Fall back to the chosen model until the picker's input has populated.
    shown_fit <- shiny::reactive({
      shiny::req(res_r())
      m <- input$model
      if (is.null(m) || !m %in% names(res_r()$fits)) m <- res_r()$chosen
      res_r()$fits[[m]]
    })

    # Per-model explanation: tracks the header reference and the displayed model.
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen(), res_r())
      interaction_help(input$reference, shown_fit()$deviation)
    })
```

with this:

```r
    # Checkpoint state. `frozen` gates Stages 2-3 (exposed to the UI as an output).
    frozen <- shiny::reactiveVal(FALSE)
    output$frozen <- shiny::reactive(isTRUE(frozen()))
    shiny::outputOptions(output, "frozen", suspendWhenHidden = FALSE)

    # Frozen curve-parameter vector (shared max = average of the two fits).
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2())
      assemble_curve_params(fit1(), fit2())
    })

    # Per-model interaction fits, keyed by model name. Autofit/Simulate write one
    # entry; Find best writes all four. `current_fit` is whatever is stored for the
    # selected model (NULL if none yet). Keeping a store -- rather than a single
    # reactiveVal that Find best would stomp when it switches the picker -- means
    # the programmatic picker switch never re-triggers or clears a fit.
    fits_store   <- shiny::reactiveVal(list())
    last_compare <- shiny::reactiveVal(NULL)
    current_fit  <- shiny::reactive({
      m <- input$model
      if (is.null(m)) return(NULL)
      fits_store()[[m]]
    })

    engine_response <- shiny::reactive(
      if (input$response == "quantal") "binary" else "continuous")
    n_starts_eff <- shiny::reactive(
      if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts)

    output$freeze_note <- shiny::renderUI({
      if (is.null(fit1()) || is.null(fit2()))
        shiny::div(class = "text-muted",
                   shiny::tags$small(
                     "Fit both single curves (Autofit or Simulate) before freezing."))
    })

    # Freeze only locks the curves and reveals Stage 2 -- it does not fit.
    shiny::observeEvent(input$freeze, {
      if (is.null(fit1()) || is.null(fit2())) {
        shiny::showNotification("Fit both single curves before freezing.", type = "warning")
        return()
      }
      frozen(TRUE)
    })

    # A changed curve or header model invalidates the freeze and clears all fits.
    shiny::observeEvent(
      list(fit1(), fit2(), input$reference, input$response),
      {
        if (isTRUE(frozen())) {
          frozen(FALSE)
          fits_store(list())
          last_compare(NULL)
        }
      },
      ignoreInit = TRUE)

    # Read the entered a/b as a named numeric (blank -> NA).
    read_ab <- function() {
      raw <- list(a = input$val_a, b = input$val_b)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

    # Autofit: fit only the selected model's interaction params, curves fixed.
    shiny::observeEvent(input$autofit, {
      shiny::req(frozen(), curve_params())
      dev <- input$model
      fit <- tryCatch(
        shiny::withProgress(message = "Fitting interaction...", value = 0.5,
          fit_model(engine_df(), input$reference, dev, engine_response(),
                    start = curve_params(), fixed = names(curve_params()),
                    n_starts = n_starts_eff(), time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      if ("a" %in% names(fit$par))
        shiny::updateNumericInput(session, "val_a", value = round(fit$par[["a"]], 4))
      if ("b" %in% names(fit$par))
        shiny::updateNumericInput(session, "val_b", value = round(fit$par[["b"]], 4))
      s <- fits_store(); s[[dev]] <- fit; fits_store(s)
    })

    # Simulate: evaluate the selected model with the entered a/b (no refit).
    shiny::observeEvent(input$simulate, {
      shiny::req(frozen(), curve_params())
      dev <- input$model
      if (dev == "reference") {
        shiny::showNotification(
          "The reference model has no interaction parameters to simulate.", type = "message")
        return()
      }
      ab <- read_ab()
      need <- if (dev == "SA") "a" else c("a", "b")
      if (any(is.na(ab[need]))) {
        shiny::showNotification("Enter a (and b) to simulate.", type = "warning")
        return()
      }
      fit <- tryCatch(
        eval_mixture(engine_df(), input$reference, dev, engine_response(),
                     curve_params(), interaction = ab[need]),
        error = function(e) {
          shiny::showNotification(paste("Simulate failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      s <- fits_store(); s[[dev]] <- fit; fits_store(s)
    })

    # Find best model: fit all four + select; store every fit, land on the chosen.
    shiny::observeEvent(input$find_best, {
      shiny::req(frozen(), curve_params())
      res <- tryCatch(
        shiny::withProgress(message = "Comparing interaction models...", value = 0.5,
          analyse_mixture(engine_df(), reference = input$reference,
                          response = engine_response(), start = curve_params(),
                          alpha = input$alpha, n_starts = n_starts_eff(),
                          time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(res)) return()
      last_compare(res)
      s <- fits_store()
      for (m in names(res$fits)) s[[m]] <- res$fits[[m]]
      fits_store(s)
      shiny::updateSelectInput(session, "model", selected = res$chosen)
      ch <- res$fits[[res$chosen]]
      shiny::updateNumericInput(session, "val_a",
        value = if ("a" %in% names(ch$par)) round(ch$par[["a"]], 4) else NA)
      shiny::updateNumericInput(session, "val_b",
        value = if ("b" %in% names(ch$par)) round(ch$par[["b"]], 4) else NA)
    })

    # Switching the picker syncs the grid to that model's stored a/b (blank if none).
    shiny::observeEvent(input$model, {
      f <- fits_store()[[input$model]]
      shiny::updateNumericInput(session, "val_a",
        value = if (!is.null(f) && "a" %in% names(f$par)) round(f$par[["a"]], 4) else NA)
      shiny::updateNumericInput(session, "val_b",
        value = if (!is.null(f) && "b" %in% names(f$par)) round(f$par[["b"]], 4) else NA)
    }, ignoreInit = TRUE)

    # Per-model explanation: tracks the reference and the selected model.
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen())
      interaction_help(input$reference, input$model)
    })

    # Fit-objective readout for the displayed model.
    output$objective <- shiny::renderUI({
      shiny::req(current_fit())
      f <- current_fit()
      lab <- if (identical(f$response, "binary")) "Deviance" else "SSR"
      shiny::tags$p(
        shiny::tags$b(paste0(lab, ": ")), round(f$objective, 2),
        "   |   ", shiny::tags$b("n: "), f$n,
        if (isTRUE(f$simulated)) shiny::tags$em(" (simulated)"))
    })
```

- [ ] **Step 4: Rewire Stage 3 outputs to `current_fit()` / `last_compare()`**

In `R/app-binary.R`, replace the Stage 2/3 output block — i.e. replace this:

```r
    output$surface <- plotly::renderPlotly({
      shiny::req(frozen(), res_r()); plot_surface(shown_fit(), engine_df())
    })
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen(), res_r())
      plot_isobole(shown_fit(), engine_df(), reference_fit = res_r()$fits$reference)
    })
    output$op <- plotly::renderPlotly({
      shiny::req(frozen(), res_r()); plot_obs_pred(shown_fit(), engine_df())
    })

    output$results <- DT::renderDT({
      shiny::req(frozen(), res_r())
      tab <- round(result_table(res_r()), 4)
      DT::datatable(as.data.frame(tab), options = list(dom = "t"))
    })
    output$comparison <- DT::renderDT({
      shiny::req(frozen(), res_r())
      DT::datatable(res_r()$comparison, rownames = FALSE, options = list(dom = "t"))
    })
    output$cis <- DT::renderDT({
      shiny::req(frozen(), shown_fit())
      f <- shown_fit()
      DT::datatable(param_ci(f, engine_df(), f$reference, f$deviation, f$response),
                    rownames = FALSE, options = list(dom = "t"))
    })

    res_r   # return for testability
```

with this:

```r
    output$surface <- plotly::renderPlotly({
      shiny::req(frozen(), current_fit()); plot_surface(current_fit(), engine_df())
    })
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen(), current_fit())
      ref <- if (!is.null(last_compare())) last_compare()$fits$reference else NULL
      plot_isobole(current_fit(), engine_df(), reference_fit = ref)
    })
    output$op <- plotly::renderPlotly({
      shiny::req(frozen(), current_fit()); plot_obs_pred(current_fit(), engine_df())
    })

    # Table-2 style all-models matrix -- only meaningful after Find best.
    output$results <- DT::renderDT({
      shiny::req(last_compare())
      tab <- round(result_table(last_compare()), 4)
      DT::datatable(as.data.frame(tab), options = list(dom = "t"))
    })
    output$comparison <- DT::renderDT({
      shiny::req(last_compare())
      DT::datatable(last_compare()$comparison, rownames = FALSE, options = list(dom = "t"))
    })
    # CIs for a fitted model; a simulated (hand-entered) set has no CIs, so show
    # its entered values instead (honest -- those are the numbers you set).
    output$cis <- DT::renderDT({
      shiny::req(current_fit())
      f <- current_fit()
      tab <- if (isTRUE(f$simulated))
        data.frame(parameter = names(f$par), value = round(unname(f$par), 4))
      else
        param_ci(f, engine_df(), f$reference, f$deviation, f$response)
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"))
    })

    list(current_fit = current_fit, last_compare = last_compare)  # return for testability
```

- [ ] **Step 5: Run the testServer test to verify it passes**

Run: `Rscript -e "Sys.setenv(NOT_CRAN='true'); devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R')"`
Expected: PASS — the workspace test, the `binary_ui` smoke test, the
`interaction_param_row` smoke test, and the unchanged single/intro tests all
green (0 fail). Pre-existing `param_ci`/crosstalk warnings may appear; they are
not failures.

- [ ] **Step 6: Run the full affected set + a load check**

Run: `Rscript -e "Sys.setenv(NOT_CRAN='true'); devtools::load_all('.'); testthat::test_file('tests/testthat/test-eval-mixture.R'); testthat::test_file('tests/testthat/test-interaction-help.R'); testthat::test_file('tests/testthat/test-app-modules.R'); cat('load OK\n')"`
Expected: all three files PASS (0 fail); `load OK` printed.

- [ ] **Step 7: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): per-model interaction workspace server (autofit/simulate/find-best)"
```

---

## Notes for the implementer

- **Run tests filtered, never the full suite** (per project convention). testServer
  tests are `skip_on_cran`; set `NOT_CRAN=true` to run them.
- **`current_fit` is derived from a per-model store**, refining spec §6's two
  reactiveVals. This is deliberate: a single `current_fit` reactiveVal would be
  cleared by the `input$model` sync observer when Find-best programmatically
  switches the picker, wiping the result Find-best just set. The store avoids the
  loop and, as a bonus, lets the user flip between all four models after Find best
  and see each stored fit. Observable behaviour still matches the spec.
- **The CI panel for a simulated fit** shows the entered `parameter`/`value` pairs
  (no CI columns) rather than the spec's exact note sentence — this realises the
  spec's intent (don't present CIs on hand-entered values) in a DT-friendly way.
- **No `Co-Authored-By` trailer. Do not push.** Stage only the files named per task.
```
