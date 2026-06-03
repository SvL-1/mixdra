# Binary Per-Model Joint Refine (Post-Selection Polish) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let the user joint-refine (or hand-edit) the parameters of a selected interaction model — Excel-style, seeded from its staged fit — and see it in the diagnostics, **without ever changing the staged "which model wins" verdict** in the Stage 2 comparison table.

**Architecture:** Keep the staged fits in `fits_store` (drives the comparison table + frozen verdict). Add a SECOND store `refined_fits` (joint-refined fits keyed by model) that the table NEVER reads; the diagnostics/CIs render `display_fit()` = refined-when-present-else-staged. A new "Optimize all params (joint)" action calls the existing `refine_joint()` engine fn and writes only `refined_fits`. The manual a/b controls + the new joint button move into a "Tune / refine selected model" panel under the Stage 2 table; Stage 3 becomes pure diagnostics.

**Tech Stack:** R package (`mixdra`), Shiny + bslib + DT + plotly, testthat (`devtools::test`). Engine fn `refine_joint()` already exists in `R/fit.R`.

**Spec:** `docs/superpowers/specs/2026-06-03-binary-joint-refine-polish-design.md`

---

## File Structure

- **`R/app-binary.R`** — all changes. Server: add `refined_fits`/`refine_pre`/`refine_post` reactiveVals, `display_fit()` reactive, an `optimize_all` observer, a `refine_readout` + `refined_badge` output; point the diagnostic renderers (`surface`/`isobole`/`op`/`cis`/`objective`) at `display_fit()`; extend the three invalidation paths to clear `refined_fits`; clear a model's refined entry on manual Autofit/Simulate; reset the readout on model switch. UI: relocate the picker + param grid + Autofit/Simulate into a new "Tune / refine selected model" card under the Stage 2 table, add the "Optimize all params (joint)" button + readout + badge; Stage 3 keeps only the diagnostics.
- **`tests/testthat/test-app-modules.R`** — add a server test for verdict-isolation + refined-wins; update the `binary_ui` HTML test for the relocated tune panel + new button.

No engine changes. No NAMESPACE changes.

---

## Task 1: Server — `refined_fits` store, joint refine, display switch

**Files:**
- Modify: `R/app-binary.R` (`binary_server`)
- Test: `tests/testthat/test-app-modules.R`

This task is server-only; the `optimize_all` button doesn't exist in the UI yet (Task 2 adds it), but `testServer` can fire `input$optimize_all` regardless.

- [ ] **Step 1: Write the failing test.** Append to `tests/testthat/test-app-modules.R`:

```r
test_that("binary_server: joint refine polishes a model but never moves the staged verdict", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30, model = "reference",
                      file = list(datapath = csv, name = "binary.csv"))
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    # staged compare-all, then snapshot the verdict
    session$setInputs(compare_all = 1)
    for (i in 1:8) session$elapse(5)
    chosen        <- last_compare()$chosen
    staged_cmp    <- last_compare()$comparison
    staged_obj    <- fits_store()[[chosen]]$objective
    staged_df     <- fits_store()[[chosen]]$df

    # joint-refine the chosen model
    session$setInputs(model = chosen)
    session$setInputs(optimize_all = 1)

    rf <- refined_fits()[[chosen]]
    expect_true(isTRUE(rf$joint))                       # it is a joint fit
    expect_lte(rf$objective, staged_obj + 1e-6)         # seeded from staged -> no worse
    expect_gt(rf$df, staged_df)                          # joint frees more params

    # VERDICT FROZEN: comparison frame, chosen, and the staged fit are unchanged
    expect_identical(last_compare()$chosen, chosen)
    expect_equal(last_compare()$comparison, staged_cmp)
    expect_equal(fits_store()[[chosen]]$objective, staged_obj)
    expect_equal(fits_store()[[chosen]]$df, staged_df)

    # diagnostics now follow the refined fit for the selected model
    expect_true(isTRUE(display_fit()$joint))

    # a manual Autofit on the same model drops its refined entry (staged supersedes)
    session$setInputs(autofit = 1)
    expect_null(refined_fits()[[chosen]])
    expect_false(isTRUE(display_fit()$joint))

    # a curve edit clears all refined fits
    session$setInputs(model = chosen, optimize_all = 2)
    expect_false(is.null(refined_fits()[[chosen]]))
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    expect_equal(length(refined_fits()), 0)
  })
})
```

- [ ] **Step 2: Run test to verify it fails.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: FAIL — `could not find ... refined_fits` / `optimize_all` does nothing.

- [ ] **Step 3: Add the new reactiveVals.** In `R/app-binary.R`, immediately after the `loop_queue`/`stepper_on` declarations (the lines `loop_queue <- shiny::reactiveVal(NULL)` / `stepper_on <- shiny::reactiveVal(FALSE)`), add:

```r
    # Joint-refined fits, keyed by model name -- a POST-SELECTION polish. The
    # comparison table NEVER reads this; only the selected model's diagnostics do.
    # This is what keeps the staged verdict (table p-values + winner) untouched.
    refined_fits <- shiny::reactiveVal(list())
    refine_pre   <- shiny::reactiveVal(NULL)   # objective before the last joint refine
    refine_post  <- shiny::reactiveVal(NULL)   # objective after
```

- [ ] **Step 4: Add `display_fit()`.** Immediately after the existing `current_fit <- shiny::reactive({...})` block, add:

```r
    # The fit shown in the diagnostics / CIs for the selected model: the
    # joint-refined fit when one exists, else the staged fit. Never read by the
    # comparison table (which stays purely staged).
    display_fit <- shiny::reactive({
      m <- input$model
      if (is.null(m)) return(NULL)
      r <- refined_fits()[[m]]
      if (!is.null(r)) r else fits_store()[[m]]
    })
```

- [ ] **Step 5: Extend the two invalidation observers to clear refined fits.** In the `observeEvent(list(input$reference, input$response), {...})` block, add inside the braces (after `last_compare(NULL)`):

```r
        refined_fits(list()); refine_pre(NULL); refine_post(NULL)
```

In the `observeEvent(list(fit1(), fit2()), {...})` block (the curve re-eval observer), add as the LAST statement inside the braces (after `last_compare(NULL)`):

```r
        refined_fits(list()); refine_pre(NULL); refine_post(NULL)
```

- [ ] **Step 6: Clear refined fits when "Compare all" is kicked.** In the `observeEvent(input$compare_all, {...})` block, add after `last_compare(NULL)`:

```r
      refined_fits(list()); refine_pre(NULL); refine_post(NULL)
```

- [ ] **Step 7: Add the joint-refine observer.** Insert immediately after the `read_ab <- function() {...}` helper (before the `observeEvent(input$autofit, ...)` block):

```r
    # Optimize all params (joint): re-fit EVERY parameter of the selected model
    # at once, seeded from its staged fit (Excel-style). Stored in refined_fits
    # ONLY -- never in fits_store -- so the staged comparison verdict (table
    # p-values + highlighted winner) is never affected. This is a post-selection
    # parameter polish, not model selection.
    shiny::observeEvent(input$optimize_all, {
      shiny::req(frozen(), current_fit())
      f   <- current_fit()
      pre <- f$objective
      newfit <- tryCatch(
        shiny::withProgress(message = "Optimizing all parameters...", value = 0.5,
          refine_joint(f, engine_df(),
                       n_starts = n_starts_eff(), time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Optimize failed:", conditionMessage(e)),
                                  type = "error")
          NULL
        })
      if (is.null(newfit)) return()
      r <- refined_fits(); r[[input$model]] <- newfit; refined_fits(r)
      refine_pre(pre); refine_post(newfit$objective)
    })
```

- [ ] **Step 8: Make manual Autofit/Simulate drop the model's refined entry.** In the `observeEvent(input$autofit, {...})` block, after the line `s <- fits_store(); s[[dev]] <- fit; fits_store(s)`, add:

```r
      r <- refined_fits(); r[[dev]] <- NULL; refined_fits(r)  # staged edit supersedes a prior joint refine
      refine_pre(NULL); refine_post(NULL)
```

In the `observeEvent(input$simulate, {...})` block, after its `s <- fits_store(); s[[dev]] <- fit; fits_store(s)` line, add the same two lines.

- [ ] **Step 9: Reset the refine readout on model switch.** In the `observeEvent(input$model, {...})` block, add as the last statement inside the braces:

```r
      refine_pre(NULL); refine_post(NULL)
```

- [ ] **Step 10: Point the diagnostic renderers at `display_fit()`.** Make these exact replacements in `R/app-binary.R`:

In `output$objective`:
```r
    output$objective <- shiny::renderUI({
      shiny::req(display_fit())
      f <- display_fit()
      lab <- if (identical(f$response, "binary")) "Deviance" else "SSR"
      shiny::tags$p(
        shiny::tags$b(paste0(lab, ": ")), round(f$objective, 2),
        "   |   ", shiny::tags$b("n: "), f$n,
        if (isTRUE(f$simulated)) shiny::tags$em(" (simulated)")
        else if (isTRUE(f$joint)) shiny::tags$em(" (joint-refined)"))
    })
```

`output$surface`:
```r
    output$surface <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit()); plot_surface(display_fit(), engine_df())
    })
```

`output$isobole`:
```r
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit())
      ref <- if (!is.null(last_compare())) last_compare()$fits$reference else NULL
      plot_isobole(display_fit(), engine_df(), reference_fit = ref)
    })
```

`output$op`:
```r
    output$op <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit()); plot_obs_pred(display_fit(), engine_df())
    })
```

`output$cis` — change the single line `shiny::req(current_fit())` to `shiny::req(display_fit())` and `f <- current_fit()` to `f <- display_fit()` (leave the rest of the block, including the `blank_pinned_ci` branch, intact — a `refine_joint` result has `joint = TRUE` and `fixed = character(0)`, so `blank_pinned_ci` is a no-op).

- [ ] **Step 11: Add the readout + badge outputs.** Insert after the `output$objective` block:

```r
    # Before -> after objective for the last joint refine of the selected model.
    output$refine_readout <- shiny::renderUI({
      shiny::req(!is.null(refine_post()), display_fit())
      lab <- if (identical(display_fit()$response, "binary")) "Deviance" else "SSR"
      improved <- refine_post() <= refine_pre() + 1e-9
      shiny::tags$p(
        shiny::tags$b(paste0("Joint refine ", lab, ": ")),
        round(refine_pre(), 2), shiny::HTML(" &rarr; "), round(refine_post(), 2),
        if (improved) shiny::tags$span(style = "color:green", " ✓"))
    })

    # Badge: the diagnostics are showing the joint-refined fit for this model.
    output$refined_badge <- shiny::renderUI({
      m <- input$model
      if (!is.null(m) && !is.null(refined_fits()[[m]]))
        shiny::div(class = "text-info", shiny::tags$small(
          "Showing the joint-refined (all-parameters) fit for this model. ",
          "The comparison table above is unchanged (staged)."))
    })
```

- [ ] **Step 12: Expose the new reactives for testability.** Change the module's final return line from:
```r
    list(current_fit = current_fit, last_compare = last_compare)  # return for testability
```
to:
```r
    list(current_fit = current_fit, last_compare = last_compare,
         refined_fits = refined_fits, display_fit = display_fit)  # return for testability
```

- [ ] **Step 13: Run the test to verify it passes.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: the new `joint refine` test PASSES and the existing `compare-all loop` test still PASSES. The `binary_ui:` HTML test may still pass (UI unchanged this task) — confirm no NEW failures; report any.

- [ ] **Step 14: Commit.**
```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): per-model joint refine into refined_fits (verdict stays staged)"
```

---

## Task 2: UI — "Tune / refine selected model" panel under the table

**Files:**
- Modify: `R/app-binary.R` (`binary_ui`)
- Test: `tests/testthat/test-app-modules.R` (the `binary_ui` HTML test)

- [ ] **Step 1: Update the failing UI test.** In `tests/testthat/test-app-modules.R`, REPLACE the whole `test_that("binary_ui: 3-stage layout with compare-all hero and explore accordion", {...})` block with:

```r
test_that("binary_ui: comparison table + tune/refine panel with joint optimize, then diagnostics", {
  html <- as.character(binary_ui("binary"))
  expect_false(grepl("Freeze curves", html, fixed = TRUE))
  # Stage 2: staged verdict hero (compare-all + alpha + table)
  expect_match(html, "binary-compare_all", fixed = TRUE)
  expect_match(html, "Fit &amp; compare all models", fixed = TRUE)
  expect_match(html, "binary-results", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "significance threshold for the model comparison", fixed = TRUE)
  # Tune / refine panel (under the table): picker, help, a/b grid, three actions
  expect_match(html, "Tune / refine selected model", fixed = TRUE)
  expect_match(html, "binary-model", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
  expect_match(html, "binary-val_a", fixed = TRUE)
  expect_match(html, "Autofit (a, b)", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "binary-optimize_all", fixed = TRUE)          # NEW joint button
  expect_match(html, "Optimize all params (joint)", fixed = TRUE)
  expect_match(html, "binary-refine_readout", fixed = TRUE)
  expect_match(html, "binary-refined_badge", fixed = TRUE)
  # Stage 3: diagnostics
  expect_match(html, "binary-surface", fixed = TRUE)
  expect_match(html, "binary-isobole", fixed = TRUE)
  expect_match(html, "binary-op", fixed = TRUE)
  expect_match(html, "binary-cis", fixed = TRUE)
  expect_match(html, "binary-n_starts", fixed = TRUE)
  # the old staged Find-best button is still gone
  expect_false(grepl("binary-find_best", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: FAIL — "Tune / refine selected model" / `binary-optimize_all` not found.

- [ ] **Step 3: Replace the Stage 2 + Stage 3 UI.** In `R/app-binary.R` `binary_ui`, replace the two `conditionalPanel(...)` blocks — from the comment `# Stage 2 -- the hero: ...` through the closing `)` of the Stage 3 `conditionalPanel` (i.e. the entire current Stage 2 and Stage 3) — with:

```r
    # Stage 2 -- staged verdict (compare-all + table) PLUS a tune/refine panel.
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Compare interaction models"),
        shiny::p("Fits every interaction model with the curves held fixed from ",
                 "the single compounds (only the interaction terms are fitted), ",
                 "reference → S/A → DR/DL, and compares them. Rows fill in as each ",
                 "model finishes; the best model is highlighted."),
        bslib::layout_columns(
          col_widths = c(7, 5),
          shiny::div(class = "mt-4",
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
      ),

      # Tune / refine the SELECTED model. This is a post-selection polish: it
      # never changes the staged verdict above. Pick a model, then optimise its
      # parameters jointly (Excel-style, all params at once) or edit a/b by hand.
      bslib::card(
        bslib::card_header("Tune / refine selected model"),
        shiny::helpText("Refine the parameters of one model. ",
                        "\"Optimize all params (joint)\" re-fits every parameter ",
                        "(curves + interaction) at once, seeded from the staged ",
                        "fit — the Excel-style fit. The comparison table above ",
                        "stays staged and unchanged."),
        shiny::selectInput(
          ns("model"), "Model",
          choices = c("No interaction (reference)" = "reference",
                      "Similar action (S/A)"       = "SA",
                      "Dose-ratio (DR)"            = "DR",
                      "Dose-level (DL)"            = "DL")),
        shiny::uiOutput(ns("interaction_help")),
        shiny::uiOutput(ns("refined_badge")),
        shiny::uiOutput(ns("objective")),
        shiny::conditionalPanel(
          condition = "input.model != 'reference'", ns = ns,
          bslib::card(
            bslib::card_header("Interaction parameters"),
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
          shiny::actionButton(ns("simulate"), "Simulate"),
          shiny::actionButton(ns("optimize_all"), "Optimize all params (joint)")
        ),
        shiny::uiOutput(ns("refine_readout"))
      )
    ),

    # Stage 3 -- diagnostics for the selected model (joint-refined if present).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 3 · Inspect a model"),
        bslib::layout_columns(
          bslib::card(bslib::card_header("3-D response surface"),
                      plotly::plotlyOutput(ns("surface"), height = "520px")),
          bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                      plotly::plotlyOutput(ns("isobole"), height = "520px"))
        ),
        bslib::card(bslib::card_header("Observed vs predicted"),
                    plotly::plotlyOutput(ns("op"))),
        bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                    DT::DTOutput(ns("cis")))
      )
    )
```

(Structural note: the Stage 2 `conditionalPanel` now contains TWO `bslib::card`s — the comparison card and the tune card — separated by a comma. Stage 3 is the final element before the closing `)` of `layout_sidebar`; no trailing comma after it.)

- [ ] **Step 4: Run the UI test to verify it passes.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: the new `binary_ui` test PASSES and both server tests still PASS.

- [ ] **Step 5: Boot-check the UI.** Run: `Rscript -e "devtools::load_all(); ui <- mixdra:::binary_ui('binary'); stopifnot(inherits(ui, c('shiny.tag','shiny.tag.list','bslib_fragment'))); cat('UI OK\n')"`. Expected: `UI OK`.

- [ ] **Step 6: Commit.**
```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): tune/refine panel under the table (joint optimize + manual)"
```

---

## Task 3: Docs, memory, verification

**Files:**
- Modify: `R/app-binary.R` (module header comment)
- Modify: memory files under `C:\Users\jelle\.claude\projects\D--sam\memory\`

- [ ] **Step 1: Update the module header comment.** Replace the top-of-file comment block in `R/app-binary.R` with:

```r
# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows), fit each chemical's curve, then "Fit & compare all models" fits
# CA/IA reference + SA/DR/DL via the staged method (curves fixed from the single
# compounds, only a/b fitted per model -- one model per tick for the live table),
# compares them by LR test, and selects the parsimonious winner. A selected model
# can then be JOINT-refined (refine_joint, all params at once, Excel-style) as a
# post-selection polish -- stored separately so the staged verdict never changes.
```

- [ ] **Step 2: Update the `excel-joint-fitting` memory.** Edit `C:\Users\jelle\.claude\projects\D--sam\memory\excel-joint-fitting.md`: append a note that on branch `binary-joint-selection-loop` the joint fit was **re-added as a per-model post-selection polish** ("Optimize all params (joint)" in a "Tune / refine selected model" panel under the staged comparison table), stored in a separate `refined_fits` so the staged verdict is never affected; the standalone Stage-3 "Explore by hand" accordion was folded into this panel.

- [ ] **Step 3: Update the `mixdra-staged-fitting` memory.** Edit `C:\Users\jelle\.claude\projects\D--sam\memory\mixdra-staged-fitting.md`: in the 2026-06-03 rebuild note, add that a per-model **joint refine** was layered back on as a *post-selection* polish (verdict stays staged; diagnostics show the refined fit via `display_fit`), so staged-selection + joint-reporting now coexist deliberately.

- [ ] **Step 4: Run the affected suites.** Run each (filtered, never the full suite — the maintainer finds it too slow):
  - `Rscript -e "devtools::test(filter='app-modules')"`
  - `Rscript -e "devtools::test(filter='summary')"`  (param_ci on the joint fit)
  Expected: PASS for both.

- [ ] **Step 5: Commit.**
```bash
git add R/app-binary.R
git commit -m "docs(binary): module comment for the joint-refine post-selection polish"
```

---

## Self-Review Notes

- **Spec coverage:** two stores with the table never reading `refined_fits` (Task 1 Steps 3,12 + test Step 1) → verdict isolation; `display_fit` refined-wins (Step 4,10); `optimize_all` via `refine_joint` seeded from staged (Step 7); before→after readout + badge (Step 11); invalidation on reference/response, curve edit, compare-all (Steps 5,6) and manual-supersedes (Step 8); tune panel under the table + Stage 3 diagnostics (Task 2); module/memory docs (Task 3). All spec sections covered.
- **Non-goals honored:** the comparison table renderer (`output$results`) is NOT modified — it still reads only `fits_store`/`last_compare`; no curve write-back to Stage 1 (no `inject` reintroduced); no new engine fn.
- **Type/name consistency:** `refined_fits`/`refine_pre`/`refine_post`/`display_fit`/`output$refine_readout`/`output$refined_badge`/`input$optimize_all` are used identically across server, UI, and tests. `refine_joint(fit, df, n_starts, time_limit)` matches the existing `R/fit.R` signature; it returns a fit with `joint = TRUE`, `fixed = character(0)` (so the existing `blank_pinned_ci` branch in `output$cis` is a safe no-op).
- **Reuse:** `refine_joint`, `param_ci`, `interaction_param_row`, `plot_*` reused unchanged.
```
