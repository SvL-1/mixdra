# Binary Table-Centric Stage 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Stage 2 comparison table the single hub: it auto-fills when both single curves exist (no "Compare all" button), clicking a row shows that model in the diagnostics (no model picker), and one "Optimize all params (joint)" button refines the selected row (a display-only joint column appears) — all while the p-values / winner stay computed from the staged fits only.

**Architecture:** Reuse the existing `refined_fits`/`display_fit`/`optimize_all` server core and the live stepper. Replace the three triggering observers (compare-all kick, reference/response clear, curve re-eval) with ONE auto-trigger on `curve_params()`/`reference`/`response`. Derive the displayed model from DT single-row selection (`sel_row`) instead of a `selectInput`. Add a display-only joint-SSR column to the results table. Remove the manual a/b path (Autofit/Simulate/grid).

**Tech Stack:** R package (`mixdra`), Shiny + bslib + DT + plotly, testthat (`devtools::test`). DT single-row selection via `selection=list(mode="single")` + `input$results_rows_selected`.

**Spec:** `docs/superpowers/specs/2026-06-03-binary-table-centric-stage2-design.md`

---

## File Structure

- **`R/app-binary.R`** — all changes (server + `binary_ui`).
- **`tests/testthat/test-app-modules.R`** — evolve the binary server tests (auto-fill, row-select, verdict isolation) and the `binary_ui` HTML test.

No engine changes. No NAMESPACE changes. The `interaction_param_row` helper becomes unused by `binary_ui` but is kept (it still has a unit test and is harmless).

---

## Task 1: Auto-fill the table (replace 3 trigger observers with 1)

**Files:**
- Modify: `R/app-binary.R` (`binary_server`)
- Test: `tests/testthat/test-app-modules.R`

The table should fill in on its own once both curves are fit, and re-fill when a curve / reference / response changes — no button. This replaces the compare-all kick, the reference/response clear observer, and the curve re-eval observer with one auto-trigger. The UI still has the (now-orphan) compare_all button until Task 3.

- [ ] **Step 1: Update the failing test.** In `tests/testthat/test-app-modules.R`, in the test `"binary_server: compare-all loop fills the store live and selects a model"`, make two edits:

(a) Remove the line `session$setInputs(compare_all = 1)` (keep the `for (i in 1:8) session$elapse(5)` line right after it). Update the comment above it to:
```r
    # Both curves are now fit -> the staged loop auto-runs (no button). Drain it.
```

(b) Replace the final block (currently the "editing a single curve ... expect_null(last_compare())" assertions) with:
```r
    # editing a single curve auto-re-runs the staged loop at the new curves
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    for (i in 1:8) session$elapse(5)
    expect_true(isTRUE(frozen()))
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    expect_false(is.null(last_compare()))        # comparison recomputed at the new curves
```

- [ ] **Step 2: Run test to verify it fails.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: FAIL — without `compare_all`, the loop never arms, so `fits_store()` stays empty.

- [ ] **Step 3: Replace the three observers with one auto-trigger.** In `R/app-binary.R`, delete these three consecutive blocks:
  - the `# A structural change (reference model or response type) ...` `observeEvent(list(input$reference, input$response), {...}, ignoreInit = TRUE)` block;
  - the `# A changed curve does NOT blank the interaction: ...` `observeEvent(list(fit1(), fit2()), {...}, ignoreInit = TRUE)` block (the curve re-eval observer, including its inner `for` loop);
  - the `# Kick: clear the store and arm the selection chain ...` `observeEvent(input$compare_all, {...})` block.

Replace all three with:
```r
    # Auto-fill: once both single curves are fit -- and again whenever the curves,
    # reference, or response type change -- clear the stores and re-run the staged
    # comparison loop. There is no "Compare all" button; the table fills in and
    # stays live on its own. `curve_params()` req()s both single fits, so this is
    # inert until frozen. The stepper below refits one model per tick.
    shiny::observeEvent(
      list(curve_params(), input$reference, input$response),
      {
        shiny::req(frozen())
        fits_store(list())
        last_compare(NULL)
        refined_fits(list()); refine_pre(NULL); refine_post(NULL)
        loop_queue(selection_chain_order(n_chem()))
        stepper_on(TRUE)
      },
      ignoreInit = TRUE)
```

- [ ] **Step 4: Run the test to verify it passes.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: the auto-fill test PASSES. (The `binary_ui` HTML test still passes — UI unchanged.) Report the summary line and any failure.

- [ ] **Step 5: Commit.**
```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): table auto-fills on freeze/curve change (no Compare-all button)"
```

**Notes:** the auto-trigger depends only on `curve_params()`/`reference`/`response` — none of which the stepper writes — so there is no reactive loop. Curve changes are discrete (Stage-1 Autofit/Simulate buttons), not slider drags, so re-running 4 cheap staged fits per change is fine.

---

## Task 2: Row-select drives display; joint column; remove manual path

**Files:**
- Modify: `R/app-binary.R` (`binary_server`)
- Test: `tests/testthat/test-app-modules.R`

Replace the `input$model` picker with DT row-selection, add a display-only joint-objective column, and remove the manual a/b controls (Autofit/Simulate/grid/`read_ab`/picker-sync).

- [ ] **Step 1: Rewrite the failing server test.** Replace the entire `test_that("binary_server: compare-all loop fills the store live and selects a model", {...})` block with:

```r
test_that("binary_server: table auto-fills, row-click selects, joint refine is isolated", {
  skip_on_cran()
  meta <- shiny::reactiveValues(chem1 = "A", chem2 = "B")
  shiny::testServer(binary_server, args = list(meta = meta), {
    csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                               "binary_ca_mps_cpf_imi_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA", thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      file = list(datapath = csv, name = "binary.csv"))
    session$setInputs(`chem1-val_max` = 700, `chem1-val_slope` = 2,
                      `chem1-val_ec50` = 1, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 600, `chem2-val_slope` = 1,
                      `chem2-val_ec50` = 5, `chem2-simulate` = 1)

    # auto-fill: no button, just drain the timer
    for (i in 1:8) session$elapse(5)
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    chosen     <- last_compare()$chosen
    staged_cmp <- last_compare()$comparison

    # with no row clicked, the displayed model defaults to the winner
    expect_equal(selected_model(), chosen)
    expect_equal(current_fit()$deviation, chosen)

    # clicking a row (by index in reference,SA,DR,DL order) selects that model
    ord <- intersect(c("reference", "SA", "DR", "DL"), names(fits_store()))
    session$setInputs(results_rows_selected = match("SA", ord))
    expect_equal(selected_model(), "SA")
    expect_equal(current_fit()$deviation, "SA")

    # optimise the selected (SA) model jointly -> refined_fits gets SA, verdict frozen
    staged_sa_obj <- fits_store()[["SA"]]$objective
    session$setInputs(optimize_all = 1)
    expect_true(isTRUE(refined_fits()[["SA"]]$joint))
    expect_lte(refined_fits()[["SA"]]$objective, staged_sa_obj + 1e-6)
    expect_true(isTRUE(display_fit()$joint))            # diagnostics follow the refined fit
    expect_identical(last_compare()$chosen, chosen)     # verdict unchanged
    expect_equal(last_compare()$comparison, staged_cmp)
    expect_equal(fits_store()[["SA"]]$objective, staged_sa_obj)

    # a curve edit re-runs the staged loop and clears refined fits
    session$setInputs(`chem1-val_max` = 720, `chem1-simulate` = 2)
    for (i in 1:8) session$elapse(5)
    expect_false(is.null(last_compare()))
    expect_equal(length(refined_fits()), 0)
  })
})
```

- [ ] **Step 2: Run test to verify it fails.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: FAIL — `selected_model` not found / `input$results_rows_selected` not wired.

- [ ] **Step 3: Add `results_order`, `sel_row`, `selected_model`.** In `R/app-binary.R`, replace the existing `current_fit <- shiny::reactive({...})` and `display_fit <- shiny::reactive({...})` blocks with:

```r
    # Row order in the comparison table (fixed) and the row the user clicked.
    # `sel_row` persists the clicked row across table re-renders (e.g. when the
    # joint column fills in). The displayed model defaults to the winner.
    results_order <- shiny::reactive(
      intersect(c("reference", "SA", "DR", "DL"), names(fits_store())))
    sel_row <- shiny::reactiveVal(integer(0))
    shiny::observeEvent(input$results_rows_selected, {
      sel_row(input$results_rows_selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    selected_model <- shiny::reactive({
      ord <- results_order()
      if (length(ord) == 0) return(NULL)
      r <- sel_row()
      if (length(r) >= 1 && r[[1]] >= 1 && r[[1]] <= length(ord)) return(ord[[r[[1]]]])
      ch <- if (!is.null(last_compare())) last_compare()$chosen else NULL
      if (!is.null(ch) && ch %in% ord) ch else ord[[1]]
    })

    current_fit <- shiny::reactive({
      m <- selected_model()
      if (is.null(m)) return(NULL)
      fits_store()[[m]]
    })

    # The fit shown in the diagnostics / CIs for the selected model: the
    # joint-refined fit when one exists, else the staged fit. Never read by the
    # comparison table (which stays purely staged).
    display_fit <- shiny::reactive({
      m <- selected_model()
      if (is.null(m)) return(NULL)
      r <- refined_fits()[[m]]
      if (!is.null(r)) r else fits_store()[[m]]
    })
```

- [ ] **Step 4: Default-select the winner in the stepper completion.** In the stepper's `if (length(rest) == 0) {...}` block, replace the `shiny::updateSelectInput(session, "model", ...)` line AND the two `shiny::updateNumericInput(session, "val_a"/"val_b", ...)` blocks with:
```r
          ord <- intersect(c("reference", "SA", "DR", "DL"), names(fits_store()))
          sel_row(match(cmp$chosen, ord))   # default-select the winner row
```
So the completion block reads:
```r
        if (length(rest) == 0) {
          stepper_on(FALSE)
          cmp <- compare_fits(fits_store(), nrow(engine_df()),
                              engine_response(), input$alpha)
          last_compare(list(fits = fits_store(), comparison = cmp$comparison,
                            chosen = cmp$chosen, reference = input$reference,
                            response = engine_response()))
          ord <- intersect(c("reference", "SA", "DR", "DL"), names(fits_store()))
          sel_row(match(cmp$chosen, ord))   # default-select the winner row
        }
```

- [ ] **Step 5: Point `optimize_all` at the selected model.** In the `observeEvent(input$optimize_all, {...})` block, change `r[[input$model]] <- newfit` to `r[[selected_model()]] <- newfit`. (Its `req(frozen())` + `current_fit()` guard already use the selected model via `current_fit`.)

- [ ] **Step 6: Remove the manual a/b path.** Delete entirely:
  - the `read_ab <- function() {...}` helper;
  - the `observeEvent(input$autofit, {...})` block;
  - the `observeEvent(input$simulate, {...})` block;
  - the `observeEvent(input$model, {...})` picker-sync block.

- [ ] **Step 7: Point `interaction_help` + `refined_badge` at the selected model.** Change `output$interaction_help`'s body to:
```r
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen(), selected_model())
      interaction_help(input$reference, selected_model())
    })
```
Change `output$refined_badge`'s first lines to use the selected model:
```r
    output$refined_badge <- shiny::renderUI({
      m <- selected_model()
      if (!is.null(m) && !is.null(refined_fits()[[m]]))
        shiny::div(class = "text-info", shiny::tags$small(
          "Showing the joint-refined (all-parameters) fit for this model. ",
          "The comparison table above is unchanged (staged)."))
    })
```

- [ ] **Step 8: Add the joint column + single-row selection to the results table.** Replace the `output$results <- DT::renderDT({...})` block with:
```r
    # One row per fitted model: staged params + objective/df + LR p-value, winner
    # highlighted. A display-only joint-objective column is filled from
    # refined_fits -- it is NEVER fed into the verdict (p / winner come from the
    # staged `cmp`/`chosen`). Single-row selection drives the Stage 3 diagnostics;
    # `selected = sel_row()` re-applies the selection across re-renders so
    # optimising a row (which fills its joint cell) does not lose the selection.
    output$results <- DT::renderDT({
      fits <- fits_store()
      shiny::req(length(fits) > 0)
      ord    <- intersect(c("reference", "SA", "DR", "DL"), names(fits))
      cmp    <- if (!is.null(last_compare())) last_compare()$comparison else NULL
      chosen <- if (!is.null(last_compare())) last_compare()$chosen else NULL
      mat  <- round(result_table(list(fits = fits[ord], comparison = cmp)), 4)
      disp <- as.data.frame(t(mat), check.names = FALSE)        # models -> rows
      rf  <- refined_fits()
      jlab <- if (identical(engine_response(), "binary")) "Dev (joint)" else "SSR (joint)"
      disp[[jlab]] <- vapply(rownames(disp), function(m)
        if (!is.null(rf[[m]])) round(rf[[m]]$objective, 4) else NA_real_, numeric(1))
      disp <- cbind(`Interaction model` = rownames(disp), disp, stringsAsFactors = FALSE)
      rownames(disp) <- NULL
      dt <- DT::datatable(
        disp, rownames = FALSE,
        selection = list(mode = "single", target = "row", selected = sel_row()),
        options = list(dom = "t", ordering = FALSE, scrollX = TRUE))
      if (!is.null(chosen) && chosen %in% disp[["Interaction model"]])
        dt <- DT::formatStyle(dt, "Interaction model", target = "row",
                              fontWeight = DT::styleEqual(chosen, "bold"),
                              backgroundColor = DT::styleEqual(chosen, "#d8f0d8"))
      dt
    })
```

- [ ] **Step 9: Run the test to verify it passes.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: the rewritten server test PASSES. The `binary_ui` HTML test will now FAIL (it still references the removed picker/autofit and the compare_all button) — that is EXPECTED; Task 3 fixes it. Confirm the ONLY failure is the `binary_ui` test; report it.

- [ ] **Step 10: Commit.**
```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): row-click selects model; display-only joint column; drop manual a/b path"
```

---

## Task 3: UI — table-centric Stage 2 + diagnostics Stage 3

**Files:**
- Modify: `R/app-binary.R` (`binary_ui`)
- Test: `tests/testthat/test-app-modules.R` (the `binary_ui` HTML test)

- [ ] **Step 1: Rewrite the failing UI test.** Replace the `test_that("binary_ui: comparison table + tune/refine panel with joint optimize, then diagnostics", {...})` block with:
```r
test_that("binary_ui: auto-fill table (selectable) + Optimize button, then diagnostics", {
  html <- as.character(binary_ui("binary"))
  expect_false(grepl("Freeze curves", html, fixed = TRUE))
  # Stage 2: the live table + alpha + the single Optimize button + readout/badge/help
  expect_match(html, "binary-results", fixed = TRUE)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "significance threshold for the model comparison", fixed = TRUE)
  expect_match(html, "binary-optimize_all", fixed = TRUE)
  expect_match(html, "Optimize all params (joint)", fixed = TRUE)
  expect_match(html, "binary-refine_readout", fixed = TRUE)
  expect_match(html, "binary-refined_badge", fixed = TRUE)
  expect_match(html, "binary-interaction_help", fixed = TRUE)
  expect_match(html, "binary-objective", fixed = TRUE)
  # Stage 3: diagnostics
  expect_match(html, "binary-surface", fixed = TRUE)
  expect_match(html, "binary-isobole", fixed = TRUE)
  expect_match(html, "binary-op", fixed = TRUE)
  expect_match(html, "binary-cis", fixed = TRUE)
  expect_match(html, "binary-n_starts", fixed = TRUE)
  # the manual path + the compare-all button + model picker are gone
  expect_false(grepl("binary-compare_all", html, fixed = TRUE))
  expect_false(grepl("binary-model", html, fixed = TRUE))
  expect_false(grepl("binary-autofit", html, fixed = TRUE))
  expect_false(grepl("binary-simulate", html, fixed = TRUE))
  expect_false(grepl("binary-val_a", html, fixed = TRUE))
  expect_false(grepl("binary-find_best", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails.** Run: `Rscript -e "devtools::test(filter='app-modules')"`. Expected: FAIL — `binary-compare_all`/`binary-model` still present, `binary-objective` not where expected, etc.

- [ ] **Step 3: Replace the Stage 2 + Stage 3 UI blocks.** In `binary_ui`, replace everything from the comment `# Stage 2 -- staged verdict (compare-all + table) PLUS a tune/refine panel.` through the closing `)` of the Stage 3 `conditionalPanel` with:
```r
    # Stage 2 -- the live comparison table IS the hub: it auto-fills, you click a
    # row to inspect that model below, and one button refines the selected row.
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Interaction models"),
        shiny::p("Fills in automatically once both single curves are fit: every ",
                 "interaction model is fit with the curves held fixed from the ",
                 "single compounds (reference → S/A → DR/DL), and the best is ",
                 "highlighted. Click a row to inspect that model below."),
        shiny::numericInput(ns("alpha"), "alpha", value = 0.05,
                            min = 0, max = 1, step = 0.01),
        shiny::helpText(
          "alpha (α) is the significance threshold for the model comparison: ",
          "a more complex model is kept only if it improves the fit at p < α ",
          "(default 0.05). The highlighted row is the selected (best) model. ",
          "The \"… (joint)\" column is filled by Optimize all params and is ",
          "display-only — it never changes which model is selected."),
        DT::DTOutput(ns("results")),
        shiny::uiOutput(ns("interaction_help")),
        shiny::uiOutput(ns("refined_badge")),
        shiny::uiOutput(ns("objective")),
        shiny::div(class = "mt-2",
          shiny::actionButton(ns("optimize_all"), "Optimize all params (joint)",
                              class = "btn-primary")),
        shiny::helpText("Re-fits every parameter (curves + interaction) of the ",
                        "selected row at once, seeded from its staged fit — the ",
                        "Excel-style joint fit. Updates the diagnostics below and ",
                        "the row's \"(joint)\" column; the rest of the table is ",
                        "unchanged."),
        shiny::uiOutput(ns("refine_readout"))
      )
    ),

    # Stage 3 -- diagnostics for the selected model (joint-refined if present).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 3 · Inspect the selected model"),
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

- [ ] **Step 4: Run the UI test + boot-check.** Run: `Rscript -e "devtools::test(filter='app-modules')"` (whole file should be green). Then `Rscript -e "devtools::load_all(); ui <- mixdra:::binary_ui('binary'); stopifnot(inherits(ui, c('shiny.tag','shiny.tag.list','bslib_fragment'))); cat('UI OK\n')"` → expect `UI OK`.

- [ ] **Step 5: Commit.**
```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): table-centric Stage 2 UI (auto-fill, row-select, one Optimize)"
```

---

## Task 4: Docs, memory, verification

**Files:**
- Modify: `R/app-binary.R` (module header comment)
- Modify: memory files under `C:\Users\jelle\.claude\projects\D--sam\memory\`

- [ ] **Step 1: Update the module header comment.** Replace the top-of-file comment block in `R/app-binary.R` with:
```r
# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows) and fit each chemical's curve. The Stage-2 comparison table then
# fills in automatically -- CA/IA reference + SA/DR/DL fit via the staged method
# (curves fixed from the single compounds, only a/b per model), one per tick,
# compared by LR test with the parsimonious winner highlighted. Click a row to
# inspect that model (Stage 3 diagnostics); "Optimize all params (joint)" refines
# the selected row (refine_joint, Excel-style) as a post-selection polish stored
# separately, so the staged verdict (p-values + winner) never changes.
```

- [ ] **Step 2: Update the `mixdra-staged-fitting` memory.** Edit `C:\Users\jelle\.claude\projects\D--sam\memory\mixdra-staged-fitting.md`: note the Stage-2 UX is now table-centric (auto-fills on freeze/curve-change — no Compare-all button; DT row-click selects the displayed model — no picker; one Optimize-all button refines the selected row; manual Autofit/Simulate a/b path removed). Verdict still staged.

- [ ] **Step 3: Update the `excel-joint-fitting` memory.** Edit `C:\Users\jelle\.claude\projects\D--sam\memory\excel-joint-fitting.md`: note the joint refine is now triggered per selected row from the auto-filled table; a display-only "(joint)" SSR/Dev column shows it; still stored in `refined_fits`, still never feeds the staged verdict.

- [ ] **Step 4: Run the affected suites (filtered, never the full suite).**
  - `Rscript -e "devtools::test(filter='app-modules')"`
  - `Rscript -e "devtools::test(filter='summary')"`
  Expected: PASS for both.

- [ ] **Step 5: Commit.**
```bash
git add R/app-binary.R
git commit -m "docs(binary): module comment for the table-centric Stage 2"
```

---

## Self-Review Notes

- **Spec coverage:** auto-fill on freeze/curve-change, no button (Task 1); row-click selection replacing the picker via `sel_row`/`selected_model` (Task 2 Steps 3,4,8); display-only joint column not feeding the verdict (Task 2 Step 8 — `p`/`chosen` come from `cmp`/`last_compare`, the joint column from `refined_fits`); single Optimize button on the selected row (Task 2 Step 5, Task 3); removal of compare_all button / picker / manual a/b (Task 2 Step 6, Task 3); Stage 3 diagnostics on `display_fit` (unchanged); docs/memory (Task 4). All spec sections covered.
- **Verdict isolation preserved:** `output$results` reads `refined_fits` ONLY for the display-only joint column; `compare_fits`/`last_compare` (the `p`/winner) are fed exclusively from `fits_store`. The server test asserts the comparison + chosen + staged objective are unchanged after a joint refine.
- **Selection stability:** `selected = sel_row()` re-applies the row selection on re-render (so filling the joint cell doesn't drop the selection); `sel_row` is updated by user clicks (`input$results_rows_selected`) and by the stepper completion (winner).
- **Type/name consistency:** `results_order`, `sel_row`, `selected_model`, `input$results_rows_selected`, the `jlab` joint column are used consistently across server + tests; `refine_joint`/`compare_fits`/`display_fit` signatures unchanged.
- **No reactive loop:** the auto-trigger depends on `curve_params()`/`reference`/`response` (not `fits_store`), which the stepper never writes.
```
