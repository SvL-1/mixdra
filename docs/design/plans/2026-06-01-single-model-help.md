# "About this model" Help Panel (Single Chemical tab) Implementation Plan

**Goal:** Add a collapsed "About this model" panel to the Single Chemical tab that shows the fitted three-parameter log-logistic equation and a one-line gloss for each reported value.

**Architecture:** A new internal helper `model_help_single()` returns static `shiny` markup (equation + parameter glossary + diagnostics note). `single_ui()` embeds it in a collapsed `bslib::accordion` panel just above the Parameters card. Pure markup — no server change, no fit dependency, no network dependency (plain HTML, not MathJax).

**Tech Stack:** R package; `shiny` + `bslib` (Suggests); `testthat` (edition 3); `roxygen2`.

**Design spec:** `docs/design/specs/2026-06-01-single-model-help-design.md`

---

## Running R and tests here (read first)

- **In RStudio:** open `mixdra.Rproj`, `Ctrl+Shift+L` (load_all), then run the task's test command in the Console.
- **From the CLI:** every non-interactive R call must prepend the user-library path or packages are not found:
  ```r
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  ```
  e.g. `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="app-modules")'`
- **Rtools is not installed on this machine** — `R CMD check` / package install cannot run here. `load_all`/`test`/`document` work. The helper is asserted by building the tag and inspecting its rendered HTML (no installed package or browser needed).
- The new test lives in `tests/testthat/test-app-modules.R`, which already begins with `skip_if_not_installed("shiny")` / `skip_if_not_installed("bslib")` — keep those guards covering the new test.

## File structure

| File | Responsibility |
|---|---|
| `R/app-single.R` (modify) | Add internal `model_help_single()` helper; embed it in `single_ui()` inside a collapsed accordion. |
| `tests/testthat/test-app-modules.R` (modify, append) | One test asserting the helper builds a tag and its HTML contains the model name and parameters. |
| `man/model_help_single.Rd` (generated) | roxygen output for the new internal helper. |

**Anchors in the current `R/app-single.R`:**
- `axis_label()` helper: lines 5–16 (the convention to match — roxygen block + `@keywords internal`).
- `single_ui()`: lines 21–41; the Parameters card is the last child of `layout_sidebar` (line 39).

---

## Task 1: "About this model" help panel

**Files:**
- Modify: `R/app-single.R` (add helper after `axis_label`, lines 5–16; edit `single_ui` body, line 39)
- Test: `tests/testthat/test-app-modules.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-app-modules.R`:

```r
test_that("model_help_single explains the log-logistic model and its parameters", {
  tag <- model_help_single()
  expect_true(inherits(tag, c("shiny.tag", "shiny.tag.list")))
  html <- as.character(tag)
  expect_match(html, "log-logistic")
  expect_match(html, "EC50")
  expect_match(html, "slope")
  expect_match(html, "SSR")
})

test_that("single_ui embeds the model help panel", {
  expect_match(as.character(single_ui("single")), "About this model")
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `devtools::test(filter = "app-modules")`
Expected: FAIL — `could not find function "model_help_single"`.

- [ ] **Step 3: Write the helper**

In `R/app-single.R`, insert this helper immediately after the `axis_label()` function (i.e. after the closing `}` on the current line 16, before the `#' Single Chemical stage UI` block):

```r
#' Static "About this model" help markup for the Single Chemical tab
#'
#' Explains the fitted three-parameter log-logistic curve and what each reported
#' value means. Pure markup, independent of any fit, so it is unit-testable on
#' its own. Rendered as plain HTML (no MathJax) so it works offline.
#' @return A `shiny` tag list.
#' @keywords internal
model_help_single <- function() {
  shiny::tagList(
    shiny::tags$p(
      "The Single Chemical tab fits a three-parameter ",
      shiny::tags$b("log-logistic"), " dose-response curve:"
    ),
    shiny::tags$p(shiny::tags$code(
      "Y = max / (1 + (C / EC50)", shiny::tags$sup("slope"), ")"
    )),
    shiny::tags$p(shiny::tags$small(
      "The response falls from ", shiny::tags$code("max"), " at zero dose."
    )),
    shiny::tags$dl(
      shiny::tags$dt("max"),
      shiny::tags$dd("Control / baseline response at C = 0 (the upper plateau)."),
      shiny::tags$dt("slope"),
      shiny::tags$dd("Steepness of the decline (> 0 means the response decreases with dose)."),
      shiny::tags$dt("EC50"),
      shiny::tags$dd("Concentration that halves the response.")
    ),
    shiny::tags$p(shiny::tags$small(
      shiny::tags$b("SSR"), " = residual sum of squares (goodness of fit); ",
      shiny::tags$b("n"), " = number of data points. ",
      "These describe the fit, not the curve."
    ))
  )
}
```

- [ ] **Step 4: Embed the panel in `single_ui()`**

In `R/app-single.R`, change the last child of `single_ui()`'s `layout_sidebar` from:

```r
    bslib::card(bslib::card_header("Parameters"), DT::DTOutput(ns("params")))
  )
}
```

to:

```r
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel("About this model", model_help_single())
    ),
    bslib::card(bslib::card_header("Parameters"), DT::DTOutput(ns("params")))
  )
}
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `devtools::test(filter = "app-modules")`
Expected: PASS — both new tests green, and the existing `app-modules` tests (intro/single/binary servers + UI fragments) still pass.

- [ ] **Step 6: Regenerate docs and run the full fast suite**

Run: `devtools::document()` (creates `man/model_help_single.Rd`; `NAMESPACE` is unchanged because the helper is internal/unexported)
Run: `devtools::test()`
Expected: the full fast suite passes, 0 failures. (The slow `validation-binary` file is unaffected.)

- [ ] **Step 7: Commit**

```bash
git add R/app-single.R tests/testthat/test-app-modules.R man/model_help_single.Rd
git commit -m "feat: 'About this model' help panel on the Single Chemical tab"
```

---

## Manual verification (optional, after the plan)

In RStudio after `devtools::load_all()`:

```r
run_app()
```

On the Single Chemical tab, confirm an "About this model" accordion sits above the
Parameters table, is collapsed by default, and expands to show the equation
`Y = max / (1 + (C / EC50)^slope)`, the three parameter definitions, and the
SSR/n diagnostics note.

---

## Self-review notes (for the implementer)

- **`shiny::tags`, not `htmltools::tags`:** matches the rest of `R/app-single.R`
  (which uses `shiny::tags$p`, `shiny::div`). `as.character()` on the returned
  tag list yields the HTML string the tests inspect.
- **Collapsed by default** (`open = FALSE`) — mirrors the `Advanced` accordion in
  `binary_ui()`, so an expert user is not forced to scroll past it.
- **No server change:** the helper is static, so `single_server()` is untouched
  and no existing reactive (`errs`, `fit_r`, `fit_df`) names change.
- **Binary tab is out of scope** per the spec — do not touch `binary_ui()`.
```
