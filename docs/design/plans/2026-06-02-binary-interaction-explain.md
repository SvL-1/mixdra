# Binary Interaction-Fit Legibility Implementation Plan

**Goal:** Make the Binary tab's Freeze → interaction-fit step legible — relocate the four fit controls to the Freeze step, explain the step inline, and show the per-model `a`/`b` formula for the displayed deviation model.

**Architecture:** Add one pure helper `interaction_help(reference, deviation)` in `R/app-binary.R` that returns Shiny markup (no reactivity, no fit object) — this is the unit-testable core. Then move the four controls (`alpha`, `n_starts`, `thorough`, `time_limit`) from the sidebar "Advanced" accordion into a "Fit options" group inside the Stage 1 card next to the Freeze button (same input ids, so server logic is untouched), add an always-visible intro sentence at the Freeze step, and wire a `renderUI("interaction_help")` block in Stage 2 driven by the header reference and the displayed model.

**Tech Stack:** R, Shiny modules (`moduleServer`/`NS`), bslib (`layout_sidebar`, `card`), testthat (`testServer` + plain unit tests). Tests are run filtered, e.g. `devtools::test(filter = "interaction-help")` — never the full suite.

**Reference spec:** `docs/design/specs/2026-06-01-binary-interaction-explain-design.md`

---

## File Structure

- **`R/app-binary.R`** (modify) — add the pure `interaction_help()` helper near the top (alongside `binary_ui`); restructure `binary_ui` to drop the sidebar Advanced accordion and add the Fit-options group + intro sentence in Stage 1 and the `interaction_help` output in Stage 2; add `output$interaction_help` to `binary_server`. The four controls keep their exact ids (`alpha`, `n_starts`, `thorough`, `time_limit`), so `res_r()` / `analyse_mixture(...)` are unchanged.
- **`tests/testthat/test-interaction-help.R`** (create) — unit tests for the pure `interaction_help()` helper (no Shiny server).
- **`tests/testthat/test-app-modules.R`** (modify) — extend the existing `binary_server` `testServer` test to assert `output$interaction_help` follows the picker, and extend the `binary_ui` smoke test to assert the Fit-options labels + intro sentence + that the old "Advanced" accordion is gone.

**Out of scope (do not touch):** `R/app-single.R` (the dropped "About this model" panel is a separate tracked follow-up); `R/app-run.R` (the `fillable = FALSE` scroll fix is a separate working-tree change); any of the engine fitting logic. Do **not** stage or commit unrelated working-tree files (`single_cpf_continuous.csv`, or any in-progress edits to `analyse.R`/`registry.R`/`test-nchem.R`). Stage only the files named in each task.

---

## Task 1: The pure `interaction_help()` helper

Add a pure function that returns the per-model explanation markup. It takes the reference (`"CA"`/`"IA"`) and the deviation key (`"reference"`/`"SA"`/`"DR"`/`"DL"`) and returns a `shiny::tagList`. No reactivity, no fit object — fully testable on its own. Mirrors the offline-safe plain-HTML style of `single_model_equation()` / the old `model_help_single()` (HTML entities, no MathJax).

**Files:**
- Create: `tests/testthat/test-interaction-help.R`
- Modify: `R/app-binary.R` (add helper after the file header comment, before `binary_ui`)

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-interaction-help.R`:

```r
# interaction_help() is a pure helper: (reference, deviation) -> Shiny markup.
# These tests render the tagList to an HTML string and assert the formula
# fragment + a/b wording for each reference x deviation combination.

html_of <- function(reference, deviation) {
  as.character(interaction_help(reference, deviation))
}

test_that("reference deviation states there is no interaction term", {
  h <- html_of("CA", "reference")
  expect_match(h, "No interaction", ignore.case = TRUE)
  expect_match(h, "CA", fixed = TRUE)
})

test_that("reference deviation names the IA baseline when reference is IA", {
  expect_match(html_of("IA", "reference"), "IA", fixed = TRUE)
})

test_that("SA shows the F = a.product formula and the a meaning", {
  h <- html_of("CA", "SA")
  expect_match(h, "F = a&middot;&prod;z", fixed = TRUE)
  expect_match(h, "strength and direction", fixed = TRUE)
})

test_that("DR shows a summation formula and the mixture-ratio meaning", {
  h <- html_of("CA", "DR")
  expect_match(h, "&Sigma;", fixed = TRUE)
  expect_match(h, "mixture ratio", fixed = TRUE)
})

test_that("DL formula differs between CA (SigmaTU) and IA (P)", {
  ca <- html_of("CA", "DL")
  ia <- html_of("IA", "DL")
  expect_match(ca, "&Sigma;TU", fixed = TRUE)
  expect_match(ia, "b&middot;P)", fixed = TRUE)
  expect_match(ia, "IA-predicted", fixed = TRUE)
})

test_that("the sign convention is stated for models with an interaction term", {
  for (dev in c("SA", "DR", "DL")) {
    h <- html_of("CA", dev)
    expect_match(h, "antagonism", fixed = TRUE)
    expect_match(h, "synergism", fixed = TRUE)
  }
})

test_that("the sign convention is omitted for the reference (no a term)", {
  h <- html_of("CA", "reference")
  expect_false(grepl("antagonism", h, fixed = TRUE))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-interaction-help.R')"`
Expected: FAIL — `could not find function "interaction_help"`.

- [ ] **Step 3: Implement the helper**

In `R/app-binary.R`, immediately after the top-of-file comment block and before `#' Binary Mixture stage UI`, add:

```r
#' Per-model explanation of the interaction fit (pure, for the Binary tab)
#'
#' Returns the formula and `a`/`b` meaning for one reference x deviation
#' combination. Plain HTML (no MathJax), matching the offline-safe style of
#' [single_model_equation()]. No reactivity or fit object -- unit-testable.
#'
#' The interaction factor `F` adjusts the chosen reference baseline on the
#' toxic-unit shares `z = TU / SigmaTU` (`TU = C / EC50`; `prod z` is the
#' product across chemicals).
#'
#' @param reference Reference model: `"CA"` or `"IA"`.
#' @param deviation Deviation key: `"reference"`, `"SA"`, `"DR"`, or `"DL"`.
#' @return A [shiny::tagList()] of formula + parameter meanings.
#' @keywords internal
interaction_help <- function(reference, deviation) {
  ref <- if (identical(reference, "IA")) "IA" else "CA"

  intro <- shiny::tags$p(shiny::tags$small(shiny::HTML(paste0(
    "The interaction factor F adjusts the ", ref,
    " baseline on the toxic-unit shares ",
    "(z = TU / &Sigma;TU, TU = C / EC50; &prod;z = product across chemicals)."))))

  sign_note <- shiny::tags$p(shiny::tags$small(
    shiny::tags$b("Sign of a: "),
    "a > 0 → antagonism (mixture less toxic than the reference predicts); ",
    "a < 0 → synergism (more toxic)."))

  body <- switch(
    deviation,
    reference = shiny::tags$p(shiny::HTML(paste0(
      "<b>No interaction term.</b> The mixture follows the ", ref,
      " reference exactly."))),
    SA = shiny::tagList(
      shiny::tags$p(shiny::tags$b("Similar action (S/A): "),
                    shiny::tags$code(shiny::HTML("F = a&middot;&prod;z"))),
      shiny::tags$p(shiny::tags$small(
        shiny::tags$code("a"),
        " sets the overall strength and direction of the interaction."))),
    DR = shiny::tagList(
      shiny::tags$p(shiny::tags$b("Dose-ratio dependent (DR): "),
                    shiny::tags$code(shiny::HTML(
                      "F = (a + &Sigma;b<sub>i</sub>&middot;z<sub>i</sub>)&middot;&prod;z"))),
      shiny::tags$p(shiny::tags$small(
        shiny::tags$code("a"), " = overall interaction; ", shiny::tags$code("b"),
        " = how it shifts with the mixture ratio (which chemical dominates). ",
        "For a binary mixture there is one free b (Jonker Eq. 8)."))),
    DL = {
      dl_form <- if (ref == "IA")
        "F = a&middot;(1 &minus; b&middot;P)&middot;&prod;z"
      else
        "F = a&middot;(1 &minus; b&middot;&Sigma;TU)&middot;&prod;z"
      dl_b <- if (ref == "IA")
        "how it shifts with the dose level (P = the IA-predicted effect at the mixture point)."
      else
        "how it shifts with the dose level (ΣTU = the summed toxic units at the mixture point)."
      shiny::tagList(
        shiny::tags$p(shiny::tags$b("Dose-level dependent (DL): "),
                      shiny::tags$code(shiny::HTML(dl_form))),
        shiny::tags$p(shiny::tags$small(
          shiny::tags$code("a"), " = overall interaction; ",
          shiny::tags$code("b"), " = ", dl_b)))
    },
    shiny::tags$p()  # unknown deviation -> empty
  )

  shiny::tagList(intro, body,
                 if (!identical(deviation, "reference")) sign_note)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-interaction-help.R')"`
Expected: PASS — all assertions green, no warnings.

- [ ] **Step 5: Commit**

```bash
git add R/app-binary.R tests/testthat/test-interaction-help.R
git commit -m "feat(binary): interaction_help() per-model fit explanation helper"
```

---

## Task 2: Relocate the four controls + intro sentence into the Stage 1 card; wire the Stage 2 output

Restructure `binary_ui`: remove the sidebar `bslib::accordion("Advanced", …)` block and the sidebar `thorough` checkbox + `thorough_note`, and add — in the Stage 1 card, after the two curve columns and before the Freeze button — a static intro sentence and a "Fit options" group holding all four controls (same ids and defaults as before). Add a `shiny::uiOutput(ns("interaction_help"))` block in the Stage 2 card next to the model picker. The sidebar keeps response, reference, template, file, and errors.

**Files:**
- Modify: `R/app-binary.R` (`binary_ui`)
- Modify: `tests/testthat/test-app-modules.R` (the `binary_ui` smoke test)

- [ ] **Step 1: Write the failing UI smoke test**

In `tests/testthat/test-app-modules.R`, replace the existing one-line `binary_ui` smoke test:

```r
test_that("binary_ui builds a Shiny UI fragment", {
  expect_true(inherits(binary_ui("binary"), c("shiny.tag", "shiny.tag.list", "bslib_fragment")))
})
```

with:

```r
test_that("binary_ui builds and exposes the relocated fit options + intro", {
  html <- as.character(binary_ui("binary"))
  # the four controls now live in the Stage 1 Fit-options group (ids unchanged)
  expect_match(html, "binary-alpha", fixed = TRUE)
  expect_match(html, "binary-n_starts", fixed = TRUE)
  expect_match(html, "binary-thorough", fixed = TRUE)
  expect_match(html, "binary-time_limit", fixed = TRUE)
  # the group label and the always-visible intro sentence are present
  expect_match(html, "Fit options", fixed = TRUE)
  expect_match(html, "fits only the interaction terms", fixed = TRUE)
  # the old sidebar "Advanced" accordion is gone
  expect_false(grepl("Advanced", html, fixed = TRUE))
  # Stage 2 hosts the per-model explanation output
  expect_match(html, "binary-interaction_help", fixed = TRUE)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R')"`
Expected: FAIL — `Advanced` still present / `Fit options` and `binary-interaction_help` not found.

- [ ] **Step 3: Restructure the sidebar (remove Advanced + thorough)**

In `R/app-binary.R`, replace the sidebar block — change this:

```r
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
      shiny::uiOutput(ns("thorough_note")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced",
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1)
        )
      ),
      shiny::uiOutput(ns("errors"))
    ),
```

to this:

```r
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::uiOutput(ns("errors"))
    ),
```

- [ ] **Step 4: Add the intro sentence + Fit-options group to the Stage 1 card**

In `R/app-binary.R`, in the Stage 1 `bslib::card(...)`, replace the Freeze button line — change this:

```r
      shiny::actionButton(ns("freeze"), "Freeze curves → fit interactions",
                          class = "btn-primary"),
      shiny::uiOutput(ns("freeze_note"))
```

to this:

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
```

- [ ] **Step 5: Add the interaction_help output to the Stage 2 card**

In `R/app-binary.R`, in the Stage 2 `conditionalPanel(...)` card, replace this:

```r
        shiny::selectInput(ns("model"), "Model to display", choices = NULL),
        shiny::helpText("The most parsimonious model is pre-selected; ",
                        "override above to inspect another."),
        DT::DTOutput(ns("comparison"))
```

with this:

```r
        shiny::selectInput(ns("model"), "Model to display", choices = NULL),
        shiny::helpText("The most parsimonious model is pre-selected; ",
                        "override above to inspect another."),
        shiny::uiOutput(ns("interaction_help")),
        DT::DTOutput(ns("comparison"))
```

- [ ] **Step 6: Run the smoke test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R')"`
Expected: The `binary_ui` smoke test PASSES. (The `binary_server` `testServer` test in the same file still passes — control ids are unchanged. It does not yet assert `interaction_help`; Task 3 adds that.)

- [ ] **Step 7: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): relocate fit controls to the Freeze step + intro sentence"
```

---

## Task 3: Render the per-model explanation in `binary_server` and assert it follows the picker

Add `output$interaction_help` to `binary_server`, rendering `interaction_help(input$reference, shown_fit()$deviation)` once frozen. `shown_fit()` already resolves the displayed deviation (picker value, else `res_r()$chosen`) and is gated by `req(res_r())`, so the block renders only after a successful freeze and updates when the reference or the picker changes.

**Files:**
- Modify: `R/app-binary.R` (`binary_server`)
- Modify: `tests/testthat/test-app-modules.R` (the `binary_server` `testServer` test)

- [ ] **Step 1: Extend the failing testServer assertions**

In `tests/testthat/test-app-modules.R`, inside the `binary_server` `testServer` block, after the existing picker-override assertions:

```r
    # picker overrides the displayed model
    session$setInputs(model = "DR")
    expect_equal(shown_fit()$deviation, "DR")
```

add (before the invalidation assertions that follow):

```r
    # the per-model explanation follows the picker. With DR shown, it is the
    # DR block and the DL-only fragment (&Sigma;TU) is absent.
    expect_match(output$interaction_help$html, "Dose-ratio dependent", fixed = TRUE)
    expect_false(grepl("&Sigma;TU", output$interaction_help$html, fixed = TRUE))
    # switching the picker to DL switches the block; with reference CA the DL
    # form is the &Sigma;TU variant.
    session$setInputs(model = "DL")
    expect_equal(shown_fit()$deviation, "DL")
    expect_match(output$interaction_help$html, "Dose-level dependent", fixed = TRUE)
    expect_match(output$interaction_help$html, "&Sigma;TU", fixed = TRUE)
```

Note: `output$interaction_help` in `testServer` is a list whose `$html` field holds the rendered HTML string (the standard shape for a `renderUI` output under `testServer`). The reference here is `"CA"` (set in the existing `session$setInputs`), so DL renders the `&Sigma;TU` form.

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R')"`
Expected: FAIL — `output$interaction_help` is `NULL` (output not defined yet), so `$html` is `NULL` and `expect_match` errors.

- [ ] **Step 3: Add the renderUI output to `binary_server`**

In `R/app-binary.R`, immediately after the `shown_fit <- shiny::reactive({...})` block (which ends with `res_r()$fits[[m]]` then `})`), add:

```r
    # Per-model explanation: tracks the header reference and the displayed model.
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen(), res_r())
      interaction_help(input$reference, shown_fit()$deviation)
    })
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-app-modules.R')"`
Expected: PASS — both the `binary_server` `testServer` test and the `binary_ui` smoke test green, no warnings.

- [ ] **Step 5: Run the interaction-help unit tests + a load check to confirm nothing regressed**

Run: `Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-interaction-help.R'); cat('load OK\n')"`
Expected: PASS — all unit tests green; `load OK` printed.

- [ ] **Step 6: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): render per-model interaction explanation in Stage 2"
```

---

## Notes for the implementer

- **Run tests filtered, never the full suite.** Use `testthat::test_file(...)` on the specific files as shown, or `devtools::test(filter = "interaction-help")` / `devtools::test(filter = "app-modules")`. The full suite is too slow for this workflow.
- **Commit scope.** Stage only the files named in each task's commit step. The working tree contains unrelated in-progress files (`single_cpf_continuous.csv` and possibly edits to `analyse.R`/`registry.R`/`test-nchem.R`) and the uncommitted `R/app-run.R` scroll fix — do **not** add them.
- **No behaviour change to the fit.** The four control input ids (`alpha`, `n_starts`, `thorough`, `time_limit`) and their defaults are identical to before; `res_r()` reads them exactly as it did. Only their UI location moves. `output$thorough_note` already exists in the server and is unchanged — only its UI position moved (it now sits in the Fit-options group).
- **No Co-Authored-By trailer** on any commit.
- Do not push. All work stays local on branch `binary-explain-fit`.
