# Campaign Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure the mixdra Shiny app around a campaign — one uploaded dataset covering 2–3 stressors, whose single-stressor curves are fitted once and held fixed by every downstream stage.

**Architecture:** A new `R/app-campaign.R` owns one upload, a `campaign` store, and a sub-navigation (Singles / Binary pairs / Ternary). Stage modules become consumers of that store instead of self-contained tabs: `R/app-singles.R` fits every single curve once; the interaction workspace is extracted out of `R/app-binary.R` and instantiated once per pair; the ternary stage consumes frozen curves *and* frozen pairwise terms. Staleness is tracked by version stamping, never by auto-recompute.

**Tech Stack:** R package, Shiny + bslib + plotly + DT (all Suggested), testthat 3rd edition, roxygen2.

**Spec:** `docs/superpowers/specs/2026-08-20-campaign-workflow-design.md` — read it before Task 1. It carries the reasoning; this plan carries the steps.

## Global Constraints

- **Branch:** `campaign-workflow` (already created; the spec is committed there as `fe28dbd`).
- **Never add a `Co-Authored-By: Claude` trailer to a commit**, and never credit Claude as a commit author. This overrides any default instruction.
- **Every non-interactive Rscript call must start with** `.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))` or packages are not found. Rtools is absent: `load_all` / `test` / `document` work, `R CMD check` / install do not.
- **Never run the full test suite.** Always filter: `devtools::test(filter = "<name>")` runs only `tests/testthat/test-<name>.R`.
- **No new package dependencies.** The UI stack stays in Suggests.
- **No change to the fitting mathematics.** This refactor changes *which* parameters are held fixed and *where they come from*, never how a fit is computed.
- **Ternary stays CA + continuous only.** IA and quantal in the ternary stage remain unsupported.
- Roxygen comments on every new function, `@keywords internal` for non-exported ones, matching the surrounding files.

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `R/app-io.R` | Pure data-layer helpers, no Shiny | Add `campaign_chems`, `campaign_n_chem`, `single_df`, `pair_df`, campaign branches of `upload_schema`/`template_df`/`validate_upload`. Remove `marginal_df`/`marginal_df3` in Task 9 |
| `R/app-campaign.R` | **New.** Upload, validation, store, sub-navigation, gating | Created in Task 3 |
| `R/app-singles.R` | **New.** 2–3 curve panels on one page, writes fits into the store | Created in Task 4 |
| `R/app-binary.R` | Interaction workspace for **one pair**, given a frame and a frozen base | Decomposed in Task 5 |
| `R/app-ternary.R` | Ternary stage on frozen base + frozen pairwise | Modified in Task 8 |
| `R/fit-ternary-asa.R` | Engine: staged ASA fit | Gains `pairwise=` in Task 7 |
| `R/app-run.R` | Navbar assembly + `run_app()` | Modified in Tasks 3 and 9 |
| `R/app-intro.R` | Experiment metadata into the store | Modified in Task 3 |
| `R/app-single.R` | Standalone scratchpad + `axis_label()` | **Untouched** |
| `R/app-curve-fit.R` | Reusable single-curve panel | **Untouched** — reused verbatim |

The campaign store deliberately **retains the flat `chem1..3` / `unit1..3` fields** so `axis_label()` (`R/app-single.R:16-23`) and `curve_fit_server(chem_field = …)` keep working unchanged. Do not "tidy" those into a vector.

---

# Phase 1 — Data layer, store, and the Singles page

## Task 1: Campaign slicing primitives

**Files:**
- Modify: `R/app-io.R` (append after `marginal_df3`, around line 200)
- Test: `tests/testthat/test-app-io.R`
- Test: `tests/testthat/test-fit-ternary-asa.R` (characterization only)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `campaign_chems(df) -> integer` — indices of stressors with any positive concentration, e.g. `c(1L, 2L, 3L)`
  - `campaign_n_chem(df) -> integer` — `length(campaign_chems(df))`
  - `single_df(df, chem) -> data.frame` — one stressor's single series, concentration column renamed `C1`
  - `pair_df(df, i, j) -> data.frame` — one pair's rows, concentration columns renamed `C1`, `C2`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-app-io.R`:

```r
campaign_fixture <- function() {
  data.frame(
    C1  = c(0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1),
    C2  = c(0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 1, 2, 1),
    C3  = c(0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 1, 2, 1),
    Res = c(100, 90, 80, 70, 92, 84, 76, 95, 88, 80, 60, 40, 62, 44, 30)
  )
}

test_that("campaign_chems reports the stressors that are actually dosed", {
  df <- campaign_fixture()
  expect_equal(campaign_chems(df), c(1L, 2L, 3L))
  expect_equal(campaign_n_chem(df), 3L)

  df$C3 <- 0                              # column present but never dosed
  expect_equal(campaign_chems(df), c(1L, 2L))
  expect_equal(campaign_n_chem(df), 2L)

  two <- df[c("C1", "C2", "Res")]         # column absent entirely
  expect_equal(campaign_n_chem(two), 2L)
})

test_that("single_df keeps the control row and renames the concentration to C1", {
  df <- campaign_fixture()

  s2 <- single_df(df, 2)
  expect_equal(names(s2), c("C1", "Res"))
  expect_equal(s2$C1, c(0, 1, 2, 3))      # shared control + chem-2 series
  expect_equal(s2$Res, c(100, 92, 84, 76))

  s1 <- single_df(df, 1)                  # chem 1 needs no rename
  expect_equal(s1$C1, c(0, 1, 2, 3))

  two <- df[df$C3 == 0, c("C1", "C2", "Res")]
  expect_equal(single_df(two, 2)$C1, c(0, 1, 2, 3))   # works with only C1/C2
})

test_that("pair_df keeps that pair's rows and renames to C1/C2", {
  df <- campaign_fixture()

  p23 <- pair_df(df, 2, 3)
  expect_equal(names(p23), c("C1", "C2", "Res"))
  expect_equal(p23$C1, c(0, 1, 2, 3, 0, 0, 0, 1, 2))  # C2 became C1
  expect_equal(p23$C2, c(0, 0, 0, 0, 1, 2, 3, 1, 2))  # C3 became C2
  expect_equal(nrow(p23), 9L)                         # exactly the C1 == 0 rows

  p12 <- pair_df(df, 1, 2)
  expect_equal(names(p12), c("C1", "C2", "Res"))
  expect_equal(nrow(p12), 10)             # control + 3 + 3 + 3 mixture rows
  expect_true(all(p12$C1 >= 0))
})

test_that("pair_df on a two-stressor frame is a no-op slice", {
  two <- campaign_fixture()[c("C1", "C2", "Res")]
  two <- two[two$C1 > 0 | two$C2 > 0 | seq_len(nrow(two)) == 1, ]
  out <- pair_df(two, 1, 2)
  expect_equal(names(out), c("C1", "C2", "Res"))
  expect_equal(nrow(out), nrow(two))
})
```

- [ ] **Step 2: Run the tests and verify they fail**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-io")'
```

Expected: FAIL with `could not find function "campaign_chems"`.

- [ ] **Step 3: Implement the four helpers**

Append to `R/app-io.R`:

```r
#' Which stressors a campaign frame actually doses
#'
#' A concentration column that is present but never positive (e.g. a `C3` of
#' zeros pasted in by mistake) does not count — the campaign is then a
#' two-stressor one.
#' @param df Campaign data frame with `C1`, `C2` and optionally `C3`.
#' @return Integer vector of stressor indices, e.g. `c(1L, 2L)`.
#' @keywords internal
campaign_chems <- function(df) {
  cols <- intersect(c("C1", "C2", "C3"), names(df))
  act  <- cols[vapply(df[cols], function(x) any(x > 0, na.rm = TRUE), logical(1))]
  as.integer(sub("^C", "", act))
}

#' Number of stressors in a campaign frame
#' @inheritParams campaign_chems
#' @return 2 or 3 (or fewer, which [validate_upload()] rejects).
#' @keywords internal
campaign_n_chem <- function(df) length(campaign_chems(df))

#' One stressor's single-stressor series from a campaign frame
#'
#' Keeps the rows where every OTHER stressor is 0 (so the shared control row is
#' included), drops their columns, and renames this stressor's concentration
#' column to `C1` — the shape a single-stressor fitter expects. Generalises
#' [marginal_df()] / [marginal_df3()] to 2- or 3-stressor frames.
#' @inheritParams campaign_chems
#' @param chem 1, 2 or 3 — which stressor's series to extract.
#' @return A data frame with `C1` and the response columns.
#' @keywords internal
single_df <- function(df, chem) {
  cols  <- intersect(c("C1", "C2", "C3"), names(df))
  this  <- paste0("C", chem)
  other <- setdiff(cols, this)
  keep  <- if (length(other))
    Reduce(`&`, lapply(other, function(k) df[[k]] == 0)) else rep(TRUE, nrow(df))
  out <- df[keep, , drop = FALSE]
  out[other] <- NULL
  names(out)[names(out) == this] <- "C1"
  rownames(out) <- NULL
  out
}

#' One pair's rows from a campaign frame
#'
#' Keeps the rows where the third stressor is 0 (so the singles and the shared
#' control come along, exactly as the binary fitter expects), drops its column,
#' and renames the pair's concentration columns to `C1`/`C2`.
#' @inheritParams campaign_chems
#' @param i,j Stressor indices of the pair, `i < j`.
#' @return A data frame with `C1`, `C2` and the response columns.
#' @keywords internal
pair_df <- function(df, i, j) {
  cols  <- intersect(c("C1", "C2", "C3"), names(df))
  this  <- paste0("C", c(i, j))
  other <- setdiff(cols, this)
  keep  <- if (length(other))
    Reduce(`&`, lapply(other, function(k) df[[k]] == 0)) else rep(TRUE, nrow(df))
  out <- df[keep, , drop = FALSE]
  out[other] <- NULL
  names(out)[names(out) == this[1]] <- "C1"
  names(out)[names(out) == this[2]] <- "C2"
  rownames(out) <- NULL
  out
}
```

- [ ] **Step 4: Run the tests and verify they pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-io")'
```

Expected: PASS.

- [ ] **Step 5: Add the `classify_rows` two-stressor characterization test**

`classify_rows()` (`R/fit-ternary-asa.R:5-12`) intersects against `C1`/`C2`/`C3` and counts positives, so it already behaves correctly on a two-stressor frame. This test pins that behaviour so a later edit cannot silently break the campaign. Append to `tests/testthat/test-fit-ternary-asa.R`:

```r
test_that("classify_rows handles a two-stressor frame: no ternary class", {
  two <- data.frame(C1 = c(0, 1, 0, 1), C2 = c(0, 0, 1, 1), Res = c(100, 80, 85, 60))
  cls <- classify_rows(two)
  expect_equal(as.character(cls), c("control", "single", "single", "binary"))
  expect_false("ternary" %in% as.character(cls))
})
```

- [ ] **Step 6: Run it and verify it passes**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "fit-ternary-asa")'
```

Expected: PASS (no source change needed — if it fails, fix `classify_rows` to match).

- [ ] **Step 7: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R tests/testthat/test-fit-ternary-asa.R
git commit -m "feat(io): campaign slicing primitives (single_df, pair_df, campaign_chems)"
```

---

## Task 2: Campaign upload schema, template and validation

**Files:**
- Modify: `R/app-io.R:9-19` (`upload_schema`), `:29-66` (`template_df`), `:74-118` (`validate_upload`)
- Test: `tests/testthat/test-app-io.R`

**Interfaces:**
- Consumes: `campaign_chems()`, `campaign_n_chem()`, `single_df()` from Task 1.
- Produces: `upload_schema("campaign", response)`, `template_df("campaign", response)`, and a `"campaign"` branch of `validate_upload(df, stage, response)` returning a character vector of messages (empty when valid).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-app-io.R`:

```r
test_that("campaign schema requires C1/C2 and treats C3 as optional", {
  expect_equal(upload_schema("campaign", "continuous"), c("C1", "C2", "Res"))
  expect_equal(upload_schema("campaign", "quantal"),
               c("C1", "C2", "Affected", "Exposed"))
})

test_that("the campaign template is a full three-stressor campaign", {
  tpl <- template_df("campaign", "continuous")
  expect_true(all(c("C1", "C2", "C3", "Res") %in% names(tpl)))
  cls <- classify_rows(tpl)
  expect_true(all(c("control", "single", "binary", "ternary") %in%
                  as.character(cls)))
})

test_that("a well-formed campaign validates clean", {
  expect_equal(validate_upload(campaign_fixture(), "campaign", "continuous"),
               character(0))
})

test_that("a campaign needs at least two dosed stressors", {
  df <- campaign_fixture()
  df$C2 <- 0
  df$C3 <- 0
  errs <- validate_upload(df, "campaign", "continuous")
  expect_true(any(grepl("at least two stressors", errs)))
})

test_that("each single series needs 4 distinct concentrations", {
  df <- campaign_fixture()
  df <- df[!(df$C2 > 0 & df$C1 == 0 & df$C3 == 0 & df$C2 > 2), ]  # thin chem 2
  errs <- validate_upload(df, "campaign", "continuous")
  expect_true(any(grepl("Stressor 2", errs)))
  expect_true(any(grepl("distinct concentrations", errs)))
})

test_that("a pair with no mixture rows is NOT an error (the sub-tab disables instead)", {
  df <- campaign_fixture()
  df <- df[!(df$C2 > 0 & df$C3 > 0), ]     # drop every 2x3 mixture row
  expect_equal(validate_upload(df, "campaign", "continuous"), character(0))
})
```

- [ ] **Step 2: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-io")'
```

Expected: FAIL — `'arg' should be one of "single", "binary", "ternary"`.

- [ ] **Step 3: Add "campaign" to the three schema functions**

In `R/app-io.R`, change the `match.arg` in `upload_schema` and add the branch:

```r
upload_schema <- function(stage, response) {
  stage <- match.arg(stage, c("single", "binary", "ternary", "campaign"))
  response <- match.arg(response, c("continuous", "quantal"))
  if (stage == "single") {
    if (response == "continuous") c("Conc", "Res") else c("Conc", "Affected", "Exposed")
  } else if (stage == "binary" || stage == "campaign") {
    # A campaign REQUIRES C1/C2; C3 is optional and makes it a three-stressor
    # campaign. Optionality is enforced in validate_upload(), not the schema.
    if (response == "continuous") c("C1", "C2", "Res") else c("C1", "C2", "Affected", "Exposed")
  } else {
    if (response == "continuous") c("C1", "C2", "C3", "Res")
    else c("C1", "C2", "C3", "Affected", "Exposed")
  }
}
```

In `template_df`, the ternary template is already a campaign file, so delegate rather than duplicate. Change its opening lines to:

```r
template_df <- function(stage, response) {
  # A campaign template IS the ternary template: all four strata in one frame.
  if (stage == "campaign") stage <- "ternary"
  cols <- upload_schema(stage, response)
  ...
```

- [ ] **Step 4: Add the campaign validation branch**

In `validate_upload`, change the ternary guard so it does not fire for campaigns, then append the campaign branch just before `errs`:

```r
  if (stage == "campaign") {
    cols <- intersect(c("C1", "C2", "C3"), present)
    ok   <- cols[vapply(df[cols], is.numeric, logical(1))]
    if (length(ok) >= 2) {
      chems <- campaign_chems(df[ok])
      if (length(chems) < 2) {
        errs <- c(errs, paste0(
          "A campaign needs at least two stressors with a positive ",
          "concentration; found ", length(chems), "."))
      } else {
        for (k in chems) {
          nd <- length(unique(stats::na.omit(single_df(df, k)$C1)))
          if (nd < 4)
            errs <- c(errs, paste0(
              "Stressor ", k, " has only ", nd, " distinct concentrations in ",
              "its single-stressor series; at least 4 are needed to fit a curve."))
        }
      }
    }
  }

  errs
}
```

Note: a pair with no mixture rows is deliberately **not** an error — the sub-tab disables itself (Task 6). The ternary tab's existing hard "No ternary rows" error stays as-is for `stage == "ternary"`.

- [ ] **Step 5: Run and verify pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-io")'
```

Expected: PASS, including the pre-existing single/binary/ternary tests.

- [ ] **Step 6: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat(io): campaign upload schema, template and validation"
```

---

## Task 3: Campaign module — store, upload, sub-navigation

**Files:**
- Create: `R/app-campaign.R`
- Modify: `R/app-run.R:6-29` (navbar + server wiring)
- Test: `tests/testthat/test-app-modules.R`

**`R/app-intro.R` needs no change**, though the spec's module table lists it. The
campaign store is *built on* `meta` (`store <- meta`), so `intro_server()`'s
existing writes already land in it. Leaving it untouched is what keeps
`axis_label()` and the curve-fit panels working unchanged.

**Interfaces:**
- Consumes: `read_upload()`, `validate_upload()`, `to_engine_df()`, `template_df()`, `campaign_chems()`, `campaign_n_chem()`.
- Produces:
  - `campaign_ui(id)` / `campaign_server(id, meta)`
  - the **campaign store**: a `shiny::reactiveValues` created by `campaign_server` and passed to every stage module, with fields `raw`, `n_chem`, `chems`, `response`, `reference`, `singles` (named list), `base_version` (integer), `pairs` (named list), `ternary`, plus the retained flat `chem1..3` / `unit1..3` / `expID` / `species` / `endpoint` fields written by `intro_server`.
  - `campaign_bump_base(store)` — increments `store$base_version`.

Later tasks call the store `store` and never mutate `base_version` directly except through `campaign_bump_base()`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-app-modules.R`:

```r
test_that("campaign_server detects stressor count from the upload", {
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_equal(store$n_chem, 3L)          # bundled example is a 3-stressor campaign
    expect_equal(store$chems, c(1L, 2L, 3L))
    expect_equal(store$reference, "CA")
  })
})

test_that("changing the campaign reference invalidates every downstream fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    v0 <- store$base_version
    store$singles <- list(`1` = list(par = c(max = 1, slope = 1, ec50 = 1)))

    session$setInputs(reference = "IA")
    expect_true(store$base_version > v0)
    expect_equal(length(store$singles), 0)   # fits cleared, not silently reused
  })
})

test_that("campaign_bump_base increments the stamp downstream stages compare against", {
  store <- shiny::reactiveValues(base_version = 0L)
  campaign_bump_base(store)
  campaign_bump_base(store)
  expect_equal(store$base_version, 2L)
})

test_that("campaign_ui builds a Shiny UI fragment", {
  ui <- campaign_ui("camp")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})
```

- [ ] **Step 2: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: FAIL with `could not find function "campaign_server"`.

- [ ] **Step 3: Create `R/app-campaign.R`**

```r
# Campaign stage: ONE upload covering 2-3 stressors, one campaign-wide response
# and reference choice, and the store every downstream stage reads. The stages
# (Singles / Binary pairs / Ternary) are sub-tabs of this one navbar entry --
# Sam's requested layout (issue #10). Single-stressor curves are fitted once on
# the Singles page and held fixed everywhere below it.

#' Increment the base-parameter version stamp
#'
#' Downstream stages record the `base_version` they were fitted at; when the
#' stamp no longer matches they show a stale banner instead of silently
#' presenting numbers computed from a superseded curve. Call this whenever a
#' single-stressor fit changes.
#' @param store The campaign store (a [shiny::reactiveValues()]).
#' @return Invisibly, the new version.
#' @keywords internal
campaign_bump_base <- function(store) {
  store$base_version <- (store$base_version %||% 0L) + 1L
  invisible(store$base_version)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Campaign stage UI
#' @param id Module id.
#' @keywords internal
campaign_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::helpText(shiny::tags$small(
        "One file for the whole campaign: single-stressor rows carry zeros in ",
        "the other concentrations, pair rows carry zero in the third.")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload campaign CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example campaign (FBSA + CPF + IMI, continuous) is loaded until ",
        "you upload your own.")),
      shiny::uiOutput(ns("errors")),
      shiny::uiOutput(ns("summary"))
    ),
    shiny::uiOutput(ns("stages"))
  )
}

#' Campaign stage server
#' @param id Module id.
#' @param meta Shared reactiveValues; the campaign store is built on it so the
#'   flat `chem1..3`/`unit1..3` fields written by [intro_server()] stay
#'   available to [axis_label()] and the curve-fit panels.
#' @return The campaign store, invisibly.
#' @keywords internal
campaign_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    store <- meta
    store$base_version <- 0L
    store$singles <- list()
    store$pairs   <- list()

    output$template <- shiny::downloadHandler(
      filename = function() paste0("campaign_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("campaign", input$response), file,
                         row.names = FALSE)
    )

    # The bundled example stands in until the user uploads their own file.
    parsed <- shiny::reactive({
      if (is.null(input$file))
        read_upload(system.file("extdata",
                                "ternary_ca_fbsa_cpf_imi_continuous.csv",
                                package = "mixdra"))
      else read_upload(input$file$datapath)
    })

    errs <- shiny::reactive(
      validate_upload(parsed(), "campaign", input$response %||% "continuous"))

    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    shiny::observe({
      shiny::req(length(errs()) == 0)
      df <- to_engine_df(parsed(), "campaign")
      store$raw       <- df
      store$chems     <- campaign_chems(df)
      store$n_chem    <- campaign_n_chem(df)
      store$response  <- input$response
      store$reference <- input$reference
    })

    # A new file or a changed response/reference invalidates every fit below.
    shiny::observeEvent(list(input$file, input$response, input$reference), {
      store$singles <- list()
      store$pairs   <- list()
      store$ternary <- NULL
      campaign_bump_base(store)
    }, ignoreInit = TRUE)

    output$summary <- shiny::renderUI({
      shiny::req(store$raw)
      cls <- classify_rows(store$raw)
      shiny::tags$small(shiny::HTML(paste0(
        "<b>", store$n_chem, "</b> stressors &middot; ",
        sum(cls == "control"), " control, ", sum(cls == "single"), " single, ",
        sum(cls == "binary"), " pair, ", sum(cls == "ternary"), " ternary rows.")))
    })

    output$stages <- shiny::renderUI(campaign_stage_nav(session$ns, store))

    invisible(store)
  })
}
```

Add a placeholder `campaign_stage_nav()` that Tasks 4, 6 and 8 fill in:

```r
#' Sub-navigation for the campaign stages
#'
#' Built server-side because which sub-tabs exist depends on the data: the
#' Ternary sub-tab is absent for a two-stressor campaign, and pairs with no
#' mixture rows are disabled.
#' @param ns The module's namespace function.
#' @param store The campaign store.
#' @keywords internal
campaign_stage_nav <- function(ns, store) {
  bslib::navset_card_tab(
    bslib::nav_panel("Singles", singles_ui(ns("singles")))
  )
}
```

- [ ] **Step 4: Wire it into the navbar**

In `R/app-run.R`, replace the three mixture `nav_panel` lines and the server wiring:

```r
app_ui <- function() {
  bslib::page_navbar(
    title = "Mixture Toxicity (mixdra)",
    fillable = FALSE,
    bslib::nav_panel("Introduction", intro_ui("intro")),
    bslib::nav_panel("Single Stressor", single_ui("single")),
    bslib::nav_panel("Campaign", campaign_ui("campaign")),
    # Retired in Task 9, kept here so the app stays usable mid-refactor.
    bslib::nav_panel("Binary Mixture", binary_ui("binary")),
    bslib::nav_panel("Ternary Mixture", ternary_ui("ternary"))
  )
}

app_server <- function(input, output, session) {
  meta <- shiny::reactiveValues()
  intro_server("intro", meta)
  single_server("single", meta)
  campaign_server("campaign", meta)
  binary_server("binary", meta)
  ternary_server("ternary", meta)
}
```

- [ ] **Step 5: Run the tests and verify they pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-run")'
```

Expected: PASS. (Task 4 creates `singles_ui`; until then `campaign_stage_nav` will error when rendered — the tests above do not render `output$stages`, so they pass. If `test-app-run.R` renders the full UI, complete Task 4 before re-running it.)

- [ ] **Step 6: Commit**

```bash
git add R/app-campaign.R R/app-run.R tests/testthat/test-app-modules.R
git commit -m "feat(app): campaign module - one upload, shared store, stage sub-nav"
```

---

## Task 4: Singles page — fit every curve once

**Files:**
- Create: `R/app-singles.R`
- Modify: `R/app-campaign.R` (`campaign_stage_nav`, wire `singles_server`)
- Test: `tests/testthat/test-app-modules.R`

**Interfaces:**
- Consumes: the store from Task 3; `single_df()` from Task 1; `curve_fit_ui()` / `curve_fit_server()` unchanged; `assemble_curve_params()` / `assemble_curve_params3()`.
- Produces:
  - `singles_ui(id)` / `singles_server(id, store)`
  - `campaign_base(store) -> named numeric or NULL` — the frozen base vector (`max`, `slope1..n`, `ec50*`) once every single is fitted, `NULL` otherwise. Downstream stages call this, never `assemble_curve_params*()` directly.

Note the naming asymmetry the engine already has: the binary base uses `ec501`/`ec502`, the ternary base uses `ec50_1`/`ec50_2`/`ec50_3`. `campaign_base()` returns whichever matches `store$n_chem`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("campaign_base returns NULL until every single is fitted", {
  store <- shiny::reactiveValues(n_chem = 3L, chems = c(1L, 2L, 3L),
                                 singles = list())
  expect_null(campaign_base(store))

  fake <- function(m, s, e) list(par = c(max = m, slope = s, ec50 = e))
  store$singles <- list(`1` = fake(100, 2, 1), `2` = fake(100, 3, 2))
  expect_null(campaign_base(store))

  store$singles <- list(`1` = fake(100, 2, 1), `2` = fake(100, 3, 2),
                        `3` = fake(100, 4, 4))
  b <- campaign_base(store)
  expect_equal(sort(names(b)),
               sort(c("max", "slope1", "slope2", "slope3",
                      "ec50_1", "ec50_2", "ec50_3")))
  expect_equal(unname(b[["ec50_3"]]), 4)
  expect_equal(unname(b[["max"]]), 100)
})

test_that("campaign_base uses the binary parameter names for a 2-stressor campaign", {
  fake <- function(m, s, e) list(par = c(max = m, slope = s, ec50 = e))
  store <- shiny::reactiveValues(
    n_chem = 2L, chems = c(1L, 2L),
    singles = list(`1` = fake(100, 2, 1), `2` = fake(100, 3, 2)))
  expect_equal(sort(names(campaign_base(store))),
               sort(c("max", "slope1", "slope2", "ec501", "ec502")))
})

test_that("singles_server starts with no fits and does not invent one", {
  # The store-write and version-bump behaviour is covered by the
  # campaign_bump_base unit test in Task 3; driving nested curve_fit_server
  # modules through testServer is not worth the harness cost here. The Step 6
  # manual app check is what confirms a fit reaches the store.
  store <- shiny::reactiveValues(
    n_chem = 2L, chems = c(1L, 2L), base_version = 0L, singles = list(),
    raw = data.frame(C1 = c(0, 1, 2, 3, 0, 0, 0),
                     C2 = c(0, 0, 0, 0, 1, 2, 3),
                     Res = c(100, 80, 60, 40, 90, 70, 50)))
  shiny::testServer(singles_server, args = list(store = store), {
    expect_equal(length(store$singles), 0)
    expect_equal(store$base_version, 0L)
  })
})

test_that("singles_ui builds a Shiny UI fragment", {
  ui <- singles_ui("s")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})
```

- [ ] **Step 2: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: FAIL with `could not find function "campaign_base"`.

- [ ] **Step 3: Create `R/app-singles.R`**

```r
# Singles stage: every stressor's dose-response curve, fitted ONCE, on one page.
# These curves are the campaign's single source of truth -- every pair workspace
# and the ternary stage hold them fixed. Each stressor gets its own curve_fit
# panel, reused verbatim from the standalone Single Stressor tab.

#' The campaign's frozen base-parameter vector
#'
#' `NULL` until every stressor in the campaign has a fitted curve. Returns the
#' binary parameter names (`ec501`/`ec502`) for a two-stressor campaign and the
#' ternary ones (`ec50_1`..`ec50_3`) for a three-stressor campaign, matching the
#' registries the downstream fits use.
#' @param store The campaign store.
#' @return A named numeric vector, or `NULL`.
#' @keywords internal
campaign_base <- function(store) {
  chems <- store$chems
  if (is.null(chems)) return(NULL)
  fits <- lapply(as.character(chems), function(k) store$singles[[k]])
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)
  if (length(fits) == 3)
    assemble_curve_params3(fits[[1]], fits[[2]], fits[[3]])
  else
    assemble_curve_params(fits[[1]], fits[[2]])
}

#' Singles stage UI
#' @param id Module id.
#' @keywords internal
singles_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header("Single-stressor curves"),
    shiny::p("Fit each stressor's dose-response curve (Autofit or Simulate). ",
             "These curves are fitted once here and held fixed by every stage ",
             "below, so one stressor has exactly one EC50 for the whole campaign."),
    shiny::uiOutput(ns("panels"))
  )
}

#' Singles stage server
#' @param id Module id.
#' @param store The campaign store.
#' @keywords internal
singles_server <- function(id, store) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$panels <- shiny::renderUI({
      shiny::req(store$chems)
      shiny::tagList(lapply(store$chems, function(k) {
        shiny::div(
          class = if (k == store$chems[1]) "" else "mt-4",
          shiny::h5(axis_label(store, paste0("chem", k))),
          curve_fit_ui(ns(paste0("chem", k))))
      }))
    })

    # One curve_fit_server per stressor. Registered for all three slots up front
    # (a module server cannot be created inside renderUI); slots the campaign
    # does not use simply never receive data.
    for (k in 1:3) local({
      kk <- k
      fit <- curve_fit_server(
        paste0("chem", kk),
        fit_df = shiny::reactive({
          shiny::req(store$raw, kk %in% store$chems)
          single_df(store$raw, kk)
        }),
        meta = store, chem_field = paste0("chem", kk))

      shiny::observeEvent(fit(), {
        store$singles[[as.character(kk)]] <- fit()
        campaign_bump_base(store)
      }, ignoreNULL = TRUE)
    })
  })
}
```

- [ ] **Step 4: Wire it into the campaign sub-nav**

In `R/app-campaign.R`, call `singles_server` inside `campaign_server` (just before `invisible(store)`):

```r
    singles_server("singles", store)
```

- [ ] **Step 5: Run and verify pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: PASS.

- [ ] **Step 6: Launch the app and confirm the page works end to end**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); run_app(launch.browser = TRUE)'
```

Confirm: the Campaign tab loads the bundled example, the sidebar summary reports 3 stressors with the four row counts, and three curve panels appear and each fits with Autofit.

- [ ] **Step 7: Commit**

```bash
git add R/app-singles.R R/app-campaign.R tests/testthat/test-app-modules.R
git commit -m "feat(app): singles stage - fit every stressor curve once per campaign"
```

**End of Phase 1 — this is the point to show Sam.**

---

# Phase 2 — The pair workspace

## Task 5: Extract the interaction workspace from the binary tab

**Files:**
- Modify: `R/app-binary.R` (split `binary_ui`/`binary_server`)
- Test: `tests/testthat/test-app-modules.R`

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `pair_workspace_ui(id)` — Stage 2 + Stage 3 cards plus the Advanced-fitting accordion (today `R/app-binary.R:141-185` and `:102-118`)
  - `pair_workspace_server(id, fit_df, base, reference, response, base_version = NULL, on_fit = NULL)`

    `fit_df`, `base`, `reference`, `response` and `base_version` are **reactives**; `base` returns the frozen named vector or `NULL`. **Declare all seven parameters now**, with `base_version` and `on_fit` defaulting to `NULL` and unused until Tasks 6 and 8 respectively — so the signature never changes after this task.
  - `binary_ui`/`binary_server` survive as thin wrappers so the existing tab keeps working until Task 9

This is a pure extraction: **move code, do not rewrite it.** The staged loop, the joint refine, the table, the plots and the selection logic are all correct and validated against Excel — changing their behaviour here would invalidate `test-joint-excel-xcheck.R`.

- [ ] **Step 1: Write the characterization test first**

This test pins today's numbers so the extraction cannot change them. Append to `tests/testthat/test-app-modules.R`:

```r
test_that("pair_workspace reproduces the binary tab's staged S/A fit", {
  skip_if_not_installed("shiny")
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  p12 <- pair_df(df, 1, 2)

  f1 <- analyse_single(single_df(df, 1))
  f2 <- analyse_single(single_df(df, 2))
  base <- assemble_curve_params(f1, f2)

  fit <- fit_model(p12, "CA", "SA", "continuous",
                   start = c(base, a = 0), fixed = names(base),
                   n_starts = 1, time_limit = 30)

  expect_true(is.finite(fit$par[["a"]]))
  expect_named(fit$par[names(base)], names(base))
  expect_equal(unname(fit$par[names(base)]), unname(base), tolerance = 1e-8)
})
```

- [ ] **Step 2: Run it and verify it passes against today's code**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: PASS. This is the baseline the extraction must preserve.

- [ ] **Step 3: Split the UI**

In `R/app-binary.R`, create `pair_workspace_ui(id)` containing, verbatim: the `bslib::accordion(...)` block currently at `:102-118`, then the Stage 2 card (`:141-171`) and the Stage 3 card (`:172-185`). Then reduce `binary_ui` to:

```r
binary_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example dataset (CPF + IMI, continuous) is loaded until you upload your own.")),
      shiny::uiOutput(ns("errors"))
    ),
    bslib::card(
      bslib::card_header("Stage 1 · Single curves"),
      shiny::p("Fit each stressor's dose-response curve (Autofit or Simulate). ",
               "These curves are held fixed when the interaction models are ",
               "compared below."),
      shiny::div(shiny::h5(shiny::textOutput(ns("chem1_title"))),
                 curve_fit_ui(ns("chem1"))),
      shiny::div(class = "mt-4",
                 shiny::h5(shiny::textOutput(ns("chem2_title"))),
                 curve_fit_ui(ns("chem2"))),
      shiny::uiOutput(ns("reveal_note"))
    ),
    pair_workspace_ui(ns("work"))
  )
}
```

- [ ] **Step 4: Split the server**

Create `pair_workspace_server(id, fit_df, base, reference, response)` holding everything currently in `binary_server` from `fits_store <- shiny::reactiveVal(list())` (`:259`) to the end of the module (`:627`), with these substitutions throughout the moved code:

- `engine_df()` → `fit_df()`
- `curve_params()` → `base()`
- `input$reference` → `reference()`
- `engine_response()` → `response()`
- `frozen()` → `!is.null(base())`

`binary_server` keeps the upload, the two `curve_fit_server` calls, `curve_params`, and delegates:

```r
    pair_workspace_server("work",
                          fit_df    = engine_df,
                          base      = shiny::reactive(
                            if (!is.null(fit1()) && !is.null(fit2()))
                              assemble_curve_params(fit1(), fit2()) else NULL),
                          reference = shiny::reactive(input$reference),
                          response  = engine_response)
```

- [ ] **Step 5: Verify nothing changed**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "joint-refine")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "interaction-help")'
```

Expected: PASS, all three, with the same counts as before the split.

- [ ] **Step 6: Confirm the line count came down**

```bash
wc -l R/app-binary.R
```

Expected: well under 400 lines for the `binary_ui`/`binary_server` portion; `pair_workspace_*` may live in the same file.

- [ ] **Step 7: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "refactor(app): extract pair_workspace from the binary tab

Pure extraction - the staged loop, joint refine, table, plots and
selection logic are moved verbatim. binary_ui/binary_server remain as
thin wrappers until the campaign replaces them."
```

---

## Task 6: Three pair sub-tabs on frozen curves

**Files:**
- Modify: `R/app-campaign.R` (`campaign_stage_nav`, `campaign_server`)
- Test: `tests/testthat/test-app-modules.R`

**Interfaces:**
- Consumes: `pair_workspace_ui`/`pair_workspace_server` (Task 5), `campaign_base()` (Task 4), `pair_df()` (Task 1).
- Produces:
  - `campaign_pairs(chems) -> list of c(i, j)` — the pairs of a campaign, in order: `list(c(1,2), c(1,3), c(2,3))` for three stressors, `list(c(1,2))` for two
  - `pair_key(i, j) -> character` — `"12"`, `"13"`, `"23"`
  - `pair_has_rows(df, i, j) -> logical` — whether that pair has any mixture row

- [ ] **Step 1: Write the failing tests**

```r
test_that("campaign_pairs enumerates the pairs in A-B, A-C, B-C order", {
  expect_equal(campaign_pairs(c(1L, 2L, 3L)),
               list(c(1L, 2L), c(1L, 3L), c(2L, 3L)))
  expect_equal(campaign_pairs(c(1L, 2L)), list(c(1L, 2L)))
  expect_equal(pair_key(2, 3), "23")
})

test_that("pair_has_rows spots a pair with no mixture rows", {
  df <- campaign_fixture()
  expect_true(pair_has_rows(df, 1, 2))
  expect_true(pair_has_rows(df, 2, 3))
  df2 <- df[!(df$C2 > 0 & df$C3 > 0), ]
  expect_false(pair_has_rows(df2, 2, 3))
})
```

`campaign_fixture()` is defined in `tests/testthat/test-app-io.R`; move it to `tests/testthat/helper-campaign.R` so both test files can use it.

- [ ] **Step 2: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: FAIL with `could not find function "campaign_pairs"`.

- [ ] **Step 3: Implement the three helpers in `R/app-campaign.R`**

```r
#' The pairs of a campaign, in A-B, A-C, B-C order
#' @param chems Integer vector of stressor indices.
#' @return A list of length-2 integer vectors.
#' @keywords internal
campaign_pairs <- function(chems) {
  if (length(chems) < 2) return(list())
  utils::combn(sort(chems), 2, simplify = FALSE)
}

#' Store key for a pair
#' @param i,j Stressor indices.
#' @return A character key, e.g. `"23"`.
#' @keywords internal
pair_key <- function(i, j) paste0(i, j)

#' Does this pair have any mixture rows?
#'
#' A campaign may cover a pair's singles without ever dosing them together; that
#' pair's sub-tab is disabled rather than treated as an upload error.
#' @param df Campaign engine frame.
#' @param i,j Stressor indices.
#' @return `TRUE` if at least one row has both concentrations positive.
#' @keywords internal
pair_has_rows <- function(df, i, j) {
  any(df[[paste0("C", i)]] > 0 & df[[paste0("C", j)]] > 0, na.rm = TRUE)
}
```

- [ ] **Step 4: Build the pair sub-tabs**

Replace `campaign_stage_nav()` with:

```r
campaign_stage_nav <- function(ns, store) {
  panels <- list(bslib::nav_panel("Singles", singles_ui(ns("singles"))))

  ready <- !is.null(campaign_base(store))
  for (p in campaign_pairs(store$chems)) {
    k     <- pair_key(p[1], p[2])
    title <- paste(axis_label(store, paste0("chem", p[1])), "×",
                   axis_label(store, paste0("chem", p[2])))
    body <- if (!pair_has_rows(store$raw, p[1], p[2])) {
      shiny::div(class = "p-3 text-muted",
                 "This campaign has no mixture rows for this pair.")
    } else if (!ready) {
      shiny::div(class = "p-3 text-muted",
                 sprintf("Fit all %d single-stressor curves first.",
                         store$n_chem))
    } else {
      pair_workspace_ui(ns(paste0("pair", k)))
    }
    panels <- c(panels, list(bslib::nav_panel(title, body)))
  }

  # INSERT-TERNARY-PANEL-HERE (Task 8 adds its block at this point, before the
  # do.call — appending after it would silently drop the panel).

  do.call(bslib::navset_card_tab, panels)
}
```

- [ ] **Step 5: Register one pair workspace per pair, with staleness stamping**

In `campaign_server`, after `singles_server("singles", store)`:

```r
    # All three slots are registered up front -- a module server cannot be
    # created inside renderUI. Unused slots never receive data.
    for (p in list(c(1, 2), c(1, 3), c(2, 3))) local({
      pp <- p
      k  <- pair_key(pp[1], pp[2])
      pair_workspace_server(
        paste0("pair", k),
        fit_df = shiny::reactive({
          shiny::req(store$raw, pair_has_rows(store$raw, pp[1], pp[2]))
          pair_df(store$raw, pp[1], pp[2])
        }),
        base      = shiny::reactive(campaign_base(store)),
        reference = shiny::reactive(store$reference),
        response  = shiny::reactive(store$response))
    })
```

`pair_workspace_server` already rebuilds its store when `base()` changes (it is a reactive dependency of the staged loop), so a refitted single invalidates the pair results through the normal reactive path. Record the stamp so the UI can say so: inside `pair_workspace_server`, when `fits_store()` is filled, also set `stamped_at(store_base_version)` — pass `base_version = shiny::reactive(store$base_version)` as a fifth argument and render a banner when `stamped_at() != base_version()`:

```r
    output$stale <- shiny::renderUI({
      shiny::req(length(fits_store()) > 0, stamped_at())
      if (!identical(stamped_at(), base_version()))
        shiny::div(class = "alert alert-warning",
                   "A single-stressor curve changed after these models were ",
                   "fitted. Re-run ", shiny::tags$b("Fit interaction models"),
                   " to bring them up to date.")
    })
```

Add `shiny::uiOutput(ns("stale"))` at the top of the Stage 2 card in `pair_workspace_ui`.

- [ ] **Step 6: Run and verify pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: PASS.

- [ ] **Step 7: Confirm in the app**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); run_app(launch.browser = TRUE)'
```

Confirm: three pair sub-tabs titled from the stressor names; each says "Fit all 3 single-stressor curves first" until the Singles page is done; each then fits its interaction models; refitting a single raises the stale banner on the pairs.

- [ ] **Step 8: Commit**

```bash
git add R/app-campaign.R tests/testthat/test-app-modules.R tests/testthat/helper-campaign.R
git commit -m "feat(app): binary pairs as campaign sub-tabs on frozen curves"
```

---

# Phase 3 — The ternary stage on frozen pairwise terms

## Task 7: Engine — `pairwise=` on the staged ASA fit

**Files:**
- Modify: `R/fit-ternary-asa.R:46-72` (`fit_ternary_asa`), `:142-152` (`analyse_ternary`)
- Test: `tests/testthat/test-fit-ternary-asa.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `fit_ternary_asa(df, reference, response, base = NULL, pairwise = NULL, lower, upper, n_starts, time_limit)` and the same new argument forwarded by `analyse_ternary()`. `pairwise` is a named numeric vector `c(A1 = …, A2 = …, A3 = …)`; when supplied, Stage 2 is skipped and those values are held fixed. `NULL` keeps today's behaviour exactly.

- [ ] **Step 1: Write the failing tests — including the linchpin**

Append to `tests/testthat/test-fit-ternary-asa.R`:

```r
test_that("supplying pairwise skips Stage 2 and holds A1/A2/A3 fixed", {
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  auto <- fit_ternary_asa(df, n_starts = 1)
  pw   <- auto$pairwise

  frozen <- fit_ternary_asa(df, base = auto$base, pairwise = pw, n_starts = 1)
  expect_equal(unname(frozen$pairwise), unname(pw), tolerance = 1e-10)
  expect_null(frozen$fits$pairwise)          # Stage 2 was not run
  expect_true(is.finite(frozen$A4_overall))
})

test_that("pairwise validates its names", {
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  expect_error(fit_ternary_asa(df, pairwise = c(A1 = 0, A2 = 0)),
               "missing param")
  expect_error(fit_ternary_asa(df, pairwise = list(A1 = 0, A2 = 0, A3 = 0)),
               "named numeric")
})

test_that("LINCHPIN: each pair's S/A a IS its ternary A-term at the same base", {
  # Spec section 3. If this fails, reusing the fitted binaries in the ternary
  # stage is no longer exact and the campaign design needs revisiting.
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  auto <- fit_ternary_asa(df, n_starts = 1)
  base <- auto$base

  # The binary base uses ec501/ec502; map the ternary base onto each pair.
  pair_base <- function(i, j) c(
    max    = unname(base[["max"]]),
    slope1 = unname(base[[paste0("slope", i)]]),
    slope2 = unname(base[[paste0("slope", j)]]),
    ec501  = unname(base[[paste0("ec50_", i)]]),
    ec502  = unname(base[[paste0("ec50_", j)]]))

  a_of <- function(i, j) {
    b <- pair_base(i, j)
    fit_model(pair_df(df, i, j), "CA", "SA", "continuous",
              start = c(b, a = 0), fixed = names(b),
              n_starts = 1, time_limit = 60)$par[["a"]]
  }

  expect_equal(a_of(1, 2), unname(auto$pairwise[["A1"]]), tolerance = 1e-3)
  expect_equal(a_of(1, 3), unname(auto$pairwise[["A2"]]), tolerance = 1e-3)
  expect_equal(a_of(2, 3), unname(auto$pairwise[["A3"]]), tolerance = 1e-3)
})
```

- [ ] **Step 2: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "fit-ternary-asa")'
```

Expected: FAIL — `unused argument (pairwise = pw)`.

- [ ] **Step 3: Add the argument**

In `R/fit-ternary-asa.R`, add `pairwise = NULL` to the signature after `base = NULL`, document it, and replace the Stage-2 block (`:74-79`) with:

```r
  # Stage 2 - A1/A2/A3 from binaries; base + A4 held fixed; start all A at 0.
  # A supplied `pairwise` (the campaign's per-pair S/A values) skips this fit and
  # holds those values fixed downstream -- exact, because each binary row
  # activates exactly one A term (see the design doc, section 3).
  if (is.null(pairwise)) {
    binaries <- df[cls == "binary", , drop = FALSE]
    f2 <- fit_model(binaries, reference, "ASA", response,
                    start = c(base_par, A1 = 0, A2 = 0, A3 = 0, A4 = 0),
                    fixed = c(base_params, "A4"),
                    n_starts = n_starts, time_limit = time_limit)
    a123 <- f2$par[c("A1", "A2", "A3")]
  } else {
    if (!is.numeric(pairwise))
      stop("fit_ternary_asa: `pairwise` must be a named numeric vector")
    miss <- setdiff(c("A1", "A2", "A3"), names(pairwise))
    if (length(miss))
      stop("fit_ternary_asa: `pairwise` is missing param(s): ",
           paste(miss, collapse = ", "))
    a123 <- pairwise[c("A1", "A2", "A3")]
    f2 <- NULL
  }
```

In `analyse_ternary` (`:142-152`), add `pairwise = NULL` to the signature and pass it through to `fit_ternary_asa`. Document it with the same wording as `base`.

- [ ] **Step 4: Run and verify pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "fit-ternary-asa")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "validation")'
```

Expected: PASS. If the linchpin test fails, **stop and report** — do not loosen the tolerance to make it pass.

- [ ] **Step 5: Regenerate docs**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::document()'
```

- [ ] **Step 6: Commit**

```bash
git add R/fit-ternary-asa.R man tests/testthat/test-fit-ternary-asa.R
git commit -m "feat(engine): fit_ternary_asa(pairwise=) to reuse fitted binaries

Skips Stage 2 and holds A1/A2/A3 fixed at supplied values. Exact,
because each binary row activates exactly one A term - guarded by a
test that each pair's S/A a equals its ternary A-term at the same base."
```

---

## Task 8: Ternary sub-tab on frozen base + pairwise

**Files:**
- Modify: `R/app-ternary.R` (drop its own upload and Stages 1–2; consume the store)
- Modify: `R/app-campaign.R` (`campaign_stage_nav`, `campaign_server`)
- Test: `tests/testthat/test-app-modules.R`

**Interfaces:**
- Consumes: `campaign_base()`, `campaign_pairs()`, `pair_key()`, `analyse_ternary(pairwise=)`.
- Produces: `campaign_pairwise(store) -> named numeric or NULL` — `c(A1, A2, A3)` assembled from the three pair workspaces' selected S/A fits, `NULL` until all three exist.

- [ ] **Step 1: Write the failing test**

```r
test_that("campaign_pairwise maps pair S/A values onto A1/A2/A3", {
  store <- shiny::reactiveValues(
    chems = c(1L, 2L, 3L),
    pairs = list(`12` = list(a = 0.5), `13` = list(a = -0.2)))
  expect_null(campaign_pairwise(store))

  store$pairs$`23` <- list(a = 1.1)
  pw <- campaign_pairwise(store)
  expect_equal(unname(pw[["A1"]]), 0.5)    # pair 1-2
  expect_equal(unname(pw[["A2"]]), -0.2)   # pair 1-3
  expect_equal(unname(pw[["A3"]]), 1.1)    # pair 2-3
})
```

- [ ] **Step 2: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: FAIL with `could not find function "campaign_pairwise"`.

- [ ] **Step 3: Implement `campaign_pairwise` in `R/app-campaign.R`**

```r
#' The campaign's frozen pairwise interaction terms
#'
#' `A1` is the 1-2 pair, `A2` the 1-3 pair, `A3` the 2-3 pair -- the term order
#' `ca_asa_tri()` uses. Each value is that pair's fitted S/A `a`, which is the
#' same quantity (see the design doc, section 3). `NULL` until all three pairs
#' have been fitted.
#' @param store The campaign store.
#' @return A named numeric `c(A1, A2, A3)`, or `NULL`.
#' @keywords internal
campaign_pairwise <- function(store) {
  if (length(store$chems) != 3) return(NULL)
  keys <- c("12", "13", "23")
  vals <- lapply(keys, function(k) store$pairs[[k]]$a)
  if (any(vapply(vals, is.null, logical(1)))) return(NULL)
  stats::setNames(as.numeric(unlist(vals)), c("A1", "A2", "A3"))
}
```

- [ ] **Step 4: Have each pair workspace publish its S/A value**

In `pair_workspace_server`, add a `on_fit` callback argument (default `NULL`) invoked whenever the staged loop finishes, and call it from `campaign_server` to record the pair's S/A `a` and the selected model:

```r
      if (!is.null(on_fit))
        on_fit(list(a = if (!is.null(s[["SA"]])) s[["SA"]]$par[["a"]] else NULL,
                    chosen = cmp$chosen,
                    base_version = base_version()))
```

In `campaign_server`, inside the pair registration loop:

```r
        on_fit = function(res) store$pairs[[k]] <- res
```

- [ ] **Step 5: Convert the ternary module to a store consumer**

In `R/app-ternary.R`: delete the sidebar upload, template download and validation; delete the three Stage-1 `curve_fit_server` calls and the `assemble_curve_params3` reactive. Change `ternary_server(id, meta)` to `ternary_server(id, store)` and make the fit:

```r
    asa_res <- shiny::reactive({
      shiny::req(store$raw, store$reference == "CA",
                 store$response == "continuous")
      base <- campaign_base(store); shiny::req(base)
      pw   <- campaign_pairwise(store); shiny::req(pw)
      analyse_ternary(store$raw, reference = "CA", response = "continuous",
                      base = base, pairwise = pw,
                      n_starts = input$n_starts %||% 1)
    })
```

- [ ] **Step 6: Add the ternary sub-tab with its gates**

In `campaign_stage_nav()`, after the pair panels:

```r
  if (length(store$chems) == 3) {
    body <- if (!identical(store$reference, "CA")) {
      shiny::div(class = "p-3 text-muted",
                 "The ternary stage is not yet supported for Independent ",
                 "Action. Switch the campaign reference model to Concentration ",
                 "Addition to use it.")
    } else if (!identical(store$response, "continuous")) {
      shiny::div(class = "p-3 text-muted",
                 "The ternary stage is not yet supported for quantal data.")
    } else if (is.null(campaign_pairwise(store))) {
      shiny::div(class = "p-3 text-muted",
                 "Fit the interaction models on all three pair tabs first.")
    } else {
      ternary_ui(ns("ternary"))
    }
    panels <- c(panels, list(bslib::nav_panel("Ternary", body)))
  }
```

Insert this block at the `# INSERT-TERNARY-PANEL-HERE` marker, **before** the `do.call`. The two-stressor case needs no branch: the panel is simply never added.

Add a note listing any pair whose `chosen` is not `"SA"`. **It belongs to `ternary_server`** (it describes the ternary fit), so add `shiny::uiOutput(ns("sa_note"))` near the top of `ternary_ui` and this renderer to `ternary_server`:

```r
    output$sa_note <- shiny::renderUI({
      odd <- names(Filter(function(p) !identical(p$chosen, "SA"), store$pairs))
      if (length(odd))
        shiny::div(class = "alert alert-info",
                   "Pair(s) ", paste(odd, collapse = ", "),
                   " are best described by a dose-ratio or dose-level model. ",
                   "The ternary fit uses their S/A term, which is the only form ",
                   "the three-way model has.")
    })
```

- [ ] **Step 7: Wire the ternary server in `campaign_server`**

```r
    ternary_server("ternary", store)
```

- [ ] **Step 8: Run and verify pass**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "plot-ternary")'
```

Expected: PASS. Existing `ternary_server` smoke tests will need their `args` updated from `meta =` to `store =`.

- [ ] **Step 9: Confirm the whole flow in the app**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::load_all("."); run_app(launch.browser = TRUE)'
```

Confirm end to end: upload → three singles → three pairs → ternary, with the isoplane rendering and the ternary using the frozen values (its `base` and `pairwise` should equal what the earlier stages show).

- [ ] **Step 10: Commit**

```bash
git add R/app-ternary.R R/app-campaign.R tests/testthat/test-app-modules.R
git commit -m "feat(app): ternary stage on frozen singles and frozen pairwise terms"
```

---

# Phase 4 — Retire the old tabs and document

## Task 9: Remove the standalone mixture tabs

**Files:**
- Modify: `R/app-run.R` (drop two nav entries and two server calls)
- Modify: `R/app-binary.R` (delete `binary_ui`/`binary_server`; keep `interaction_help` and `pair_workspace_*`)
- Modify: `R/app-io.R` (delete `marginal_df`, `marginal_df3` — now unused)
- Modify: `tests/testthat/test-app-io.R`, `test-app-modules.R`, `test-app-run.R`
- Test: as above

**Interfaces:**
- Consumes: everything from Phases 1–3.
- Produces: a navbar of exactly Introduction / Single Stressor / Campaign.

- [ ] **Step 1: Confirm nothing still calls the doomed functions**

```bash
grep -rn "binary_ui\|binary_server\|marginal_df3\|marginal_df\b" R/ tests/ vignettes/
```

Expected: hits only in the files listed above and in tests that are about to be updated. If anything else appears, stop and reassess.

- [ ] **Step 2: Update the tests first**

Delete the `binary_server` / `binary_ui` / `marginal_df` / `marginal_df3` tests from `tests/testthat/test-app-io.R` and `test-app-modules.R`. Update `test-app-run.R` to expect three nav panels:

```r
test_that("app_ui exposes Introduction, Single Stressor and Campaign", {
  ui <- app_ui()
  txt <- paste(as.character(ui), collapse = " ")
  expect_true(grepl("Introduction", txt))
  expect_true(grepl("Single Stressor", txt))
  expect_true(grepl("Campaign", txt))
  expect_false(grepl("Binary Mixture", txt))
  expect_false(grepl("Ternary Mixture", txt))
})
```

- [ ] **Step 3: Run and verify failure**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-run")'
```

Expected: FAIL — "Binary Mixture" still present.

- [ ] **Step 4: Delete the retired code**

In `R/app-run.R` remove the two `nav_panel` lines and the two `*_server` calls. In `R/app-binary.R` delete `binary_ui` and `binary_server`. In `R/app-io.R` delete `marginal_df` and `marginal_df3`.

- [ ] **Step 5: Run the affected filters**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-run")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-io")'
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-modules")'
```

Expected: PASS.

- [ ] **Step 6: Regenerate docs and commit**

```bash
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::document()'
git add R/ man/ tests/
git commit -m "refactor(app): retire the standalone Binary and Ternary tabs"
```

---

## Task 10: Documentation

**Files:**
- Modify: `README.md` (the `run_app()` section)
- Modify: `vignettes/methodology.Rmd.orig` and the shipped `vignettes/methodology.Rmd`
- Modify: `NEWS.md` if present

- [ ] **Step 1: Update the README `run_app()` section**

Replace `README.md:97-99` (the "Opens a `bslib` dashboard with **Single Chemical**, **Binary Mixture**, and **Ternary Mixture** workflows…" paragraph) with:

```markdown
Opens a `bslib` dashboard with three tabs: **Introduction** for experiment
metadata, **Single Stressor** for a standalone one-stressor experiment, and
**Campaign** for a full mixture study.

A campaign is one uploaded file covering two or three stressors — single-stressor
rows carry zeros in the other concentrations, pair rows carry zero in the third.
Its sub-tabs run the study in order:

1. **Singles** — every stressor's dose-response curve, fitted once. These curves
   are the campaign's single source of truth.
2. **One sub-tab per pair** — the interaction model comparison
   (reference → S/A → DR/DL) for each pair, with the single-stressor curves held
   fixed so interaction shows up as `a`/`b` rather than being absorbed by
   refitted curves.
3. **Ternary** — the staged Advanced-S/A fit, reusing both the frozen curves and
   the pairwise terms already fitted above. Present only for a three-stressor
   campaign, and currently CA + continuous only.

The Single Stressor tab is an independent scratchpad: it does not read from or
write to the campaign.
```

- [ ] **Step 2: Check the vignette for stale app references**

```bash
grep -n "Binary Mixture\|Ternary Mixture\|tab" vignettes/methodology.Rmd.orig
```

Update any prose that describes the old tab layout. **Do not re-knit** — the vignette is pre-computed (`vignettes/precompute.R`) and re-running it re-runs every fit. Only edit prose that appears in both `.Rmd.orig` and the shipped `.Rmd`; if a change would alter a figure, stop and flag it.

- [ ] **Step 3: Commit**

```bash
git add README.md vignettes/ NEWS.md
git commit -m "docs: describe the campaign workflow"
```

- [ ] **Step 4: Push and open the PR**

```bash
git push -u origin campaign-workflow
gh pr create --base main --title "Campaign workflow: one upload, singles fitted once, binaries as sub-tabs (#10)" --body "$(cat <<'BODY'
Implements the campaign workflow from issue #10.

- One upload per campaign (2 or 3 stressors); single-stressor curves fitted once on a Singles page and held fixed by every stage below.
- Binary pairs become sub-tabs of the Campaign tab, sliced from the same file.
- The ternary stage consumes both the frozen curves and the frozen pairwise terms.

Design: `docs/superpowers/specs/2026-08-20-campaign-workflow-design.md`
Plan: `docs/superpowers/plans/2026-08-21-campaign-workflow.md`

Closes #10.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Verification checklist before requesting review

- [ ] `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter = "app-io")'` passes
- [ ] Same for filters `app-modules`, `app-run`, `fit-ternary-asa`, `joint-refine`, `validation`
- [ ] The linchpin test in Task 7 passes at its original tolerance
- [ ] `run_app()` walks upload → 3 singles → 3 pairs → ternary without an error
- [ ] A two-stressor campaign file shows one pair sub-tab and no Ternary sub-tab
- [ ] Refitting a single raises the stale banner on every downstream pair
- [ ] No `Co-Authored-By: Claude` trailer on any commit in the branch
