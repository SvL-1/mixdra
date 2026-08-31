# mixdra Shiny App (single + binary) Implementation Plan

**Goal:** Add an interactive `bslib` Shiny front-end (Introduction, Single Chemical, Binary Mixture) to the `mixdra` package, launched by an exported `run_app()`, built entirely on the finished fitting engine and `plotly` plotting layer.

**Architecture:** Shiny modules live in `R/` so they are unit-testable with `shiny::testServer`; pure I/O helpers (template/validation/read) are plotly- and Shiny-free and carry the bulk of the test coverage. The app uses `bslib::page_navbar` with a `layout_sidebar` per analysis stage. All UI dependencies are `Suggests`, guarded by `run_app()`, so the engine stays lean and `R CMD check`-clean.

**Tech Stack:** R package; `shiny`, `bslib`, `plotly`, `DT` (all Suggests); `testthat` (edition 3); `roxygen2`.

**Design spec:** `docs/design/specs/2026-05-31-mixdra-shiny-app-design.md`

---

## Running R and tests here (read first)

- **In RStudio:** open `mixdra.Rproj`, `Ctrl+Shift+L` (load_all), then run each task's test command in the Console.
- **From the CLI:** every non-interactive R call must prepend the user-library path or packages are not found:
  ```r
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  ```
  e.g. `Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); devtools::test(filter="app-io")'`
- **Test commands** use `devtools::test(filter = "<name>")`, which runs `tests/testthat/test-<name>.R`. `devtools::load_all()` makes internal (non-exported) functions directly callable in tests.
- **Rtools is not installed on this machine** — `R CMD check` / package install cannot run here. `load_all`/`test`/`document` work. Do NOT add tests that require an installed package (this rules out `shinytest2::AppDriver`, which launches a separate R process that `library(mixdra)`); the smoke test in Task 5 asserts the app object *assembles* instead.
- Existing suite is 16 files / 54 fast tests (plus the slow `validation-binary`), all green — every task here must keep them green.

## File structure

| File | Responsibility |
|---|---|
| `R/app-io.R` (create) | Pure helpers: `upload_schema`, `template_df`, `validate_upload`, `read_upload`, `to_engine_df`, `collect_bounds`. No Shiny, no plotly. |
| `R/app-intro.R` (create) | `intro_ui(id)` / `intro_server(id, meta)` — experiment metadata into shared `meta`. |
| `R/app-single.R` (create) | `single_ui(id)` / `single_server(id, meta)` — one-chemical upload → `analyse_single` → plots + table. |
| `R/app-binary.R` (create) | `binary_ui(id)` / `binary_server(id, meta)` — binary upload → `analyse_mixture` → tables + plots, Advanced bounds. |
| `R/app-run.R` (create) | `app_ui()` / `app_server` (internal) assemble the navbar + wire modules; `run_app()` (exported) guards deps and launches. |
| `inst/app/app.R` (create) | Two-line launcher referencing the package, for `shiny::runApp`. |
| `DESCRIPTION` (modify) | Add `shiny`, `bslib`, `DT`, `shinytest2` to Suggests. |
| `tests/testthat/test-app-io.R` (create) | Rigorous helper unit tests (no Shiny). |
| `tests/testthat/test-app-modules.R` (create) | `testServer` tests for the three module servers. |
| `tests/testthat/test-app-run.R` (create) | App-assembles smoke test (`skip_if_not_installed`). |

**Key engine facts the modules rely on (verified against current `R/`):**
- `analyse_single(df)` reads `df$C1` + (`df$Res` or `df$Affected`/`df$Exposed`); returns a `fit_single` object tagged `kind = "single"`.
- `analyse_mixture(df, reference, response, start, alpha, lower, upper, n_starts, time_limit)` — `response` is `"continuous"` or `"binary"`; there is **no `fixed` argument**, so parameter-fixing in the UI is emulated by passing equal `lower`/`upper` for that parameter.
- `result_table(res)`, `res$comparison`, `res$chosen`, `param_ci(fit, df, reference, deviation, response, level)`.
- Plot builders accept `(fit, df)`; for single fits `df` must carry a `C1` column (so the upload's `Conc` is renamed to `C1` before fitting/plotting).
- The UI's response labels are `"continuous"`/`"quantal"`; the engine's are `"continuous"`/`"binary"`. Convert with: `engine_response <- if (response == "quantal") "binary" else "continuous"`.

---

## Task 1: Pure I/O helpers

**Files:**
- Create: `R/app-io.R`
- Test: `tests/testthat/test-app-io.R`

- [x] **Step 1: Write the failing tests**

Create `tests/testthat/test-app-io.R`:

```r
test_that("upload_schema returns the fixed columns per stage and response", {
  expect_equal(upload_schema("single", "continuous"), c("Conc", "Res"))
  expect_equal(upload_schema("single", "quantal"), c("Conc", "Affected", "Exposed"))
  expect_equal(upload_schema("binary", "continuous"), c("C1", "C2", "Res"))
  expect_equal(upload_schema("binary", "quantal"), c("C1", "C2", "Affected", "Exposed"))
})

test_that("template_df has exactly the schema columns and at least one example row", {
  d <- template_df("binary", "continuous")
  expect_equal(names(d), c("C1", "C2", "Res"))
  expect_gt(nrow(d), 0)
  expect_true(all(vapply(d, is.numeric, logical(1))))
  # binary template includes single-chemical rows (one chem at 0) for seeding
  expect_true(any(d$C2 == 0 & d$C1 > 0))
  expect_true(any(d$C1 == 0 & d$C2 > 0))
})

test_that("validate_upload returns no errors for a valid file", {
  good <- data.frame(C1 = c(0, 1, 0, 2), C2 = c(0, 0, 1, 2), Res = c(100, 50, 60, 20))
  expect_length(validate_upload(good, "binary", "continuous"), 0)
})

test_that("validate_upload reports missing columns", {
  bad <- data.frame(C1 = c(0, 1), Res = c(100, 50))   # no C2
  errs <- validate_upload(bad, "binary", "continuous")
  expect_match(paste(errs, collapse = " "), "C2")
})

test_that("validate_upload rejects non-numeric, negative conc, and Affected > Exposed", {
  expect_match(paste(validate_upload(
    data.frame(C1 = c("a", "b"), C2 = c(0, 1), Res = c(1, 2)),
    "binary", "continuous"), collapse = " "), "[Nn]on-numeric")
  expect_match(paste(validate_upload(
    data.frame(C1 = c(-1, 1), C2 = c(0, 1), Res = c(1, 2)),
    "binary", "continuous"), collapse = " "), ">= 0|negative|0")
  expect_match(paste(validate_upload(
    data.frame(C1 = c(0, 1), C2 = c(0, 1), Affected = c(2, 12), Exposed = c(10, 10)),
    "binary", "quantal"), collapse = " "), "Affected")
})

test_that("validate_upload needs >= 4 distinct concentrations for a single fit", {
  short <- data.frame(Conc = c(0, 1, 2), Res = c(100, 50, 10))
  expect_match(paste(validate_upload(short, "single", "continuous"), collapse = " "),
               "distinct")
})

test_that("read_upload reads a CSV file into a data frame", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Conc = c(0, 1), Res = c(100, 50)), path, row.names = FALSE)
  d <- read_upload(path)
  expect_equal(names(d), c("Conc", "Res"))
  expect_equal(nrow(d), 2)
})

test_that("to_engine_df renames Conc to C1 for single only", {
  s <- to_engine_df(data.frame(Conc = c(0, 1), Res = c(1, 2)), "single")
  expect_true("C1" %in% names(s))
  expect_false("Conc" %in% names(s))
  b <- to_engine_df(data.frame(C1 = 0, C2 = 1, Res = 3), "binary")
  expect_equal(names(b), c("C1", "C2", "Res"))
})

test_that("collect_bounds keeps only supplied bounds, named by parameter", {
  vals <- list(lo_max = NA, hi_max = 1, lo_slope1 = NA, hi_slope1 = NA,
               lo_slope2 = NA, hi_slope2 = NA, lo_ec501 = 0.01, hi_ec501 = NA,
               lo_ec502 = NA, hi_ec502 = NA)
  b <- collect_bounds(vals)
  expect_equal(b$upper, c(max = 1))
  expect_equal(b$lower, c(ec501 = 0.01))
})

test_that("collect_bounds returns NULL bounds when nothing supplied", {
  vals <- setNames(as.list(rep(NA, 10)),
                   c(paste0("lo_", c("max","slope1","slope2","ec501","ec502")),
                     paste0("hi_", c("max","slope1","slope2","ec501","ec502"))))
  b <- collect_bounds(vals)
  expect_null(b$lower)
  expect_null(b$upper)
})
```

- [x] **Step 2: Run the tests to verify they fail**

Run: `devtools::test(filter = "app-io")`
Expected: FAIL — `could not find function "upload_schema"`.

- [x] **Step 3: Write the implementation**

Create `R/app-io.R`:

```r
# Pure helpers for the Shiny app's data layer. No Shiny, no plotly here, so these
# are fully unit-testable and run without the UI stack installed.

#' Fixed column schema for a stage and response type
#' @param stage "single" or "binary".
#' @param response "continuous" or "quantal".
#' @return Character vector of required column names.
#' @keywords internal
upload_schema <- function(stage, response) {
  stage <- match.arg(stage, c("single", "binary"))
  response <- match.arg(response, c("continuous", "quantal"))
  if (stage == "single") {
    if (response == "continuous") c("Conc", "Res") else c("Conc", "Affected", "Exposed")
  } else {
    if (response == "continuous") c("C1", "C2", "Res") else c("C1", "C2", "Affected", "Exposed")
  }
}

#' Example template data frame for a stage and response type
#'
#' Returns a small, illustrative dataset with exactly the schema columns. The
#' binary template includes single-chemical rows (one chemical at 0) because
#' [mixdra::analyse_mixture()] seeds itself from them.
#' @inheritParams upload_schema
#' @return A data frame the user can download, fill in, and re-upload.
#' @keywords internal
template_df <- function(stage, response) {
  cols <- upload_schema(stage, response)
  if (stage == "single") {
    conc <- c(0, 0.1, 0.3, 1, 3, 10)
    if (response == "continuous") {
      data.frame(Conc = conc, Res = c(100, 96, 82, 50, 18, 4))
    } else {
      data.frame(Conc = conc, Affected = c(0, 1, 2, 5, 8, 10), Exposed = rep(10, 6))
    }
  } else {
    # single-chemical series for each chemical + a few mixture rows
    c1 <- c(0, 0.1, 0.3, 1, 0, 0, 0, 0.1, 0.3, 1)
    c2 <- c(0, 0,   0,   0, 0.1, 0.3, 1, 0.1, 0.3, 1)
    if (response == "continuous") {
      data.frame(C1 = c1, C2 = c2,
                 Res = c(100, 80, 55, 20, 88, 70, 35, 72, 45, 12))
    } else {
      data.frame(C1 = c1, C2 = c2,
                 Affected = c(0, 2, 4, 8, 1, 3, 7, 3, 6, 9), Exposed = rep(10, 10))
    }
  }
}

#' Validate an uploaded data frame against the fixed schema
#'
#' @inheritParams upload_schema
#' @param df The uploaded data frame.
#' @return Character vector of human-readable error messages; empty if valid.
#' @keywords internal
validate_upload <- function(df, stage, response) {
  errs <- character(0)
  req_cols <- upload_schema(stage, response)

  missing <- setdiff(req_cols, names(df))
  if (length(missing))
    errs <- c(errs, paste0("Missing column(s): ", paste(missing, collapse = ", "),
                           ". Expected exactly: ", paste(req_cols, collapse = ", "), "."))

  present <- intersect(req_cols, names(df))
  non_num <- present[!vapply(df[present], is.numeric, logical(1))]
  if (length(non_num))
    errs <- c(errs, paste0("Non-numeric column(s): ", paste(non_num, collapse = ", "), "."))

  # Range checks only on numeric columns that are present.
  conc_cols <- intersect(c("Conc", "C1", "C2"), present)
  conc_ok <- conc_cols[vapply(df[conc_cols], is.numeric, logical(1))]
  if (length(conc_ok) && any(unlist(df[conc_ok]) < 0, na.rm = TRUE))
    errs <- c(errs, "Concentrations must be >= 0.")

  if (response == "quantal" && all(c("Affected", "Exposed") %in% present) &&
      is.numeric(df$Affected) && is.numeric(df$Exposed)) {
    if (any(df$Affected < 0 | df$Exposed < 0, na.rm = TRUE))
      errs <- c(errs, "Affected and Exposed must be non-negative.")
    if (any(df$Affected > df$Exposed, na.rm = TRUE))
      errs <- c(errs, "Affected must be <= Exposed.")
  }

  if (stage == "single" && "Conc" %in% present && is.numeric(df$Conc)) {
    n_distinct <- length(unique(df$Conc[!is.na(df$Conc)]))
    if (n_distinct < 4)
      errs <- c(errs, "Need at least 4 distinct concentrations to fit a single-chemical curve.")
  }

  errs
}

#' Read an uploaded CSV file
#' @param path File path.
#' @return A data frame.
#' @keywords internal
read_upload <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
}

#' Map an uploaded data frame to the engine's column convention
#'
#' The single-chemical template uses `Conc`; the engine expects `C1`.
#' @param df Uploaded data frame.
#' @param stage "single" or "binary".
#' @return The data frame with engine-ready column names.
#' @keywords internal
to_engine_df <- function(df, stage) {
  if (stage == "single") names(df)[names(df) == "Conc"] <- "C1"
  df
}

#' Assemble lower/upper bound vectors from Advanced-panel inputs
#'
#' Reads `lo_<param>` / `hi_<param>` values for the five binary base parameters;
#' blank/NA entries are dropped. To FIX a parameter, set its lower and upper to
#' the same value (the engine has no `fixed` argument via `analyse_mixture`).
#' @param values Named list (e.g. a Shiny `input`) holding `lo_*`/`hi_*` numbers.
#' @return A list with `lower` and `upper` named numeric vectors (or NULL).
#' @keywords internal
collect_bounds <- function(values) {
  params <- c("max", "slope1", "slope2", "ec501", "ec502")
  pick <- function(prefix) {
    v <- vapply(params, function(p) {
      x <- values[[paste0(prefix, p)]]
      if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else as.numeric(x)
    }, numeric(1))
    v <- v[!is.na(v)]
    if (length(v)) v else NULL
  }
  list(lower = pick("lo_"), upper = pick("hi_"))
}
```

- [x] **Step 4: Run the tests to verify they pass**

Run: `devtools::test(filter = "app-io")`
Expected: PASS (10 tests).

- [x] **Step 5: Commit**

```bash
git add R/app-io.R tests/testthat/test-app-io.R
git commit -m "feat: pure I/O helpers for the Shiny app (schema, template, validate)"
```

---

## Task 2: Introduction module

**Files:**
- Create: `R/app-intro.R`
- Test: `tests/testthat/test-app-modules.R`

- [x] **Step 1: Write the failing test**

Create `tests/testthat/test-app-modules.R`:

```r
skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

test_that("intro_server writes form inputs into the shared meta store", {
  meta <- shiny::reactiveValues()
  shiny::testServer(intro_server, args = list(meta = meta), {
    session$setInputs(expID = "E1", species = "Daphnia", endpoint = "reproduction",
                      chem1 = "CPF", chem2 = "MPs", unit = "mg/L")
    expect_equal(meta$expID, "E1")
    expect_equal(meta$chem1, "CPF")
    expect_equal(meta$chem2, "MPs")
    expect_equal(meta$unit, "mg/L")
    expect_equal(meta$endpoint, "reproduction")
  })
})

test_that("intro_ui builds a Shiny UI fragment", {
  ui <- intro_ui("intro")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})
```

- [x] **Step 2: Run the test to verify it fails**

Run: `devtools::test(filter = "app-modules")`
Expected: FAIL — `could not find function "intro_server"`.

- [x] **Step 3: Write the implementation**

Create `R/app-intro.R`:

```r
# Introduction stage: collects experiment metadata into the shared `meta` store.
# `meta` is a reactiveValues created by app_server() and shared with every stage,
# which read it only for plot titles and axis labels.

#' Introduction stage UI
#' @param id Module id.
#' @keywords internal
intro_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::textInput(ns("expID"), "Experiment ID"),
      shiny::textInput(ns("species"), "Species"),
      shiny::textInput(ns("endpoint"), "Endpoint"),
      shiny::helpText("Enter chemical names in the order used throughout the app."),
      shiny::textInput(ns("chem1"), "Chemical 1 name"),
      shiny::textInput(ns("chem2"), "Chemical 2 name"),
      shiny::textInput(ns("unit"), "Concentration unit (e.g. mg/L)")
    ),
    bslib::card(
      bslib::card_header("Welcome"),
      shiny::p("Upload concentration-response data on the Single Chemical and ",
               "Binary Mixture tabs. The experiment details entered here are used ",
               "to label plots and axes.")
    )
  )
}

#' Introduction stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
intro_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      meta$expID    <- input$expID
      meta$species  <- input$species
      meta$endpoint <- input$endpoint
      meta$chem1    <- input$chem1
      meta$chem2    <- input$chem2
      meta$unit     <- input$unit
    })
  })
}
```

- [x] **Step 4: Run the test to verify it passes**

Run: `devtools::test(filter = "app-modules")`
Expected: PASS (2 tests).

- [x] **Step 5: Commit**

```bash
git add R/app-intro.R tests/testthat/test-app-modules.R
git commit -m "feat: Introduction stage module (experiment metadata)"
```

---

## Task 3: Single Chemical module

**Files:**
- Create: `R/app-single.R`
- Test: `tests/testthat/test-app-modules.R` (append)

- [x] **Step 1: Write the failing test**

Append to `tests/testthat/test-app-modules.R`:

```r
test_that("single_server validates, fits, and exposes a single fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    # Valid continuous single-chemical CSV with >= 4 distinct concentrations.
    path <- tempfile(fileext = ".csv")
    utils::write.csv(
      data.frame(Conc = c(0, 0.1, 0.3, 1, 3, 10),
                 Res  = ll3_predict(c(0, 0.1, 0.3, 1, 3, 10), 100, 2, 0.5)),
      path, row.names = FALSE)
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "single.csv"))
    expect_length(errs(), 0)               # passes validation
    session$setInputs(fit = 1)
    expect_equal(fit_r()$kind, "single")
    expect_equal(unname(round(fit_r()$par[["ec50"]], 2)), 0.5)
  })
})

test_that("single_server surfaces validation errors and withholds a fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(Conc = c(0, 1, 2), Res = c(100, 50, 10)),
                     path, row.names = FALSE)   # only 3 distinct concentrations
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "short.csv"))
    expect_gt(length(errs()), 0)
    expect_match(paste(errs(), collapse = " "), "distinct")
  })
})

test_that("single_ui builds a Shiny UI fragment", {
  expect_true(inherits(single_ui("single"), c("shiny.tag", "shiny.tag.list", "bslib_fragment")))
})
```

- [x] **Step 2: Run the test to verify it fails**

Run: `devtools::test(filter = "app-modules")`
Expected: FAIL — `could not find function "single_server"`.

- [x] **Step 3: Write the implementation**

Create `R/app-single.R`:

```r
# Single Chemical stage: upload one chemical's dose-response data, fit a
# three-parameter log-logistic curve via analyse_single(), and show the curve,
# observed-vs-predicted, and a parameter table.

#' Concentration-axis label from the shared meta store
#'
#' With `chem_field` (e.g. "chem1") uses that chemical's name; without it uses a
#' generic "Concentration". Appends the unit when present. Used by the single
#' stage (generic) and the binary stage (per chemical).
#' @keywords internal
axis_label <- function(meta, chem_field = NULL) {
  nm <- if (!is.null(chem_field)) meta[[chem_field]] else NULL
  unit <- meta$unit
  base <- if (!is.null(nm) && nzchar(nm)) nm else "Concentration"
  if (!is.null(unit) && nzchar(unit)) paste0(base, " (", unit, ")") else base
}

#' Single Chemical stage UI
#' @param id Module id.
#' @keywords internal
single_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::uiOutput(ns("errors")),
      shiny::actionButton(ns("fit"), "Fit", class = "btn-primary")
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Dose-response curve"),
                  plotly::plotlyOutput(ns("dr"))),
      bslib::card(bslib::card_header("Observed vs predicted"),
                  plotly::plotlyOutput(ns("op")))
    ),
    bslib::card(bslib::card_header("Parameters"), DT::DTOutput(ns("params")))
  )
}

#' Single Chemical stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
single_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() paste0("single_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("single", input$response), file, row.names = FALSE)
    )

    parsed <- shiny::reactive({
      shiny::req(input$file)
      read_upload(input$file$datapath)
    })

    errs <- shiny::reactive({
      shiny::req(input$file)
      validate_upload(parsed(), "single", input$response)
    })

    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    fit_df <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      to_engine_df(parsed(), "single")
    })

    fit_r <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      tryCatch(analyse_single(fit_df()), error = function(e) {
        shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
        NULL
      })
    })

    output$dr <- plotly::renderPlotly({
      shiny::req(fit_r())
      p <- plot_dose_response(fit_r(), fit_df())
      plotly::layout(p, xaxis = list(title = axis_label(meta)),
                     yaxis = list(title = if (!is.null(meta$endpoint) && nzchar(meta$endpoint))
                                            meta$endpoint else "Response"))
    })
    output$op <- plotly::renderPlotly({
      shiny::req(fit_r())
      plot_obs_pred(fit_r(), fit_df())
    })
    output$params <- DT::renderDT({
      shiny::req(fit_r())
      p <- fit_r()$par
      DT::datatable(
        data.frame(Parameter = c("max", "slope", "ec50", "SSR", "n"),
                   Value = c(round(unname(p[c("max", "slope", "ec50")]), 4),
                             round(fit_r()$ssr, 2), nrow(fit_df()))),
        rownames = FALSE, options = list(dom = "t"))
    })
  })
}
```

- [x] **Step 4: Run the test to verify it passes**

Run: `devtools::test(filter = "app-modules")`
Expected: PASS (5 tests so far in this file).

- [x] **Step 5: Commit**

```bash
git add R/app-single.R tests/testthat/test-app-modules.R
git commit -m "feat: Single Chemical stage module (upload, fit, plots)"
```

---

## Task 4: Binary Mixture module

**Files:**
- Create: `R/app-binary.R`
- Test: `tests/testthat/test-app-modules.R` (append)

- [x] **Step 1: Write the failing test**

Append to `tests/testthat/test-app-modules.R`:

```r
test_that("binary_server fits the four models and exposes the chosen model", {
  skip_on_cran()
  meta <- shiny::reactiveValues()
  shiny::testServer(binary_server, args = list(meta = meta), {
    # Use the validated continuous binary fixture (C1, C2, Res).
    csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
    skip_if_not(file.exists(csv), "binary fixture missing")
    session$setInputs(response = "continuous", reference = "CA",
                      thorough = FALSE,
                      n_starts = 1, alpha = 0.05, time_limit = 30,
                      lo_max = NA, hi_max = NA, lo_slope1 = NA, hi_slope1 = NA,
                      lo_slope2 = NA, hi_slope2 = NA, lo_ec501 = NA, hi_ec501 = NA,
                      lo_ec502 = NA, hi_ec502 = NA,
                      file = list(datapath = csv, name = "binary.csv"))
    expect_length(errs(), 0)
    session$setInputs(fit = 1)
    res <- res_r()
    expect_setequal(names(res$fits), c("reference", "SA", "DR", "DL"))
    expect_true(res$chosen %in% c("reference", "SA", "DR", "DL"))
    # the displayed fit defaults to the chosen model (deviation tag matches)
    expect_equal(shown_fit()$deviation, res$chosen)
  })
})

test_that("binary_ui builds a Shiny UI fragment", {
  expect_true(inherits(binary_ui("binary"), c("shiny.tag", "shiny.tag.list", "bslib_fragment")))
})
```

- [x] **Step 2: Run the test to verify it fails**

Run: `devtools::test(filter = "app-modules")`
Expected: FAIL — `could not find function "binary_server"`.

- [x] **Step 3: Write the implementation**

Create `R/app-binary.R`:

```r
# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows), fit CA/IA reference + SA/DR/DL deviations via analyse_mixture(),
# compare them, and visualise the chosen (or any) fit.

#' One lower/upper bound numeric-input pair for a parameter
#' @keywords internal
bound_row <- function(ns, param, label, hi_default = NA) {
  shiny::fluidRow(
    shiny::column(4, shiny::tags$small(label)),
    shiny::column(4, shiny::numericInput(ns(paste0("lo_", param)), NULL, value = NA)),
    shiny::column(4, shiny::numericInput(ns(paste0("hi_", param)), NULL, value = hi_default))
  )
}

#' Binary Mixture stage UI
#' @param id Module id.
#' @keywords internal
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
      shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
      shiny::uiOutput(ns("thorough_note")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced",
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1),
          shiny::helpText("Parameter bounds (blank = default). Set lower = upper to fix."),
          shiny::tags$div(shiny::tags$small(shiny::tags$b("param / lower / upper"))),
          bound_row(ns, "max", "max", hi_default = 1),
          bound_row(ns, "slope1", "slope1"),
          bound_row(ns, "slope2", "slope2"),
          bound_row(ns, "ec501", "ec501"),
          bound_row(ns, "ec502", "ec502")
        )
      ),
      shiny::uiOutput(ns("errors")),
      shiny::actionButton(ns("fit"), "Fit", class = "btn-primary")
    ),
    shiny::selectInput(ns("model"), "Model to display", choices = NULL),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Dose-response (Chemical 1)"),
                  plotly::plotlyOutput(ns("dr1"))),
      bslib::card(bslib::card_header("Dose-response (Chemical 2)"),
                  plotly::plotlyOutput(ns("dr2")))
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("3-D response surface"),
                  plotly::plotlyOutput(ns("surface"))),
      bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                  plotly::plotlyOutput(ns("isobole")))
    ),
    bslib::card(bslib::card_header("Observed vs predicted"),
                plotly::plotlyOutput(ns("op"))),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Results table"), DT::DTOutput(ns("results"))),
      bslib::card(bslib::card_header("Model comparison"), DT::DTOutput(ns("comparison")))
    ),
    bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                DT::DTOutput(ns("cis")))
  )
}

#' Binary Mixture stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
binary_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$template <- shiny::downloadHandler(
      filename = function() paste0("binary_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("binary", input$response), file, row.names = FALSE)
    )

    output$thorough_note <- shiny::renderUI({
      if (isTRUE(input$thorough))
        shiny::div(class = "text-warning",
                   shiny::tags$small("Multi-start fitting may take several minutes."))
    })

    parsed <- shiny::reactive({
      shiny::req(input$file); read_upload(input$file$datapath)
    })
    errs <- shiny::reactive({
      shiny::req(input$file); validate_upload(parsed(), "binary", input$response)
    })
    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    fit_df <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      to_engine_df(parsed(), "binary")
    })

    res_r <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      df <- to_engine_df(parsed(), "binary")
      engine_response <- if (input$response == "quantal") "binary" else "continuous"
      n_starts <- if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts
      b <- collect_bounds(shiny::reactiveValuesToList(input))
      shiny::withProgress(message = "Fitting models...", value = 0.5, {
        tryCatch(
          analyse_mixture(df, reference = input$reference, response = engine_response,
                          alpha = input$alpha, lower = b$lower, upper = b$upper,
                          n_starts = n_starts, time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
            NULL
          })
      })
    })

    # Refresh the (static) model picker after each fit, defaulting to the chosen model.
    shiny::observeEvent(res_r(), {
      shiny::updateSelectInput(session, "model",
                               choices = names(res_r()$fits), selected = res_r()$chosen)
    })

    # Fall back to the chosen model until the picker's input has populated (also
    # makes the reactive testable under shiny::testServer, where updateSelectInput
    # does not round-trip an input value).
    shown_fit <- shiny::reactive({
      shiny::req(res_r())
      m <- input$model
      if (is.null(m) || !m %in% names(res_r()$fits)) m <- res_r()$chosen
      res_r()$fits[[m]]
    })

    output$dr1 <- plotly::renderPlotly({
      plotly::layout(plot_dose_response(shown_fit(), fit_df(), chem = 1),
                     xaxis = list(title = axis_label(meta, "chem1")))
    })
    output$dr2 <- plotly::renderPlotly({
      plotly::layout(plot_dose_response(shown_fit(), fit_df(), chem = 2),
                     xaxis = list(title = axis_label(meta, "chem2")))
    })
    output$surface <- plotly::renderPlotly({
      plot_surface(shown_fit(), fit_df())
    })
    output$isobole <- plotly::renderPlotly({
      plot_isobole(shown_fit(), fit_df(), reference_fit = res_r()$fits$reference)
    })
    output$op <- plotly::renderPlotly({
      plot_obs_pred(shown_fit(), fit_df())
    })

    output$results <- DT::renderDT({
      shiny::req(res_r())
      tab <- round(result_table(res_r()), 4)
      DT::datatable(as.data.frame(tab), options = list(dom = "t"))
    })
    output$comparison <- DT::renderDT({
      shiny::req(res_r())
      DT::datatable(res_r()$comparison, rownames = FALSE, options = list(dom = "t"))
    })
    output$cis <- DT::renderDT({
      shiny::req(shown_fit())
      f <- shown_fit()
      DT::datatable(param_ci(f, fit_df(), f$reference, f$deviation, f$response),
                    rownames = FALSE, options = list(dom = "t"))
    })

    res_r   # return for testability
  })
}
```

- [x] **Step 4: Run the test to verify it passes**

Run: `devtools::test(filter = "app-modules")`
Expected: PASS (7 tests in this file). The binary fit runs at `n_starts = 1` so it is fast.

- [x] **Step 5: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat: Binary Mixture stage module (fit, compare, plots, advanced bounds)"
```

---

## Task 5: App assembly, run_app(), dependencies, and docs

**Files:**
- Create: `R/app-run.R`
- Create: `inst/app/app.R`
- Modify: `DESCRIPTION`
- Test: `tests/testthat/test-app-run.R`
- Generated: `NAMESPACE`, `man/*.Rd`

- [x] **Step 1: Add UI dependencies to DESCRIPTION**

In `DESCRIPTION`, change the `Suggests:` block from:

```
Suggests:
    plotly,
    readxl,
    testthat (>= 3.0.0)
```

to:

```
Suggests:
    bslib,
    DT,
    plotly,
    readxl,
    shiny,
    shinytest2,
    testthat (>= 3.0.0)
```

- [x] **Step 2: Write the failing test**

Create `tests/testthat/test-app-run.R`:

```r
skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

test_that("app_ui assembles a bslib page", {
  ui <- app_ui()
  expect_true(inherits(ui, c("shiny.tag", "shiny.tag.list", "bslib_page", "bslib_fragment")))
})

test_that("app_server is a function of (input, output, session)", {
  expect_true(is.function(app_server))
  expect_setequal(names(formals(app_server)), c("input", "output", "session"))
})

test_that("run_app is exported and errors clearly when a dependency is missing", {
  # We cannot uninstall packages in a test; just assert the guard logic exists by
  # checking run_app is a function that references the required packages.
  expect_true(is.function(run_app))
  body_txt <- paste(deparse(body(run_app)), collapse = " ")
  expect_match(body_txt, "bslib")
  expect_match(body_txt, "shiny")
})
```

- [x] **Step 3: Run the test to verify it fails**

Run: `devtools::test(filter = "app-run")`
Expected: FAIL — `could not find function "app_ui"`.

- [x] **Step 4: Write the implementation**

Create `R/app-run.R`:

```r
# App assembly + launcher. app_ui()/app_server are internal; run_app() is the
# single exported entry point and guards the Suggested UI dependencies.

#' Assemble the navbar UI from the three stage modules
#' @keywords internal
app_ui <- function() {
  bslib::page_navbar(
    title = "Mixture Toxicity (mixdra)",
    bslib::nav_panel("Introduction", intro_ui("intro")),
    bslib::nav_panel("Single Chemical", single_ui("single")),
    bslib::nav_panel("Binary Mixture", binary_ui("binary"))
  )
}

#' Wire the three stage module servers around a shared meta store
#' @keywords internal
app_server <- function(input, output, session) {
  meta <- shiny::reactiveValues()
  intro_server("intro", meta)
  single_server("single", meta)
  binary_server("binary", meta)
}

#' Launch the mixdra Shiny app
#'
#' Starts the interactive app (Introduction, Single Chemical, Binary Mixture).
#' The UI stack (`shiny`, `bslib`, `plotly`, `DT`) is a set of Suggested
#' dependencies; this function stops with an install hint if any are missing.
#' @param ... Passed to [shiny::runApp()] (e.g. `launch.browser`, `port`).
#' @return Invisibly; called for the side effect of running the app.
#' @export
run_app <- function(...) {
  needed <- c("shiny", "bslib", "plotly", "DT")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("run_app() needs these packages: ", paste(missing, collapse = ", "),
         ". Install with install.packages(c(",
         paste0("\"", missing, "\"", collapse = ", "), ")).", call. = FALSE)
  shiny::runApp(shiny::shinyApp(ui = app_ui(), server = app_server), ...)
}
```

Create `inst/app/app.R` (for `shiny::runApp(system.file("app", ...))` users and deployment):

```r
# Launcher for the mixdra Shiny app. Assumes the mixdra package is installed.
library(mixdra)
shiny::shinyApp(ui = mixdra:::app_ui(), server = mixdra:::app_server)
```

- [x] **Step 5: Run the test to verify it passes**

Run: `devtools::test(filter = "app-run")`
Expected: PASS (3 tests).

- [x] **Step 6: Document and run the full fast suite**

Run: `devtools::document()` (regenerates `NAMESPACE` with `export(run_app)` and `man/run_app.Rd`)
Run: `devtools::test()`
Expected: the existing 54 fast tests plus the new app-io (10), app-modules (7), and app-run (3) tests pass, 0 failures. (The slow `validation-binary` file is unaffected.)

- [x] **Step 7: Commit**

```bash
git add R/app-run.R inst/app/app.R DESCRIPTION NAMESPACE man/
git commit -m "feat: assemble mixdra Shiny app and export run_app()"
```

---

## Manual verification (after the plan, optional)

Build a full-browser check is not automated here (needs an installed package + a
browser, which `R CMD check` cannot run on this machine). To verify manually in
RStudio after `devtools::load_all()`:

```r
mixdra::run_app()    # or: devtools::load_all(); run_app()
```

Then: enter chemical names on Introduction; on Single Chemical download the template,
fill/upload it, Fit, and confirm the curve + parameters; on Binary Mixture upload the
continuous fixture (`tests/testthat/fixtures/binary_mps_cpf_continuous.csv`), Fit at
default starts, and confirm the results table, comparison, surface, and isobole render
and that the model picker defaults to the chosen model.

---

## Self-review notes (for the implementer)

- **Response label mismatch is deliberate:** the UI uses `"continuous"`/`"quantal"`;
  the engine uses `"continuous"`/`"binary"`. `binary_server` converts before calling
  `analyse_mixture`. `analyse_single` infers from columns, so it needs no conversion.
- **Single fits need `C1`:** `to_engine_df(df, "single")` renames `Conc` → `C1` so the
  plot builders (which read `df$C1`) work.
- **Fixing a parameter** is done by setting equal lower/upper in the Advanced panel —
  `analyse_mixture` exposes `lower`/`upper` but not `fixed`.
- **testServer reactives** (`errs`, `fit_r`, `res_r`, `fit_df`) are referenced by name
  in tests; keep those names if you refactor.
- **No browser tests in the suite:** the smoke test asserts assembly only, because
  `shinytest2::AppDriver` requires an installed package and cannot run under
  `load_all()` on this machine.
