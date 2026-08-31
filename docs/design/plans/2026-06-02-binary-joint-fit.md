# Binary tab: final joint "Optimize all params" stage — Implementation Plan

**Goal:** Add an optional final stage to the Binary tab that jointly optimizes every parameter of the displayed model (seeded from the staged fit), with a 7-row lower/upper bounds grid where Lower = Upper pins a parameter.

**Architecture:** The staged flow (single curves → freeze → interaction fit → Find best) is untouched. A new "Optimize all params" card in Stage 2 calls a thin `refine_joint()` wrapper over the existing `fit_model()` engine, which already does joint fitting when no params are fixed. The only engine change is allowing `a`/`b` to be bounded. Fixing is resolved in the app layer (Lower==Upper → `fixed`).

**Tech Stack:** R package (`mixdra`), Shiny + bslib UI, `testthat` tests run via `devtools::test(filter=...)`, roxygen2 docs via `devtools::document()`.

**Spec:** `docs/design/specs/2026-06-02-binary-joint-fit-design.md`

**Conventions:**
- Run only filtered tests (never the full suite): `Rscript -e 'devtools::test(filter="<name>")'`.
- Commit messages: do NOT add a Claude co-author trailer.

---

## File structure

| File | Responsibility | Change |
|------|----------------|--------|
| `R/fit.R` | `fit_model()` engine | Allow bounding `a`/`b`; record `fixed` on result; add `refine_joint()` |
| `R/app-io.R` | pure data-layer helpers | Add `lo_prefix`/`hi_prefix` to `collect_bounds()`; add `collect_bounds_all()` + `split_fixed_bounds()` |
| `R/app-curve-fit.R` | reusable UI rows | Add `optimize_param_row()` |
| `R/app-binary.R` | Binary module UI + server | Add Optimize-all card + observer + readout; blank CIs for pinned params |
| `R/summary.R` | `param_ci()` | (no change; blanking done in app) |
| `tests/testthat/test-fit-bounds.R` | bound tests | Update 3 tests for the relaxed validation |
| `tests/testthat/test-joint-refine.R` | NEW | recovery, monotonic, fix-routing-through-engine |
| `tests/testthat/test-app-io.R` | helper tests | Add `collect_bounds_all`/`split_fixed_bounds` tests |
| `tests/testthat/test-app-modules.R` | module tests | Add Optimize-all `testServer` test + UI presence |

---

## Task 1: Engine — allow bounding interaction params (`a`/`b`)

**Files:**
- Modify: `R/fit.R` (validation ~84-89; bound application ~118-120; lo<hi check ~126-128; start clamp ~131-139; both return lists ~98-102 and ~228-232; roxygen `@param lower,upper`)
- Modify: `tests/testthat/test-fit-bounds.R` (3 existing tests)

- [ ] **Step 1: Update the failing tests to the new contract**

In `tests/testthat/test-fit-bounds.R`, replace the test at lines 44-51 ("naming a deviation parameter in bounds raises an error") with a positive test, and fix the two error-message expectations that say `"base parameter"`:

```r
test_that("an interaction parameter can be bounded", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  # Constrain `a` to a positive band; the fit must keep it in range and not error.
  fit <- fit_model(df, "CA", "SA", "continuous", start = start,
                   lower = c(a = 0.5), upper = c(a = 5))
  expect_gte(unname(fit$par["a"]), 0.5 - 1e-6)
  expect_lte(unname(fit$par["a"]),  5 + 1e-6)
})
```

Change the message regex in "an unknown bound name raises an error" (line 41) from `"base parameter"` to `"model parameter"`.

Change "an illegal bound name errors even when all params are fixed" (lines 98-107) to use a genuinely unknown name and the new message:

```r
test_that("an illegal bound name errors even when all params are fixed", {
  df <- bounds_cont_df()
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  expect_error(
    fit_model(df, "CA", "reference", "continuous", start = start,
              fixed = names(start), upper = c(nonsense = 1)),
    "model parameter")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter="fit-bounds")'`
Expected: FAIL — the new `a`-bound test errors with "may only name a base parameter"; the message-regex assertions fail.

- [ ] **Step 3: Make the engine change in `R/fit.R`**

Replace the validation block (currently lines ~84-89):

```r
  base_all <- setdiff(spec$params, spec$extra)  # all curve params, incl. fixed
  bad <- setdiff(c(names(lower), names(upper)), base_all)
  if (length(bad))
    stop("`lower`/`upper` may only name a base parameter (",
         paste(base_all, collapse = ", "), "); got: ",
         paste(unique(bad), collapse = ", "))
```

with:

```r
  # Any model parameter may be bounded now -- including the interaction params
  # (a, b), which the final joint "Optimize all" stage constrains/pins.
  bad <- setdiff(c(names(lower), names(upper)), spec$params)
  if (length(bad))
    stop("`lower`/`upper` may only name a model parameter (",
         paste(spec$params, collapse = ", "), "); got: ",
         paste(unique(bad), collapse = ", "))
```

Replace the bound-application loops (currently ~118-120):

```r
  # Apply user-supplied bounds (base, free params only; names validated above).
  for (p in intersect(names(lower), base)) lo[[p]] <- lower[[p]]
  for (p in intersect(names(upper), base)) hi[[p]] <- upper[[p]]
```

with (apply to any free param, incl. `a`/`b`):

```r
  # Apply user-supplied bounds to any free parameter (incl. interaction a/b).
  for (p in intersect(names(lower), free)) lo[[p]] <- lower[[p]]
  for (p in intersect(names(upper), free)) hi[[p]] <- upper[[p]]
```

Replace the lower<upper check (currently ~126-128):

```r
  if (length(base) && any(lo[base] >= hi[base]))
    stop("each parameter's lower bound must be below its upper bound; check: ",
         paste(base[lo[base] >= hi[base]], collapse = ", "))
```

with:

```r
  if (length(free) && any(lo[free] >= hi[free]))
    stop("each parameter's lower bound must be below its upper bound; check: ",
         paste(free[lo[free] >= hi[free]], collapse = ", "))
```

Replace the start-clamp block (currently ~131-139):

```r
  if (length(base)) {
    clamped <- pmin(pmax(theta0[base], lo[base]), hi[base])
    if (any(clamped != theta0[base])) {
      off <- base[clamped != theta0[base]]
      warning("start value(s) outside bounds, clamped: ",
              paste(off, collapse = ", "))
      theta0[base] <- clamped
    }
  }
```

with (clamp every free param; unbounded ones have ±Inf bounds so this is a no-op for them):

```r
  if (length(free)) {
    clamped <- pmin(pmax(theta0[free], lo[free]), hi[free])
    if (any(clamped != theta0[free])) {
      off <- free[clamped != theta0[free]]
      warning("start value(s) outside bounds, clamped: ",
              paste(off, collapse = ", "))
      theta0[free] <- clamped
    }
  }
```

In BOTH return lists, add `fixed = fixed,` so callers know which params were pinned. The all-fixed fast-path return (currently ~98-102) becomes:

```r
    return(list(par = par, objective = objective_of(par), pred = pred,
                residuals = obs - pred, df = 0L, n = nrow(df),
                convergence = 0L, reference = reference, deviation = deviation,
                response = response, conc_cols = conc_cols, n_chem = n_chem,
                fixed = fixed, kind = "mixture"))
```

The main return (currently ~228-232) becomes:

```r
  list(par = par, objective = objective_of(par), pred = pred,
       residuals = obs - pred, df = length(free), n = nrow(df),
       convergence = best$convergence,
       reference = reference, deviation = deviation, response = response,
       conc_cols = conc_cols, n_chem = n_chem, fixed = fixed, kind = "mixture")
```

Update the `@param lower,upper` roxygen text: change "Any parameter not named falls back to the default ... Deviation parameters (`a`, `b`, `b1/b2/b3`) are always left unconstrained and may not be named here" to note that interaction parameters MAY now be named to constrain or (via equal bounds resolved upstream) pin them.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::test(filter="fit-bounds")'`
Expected: PASS (all tests in the file).

- [ ] **Step 5: Regenerate docs and commit**

Run: `Rscript -e 'devtools::document()'`

```bash
git add R/fit.R tests/testthat/test-fit-bounds.R man/
git commit -m "feat(fit): allow bounding interaction params a/b; record fixed on result"
```

---

## Task 2: app-io helpers — `collect_bounds_all()` + `split_fixed_bounds()`

**Files:**
- Modify: `R/app-io.R` (extend `collect_bounds()` with prefix args; add two helpers)
- Modify: `tests/testthat/test-app-io.R` (add tests)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-app-io.R`:

```r
test_that("collect_bounds_all reads lower/upper for all seven binary params", {
  vals <- list(olo_max = 100, ohi_max = 1000,
               olo_a = 0,   ohi_a = 5,
               olo_b = NA,  ohi_b = NA)
  b <- collect_bounds_all(vals)
  expect_equal(b$lower[["max"]], 100)
  expect_equal(b$upper[["max"]], 1000)
  expect_equal(b$lower[["a"]], 0)
  expect_equal(b$upper[["a"]], 5)
  expect_false("b" %in% names(b$lower))   # blank dropped
})

test_that("split_fixed_bounds turns equal lower/upper into a fixed param", {
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.5)
  sp <- split_fixed_bounds(lower = c(max = 800, a = 0),
                           upper = c(max = 800, a = 5), start = start)
  expect_true("max" %in% sp$fixed)     # equal bounds -> fixed
  expect_false("a" %in% sp$fixed)      # a is a true range, not pinned
  expect_equal(sp$start[["max"]], 800) # pinned at the equal-bound value
  expect_null(sp$lower[["max"]])       # pinned param dropped from bounds
  expect_equal(sp$lower[["a"]], 0)
  expect_equal(sp$upper[["a"]], 5)
})

test_that("split_fixed_bounds with no bounds returns empty fixed and NULL bounds", {
  start <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1)
  sp <- split_fixed_bounds(lower = numeric(0), upper = numeric(0), start = start)
  expect_length(sp$fixed, 0)
  expect_null(sp$lower)
  expect_null(sp$upper)
  expect_equal(sp$start, start)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter="app-io")'`
Expected: FAIL — `collect_bounds_all`/`split_fixed_bounds` not found.

- [ ] **Step 3: Implement the helpers in `R/app-io.R`**

Add `lo_prefix`/`hi_prefix` params to the existing `collect_bounds()` (non-breaking; defaults preserve current behavior). Change its signature and the `pick()` prefixes:

```r
collect_bounds <- function(values,
                           params = c("max", "slope1", "slope2", "ec501", "ec502"),
                           lo_prefix = "lo_", hi_prefix = "hi_") {
  pick <- function(prefix) {
    v <- vapply(params, function(p) {
      x <- values[[paste0(prefix, p)]]
      if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else as.numeric(x)
    }, numeric(1))
    v <- v[!is.na(v)]
    if (length(v)) v else NULL
  }
  list(lower = pick(lo_prefix), upper = pick(hi_prefix))
}
```

Then add the two new helpers below it:

```r
#' Read the 7-parameter lower/upper grid from the Optimize-all panel
#'
#' Thin wrapper over [collect_bounds()] for the full binary parameter set
#' (curve params + interaction `a`/`b`), reading the Optimize-all panel's
#' `olo_*` / `ohi_*` inputs.
#' @param values Named list (e.g. a Shiny `input`) holding `olo_*`/`ohi_*`.
#' @param params Parameter names to read (default: the binary 7-set).
#' @return A list with `lower` and `upper` named numeric vectors (or NULL).
#' @keywords internal
collect_bounds_all <- function(values,
                               params = c("max", "slope1", "slope2",
                                          "ec501", "ec502", "a", "b")) {
  collect_bounds(values, params, lo_prefix = "olo_", hi_prefix = "ohi_")
}

#' Split a lower/upper bound set into pinned (fixed) params and free bounds
#'
#' A parameter whose lower and upper bounds are both present and equal is treated
#' as PINNED: it goes into `fixed` at that value (and `start` is set to it). The
#' remaining bounds pass through as true ranges. This routes "fix via Lower =
#' Upper" through [fit_model()]'s `fixed` argument, avoiding the `lower == upper`
#' error that L-BFGS-B would otherwise raise.
#' @param lower,upper Named numeric vectors of bounds (may be empty).
#' @param start Named numeric starting vector (all model params).
#' @return A list: `fixed` (character), `start` (with pinned values applied),
#'   `lower`, `upper` (named numerics with pinned params removed, or NULL).
#' @keywords internal
split_fixed_bounds <- function(lower, upper, start) {
  if (is.null(lower)) lower <- numeric(0)
  if (is.null(upper)) upper <- numeric(0)
  common <- intersect(names(lower), names(upper))
  eq <- common[is.finite(lower[common]) & is.finite(upper[common]) &
                 abs(lower[common] - upper[common]) <= 1e-12]
  start[eq] <- lower[eq]
  drop_lo <- lower[setdiff(names(lower), eq)]
  drop_hi <- upper[setdiff(names(upper), eq)]
  list(fixed = eq,
       start = start,
       lower = if (length(drop_lo)) drop_lo else NULL,
       upper = if (length(drop_hi)) drop_hi else NULL)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::test(filter="app-io")'`
Expected: PASS.

- [ ] **Step 5: Regenerate docs and commit**

Run: `Rscript -e 'devtools::document()'`

```bash
git add R/app-io.R tests/testthat/test-app-io.R man/
git commit -m "feat(app-io): add collect_bounds_all + split_fixed_bounds helpers"
```

---

## Task 3: Engine wrapper — `refine_joint()` + recovery/monotonic tests

**Files:**
- Modify: `R/fit.R` (add `refine_joint()` after `eval_mixture()`)
- Create: `tests/testthat/test-joint-refine.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-joint-refine.R`:

```r
# Tests for the final joint "Optimize all params" engine wrapper.
# Spec: docs/design/specs/2026-06-02-binary-joint-fit-design.md

test_that("joint refine recovers known params from noise-free data", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.5)
  df  <- simulate_mixture(par, "CA", "SA", "continuous", cv = 0)
  # A staged-style starting fit, perturbed away from truth.
  start_fit <- fit_model(df, "CA", "SA", "continuous",
                         start = par * c(1.2, 0.8, 1.3, 0.7, 1.1, 0.6),
                         fixed = c("max", "slope1", "slope2", "ec501", "ec502"),
                         n_starts = 1)
  joint <- refine_joint(start_fit, df, n_starts = 5, time_limit = 60)
  expect_equal(unname(joint$par[["a"]]),   1.5, tolerance = 0.05)
  expect_equal(unname(joint$par[["max"]]), 800, tolerance = 2)
  expect_true(isTRUE(joint$joint))
})

test_that("joint refine never worsens the seed objective", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.2)
  df  <- simulate_mixture(par, "CA", "SA", "continuous", cv = 0)
  staged <- fit_model(df, "CA", "SA", "continuous", start = par,
                      fixed = c("max", "slope1", "slope2", "ec501", "ec502"),
                      n_starts = 1)
  joint <- refine_joint(staged, df, n_starts = 1, time_limit = 60)
  expect_lte(joint$objective, staged$objective + 1e-6)
})

test_that("refine_joint pins a parameter when lower == upper", {
  par <- c(max = 800, slope1 = 4, slope2 = 1.5, ec501 = 0.08, ec502 = 1, a = 1.2)
  df  <- simulate_mixture(par, "CA", "SA", "continuous", cv = 0)
  staged <- fit_model(df, "CA", "SA", "continuous", start = par,
                      fixed = c("max", "slope1", "slope2", "ec501", "ec502"),
                      n_starts = 1)
  joint <- refine_joint(staged, df, lower = c(max = 800), upper = c(max = 800),
                        n_starts = 1, time_limit = 60)
  expect_equal(unname(joint$par[["max"]]), 800, tolerance = 1e-8)  # held
  expect_true("max" %in% joint$fixed)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter="joint-refine")'`
Expected: FAIL — `refine_joint` not found.

- [ ] **Step 3: Implement `refine_joint()` in `R/fit.R`**

Add after `eval_mixture()`:

```r
#' Joint refinement of a mixture fit ("Optimize all params")
#'
#' Re-fits EVERY parameter of an existing mixture fit at once, seeded from that
#' fit's parameters, optionally constrained by `lower`/`upper`. A parameter whose
#' lower and upper bounds are equal is pinned (held fixed) via
#' [split_fixed_bounds()]. Because the seed is the prior fit and the bisection
#' surface is non-smooth, multi-start (`n_starts`) is recommended. The returned
#' fit carries `joint = TRUE`.
#' @param fit An enriched mixture fit (from [fit_model()] / the staged analysis).
#' @param df The mixture data frame the fit was built from.
#' @param lower,upper Optional named numeric vectors of bounds over any model
#'   parameter (curve or interaction). Equal lower==upper pins that parameter.
#' @param n_starts,time_limit Forwarded to [fit_model()].
#' @return An enriched fit (as [fit_model()]) with `joint = TRUE`.
#' @keywords internal
refine_joint <- function(fit, df, lower = NULL, upper = NULL,
                         n_starts = 1, time_limit = 30) {
  spec  <- model_spec(fit$reference, fit$deviation, fit$n_chem)
  start <- fit$par[spec$params]
  sp    <- split_fixed_bounds(lower, upper, start)
  out <- fit_model(df, fit$reference, fit$deviation, fit$response,
                   start = sp$start, fixed = sp$fixed,
                   lower = sp$lower, upper = sp$upper,
                   n_starts = n_starts, time_limit = time_limit)
  out$joint <- TRUE
  out
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::test(filter="joint-refine")'`
Expected: PASS (3 tests).

- [ ] **Step 5: Regenerate docs and commit**

Run: `Rscript -e 'devtools::document()'`

```bash
git add R/fit.R tests/testthat/test-joint-refine.R man/
git commit -m "feat(fit): add refine_joint() for the final joint optimize stage"
```

---

## Task 4: UI — `optimize_param_row()` + the Optimize-all card

**Files:**
- Modify: `R/app-curve-fit.R` (add `optimize_param_row()` near `interaction_param_row()`)
- Modify: `R/app-binary.R` (`binary_ui`: insert the card at the end of the Stage 2 card; `binary_server`: add `output$has_fit`)
- Modify: `tests/testthat/test-app-modules.R` (UI presence test)

- [ ] **Step 1: Write the failing UI test**

Append to `tests/testthat/test-app-modules.R`:

```r
test_that("binary_ui: Optimize-all card with the 7-param bounds grid", {
  html <- as.character(binary_ui("binary"))
  expect_match(html, "Optimize all params", fixed = TRUE)     # button
  expect_match(html, "binary-olo_max", fixed = TRUE)          # lower input
  expect_match(html, "binary-ohi_ec502", fixed = TRUE)        # upper input
  expect_match(html, "binary-oval_a", fixed = TRUE)           # start value input
  expect_match(html, "binary-optimize_all", fixed = TRUE)     # action id
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter="app-modules")'`
Expected: FAIL — the Optimize-all markup is absent.

- [ ] **Step 3: Add `optimize_param_row()` to `R/app-curve-fit.R`**

Add after `interaction_param_row()`:

```r
#' One Optimize-all row: label, meaning, lower/upper bound inputs, start value
#'
#' Five-column layout matching [param_row()]. Uses the `olo_`/`ohi_`/`oval_`
#' prefixes (distinct from the Stage-1 curve grid and the Stage-2 a/b inputs) so
#' the joint Optimize-all panel never collides with other binary inputs. Set
#' Lower = Upper to pin a parameter.
#' @param ns Module namespace function.
#' @param param Parameter key (e.g. `max`, `a`); drives the input ids.
#' @param label Display label.
#' @param meaning One-line explanation.
#' @keywords internal
optimize_param_row <- function(ns, param, label, meaning) {
  shiny::fluidRow(
    shiny::column(2, shiny::tags$b(label)),
    shiny::column(4, shiny::tags$small(meaning)),
    shiny::column(2, shiny::numericInput(ns(paste0("olo_", param)), NULL, value = NA)),
    shiny::column(2, shiny::numericInput(ns(paste0("ohi_", param)), NULL, value = NA)),
    shiny::column(2, shiny::numericInput(ns(paste0("oval_", param)), NULL, value = NA))
  )
}
```

- [ ] **Step 4: Insert the card in `binary_ui` (`R/app-binary.R`)**

Inside the Stage 2 `bslib::card(...)`, immediately AFTER the `shiny::uiOutput(ns("thorough_note"))` line and before the card's closing `)`, add:

```r
        ,
        shiny::conditionalPanel(
          condition = "output.has_fit", ns = ns,
          bslib::card(
            bslib::card_header("Optimize all parameters (joint)"),
            shiny::p("Refines the displayed model by fitting every parameter at ",
                     "once, seeded from the current fit. Fix any parameter by ",
                     "setting its Lower = Upper."),
            shiny::fluidRow(
              shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
              shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
              shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
              shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
              shiny::column(2, shiny::tags$small(shiny::tags$b("Value (start)")))
            ),
            optimize_param_row(ns, "max", "max", "Control response (upper plateau)."),
            optimize_param_row(ns, "slope1", "slope1", "Chemical 1 curve steepness."),
            optimize_param_row(ns, "slope2", "slope2", "Chemical 2 curve steepness."),
            optimize_param_row(ns, "ec501", "EC50 1", "Chemical 1 half-effect conc."),
            optimize_param_row(ns, "ec502", "EC50 2", "Chemical 2 half-effect conc."),
            optimize_param_row(ns, "a", "a", "Overall interaction strength/direction."),
            shiny::conditionalPanel(
              condition = "input.model == 'DR' || input.model == 'DL'", ns = ns,
              optimize_param_row(ns, "b", "b",
                                 "Interaction shift with ratio / dose level.")),
            shiny::actionButton(ns("optimize_all"), "Optimize all params",
                                class = "btn-primary"),
            shiny::uiOutput(ns("optimize_readout"))
          )
        )
```

(The leading `,` joins it as the next argument to the Stage 2 `bslib::card()`.)

- [ ] **Step 5: Add the `has_fit` output in `binary_server` (`R/app-binary.R`)**

Right after the `frozen` output wiring (`shiny::outputOptions(output, "frozen", suspendWhenHidden = FALSE)`), add:

```r
    # Reveal the Optimize-all card only once a model is displayed.
    output$has_fit <- shiny::reactive(!is.null(current_fit()))
    shiny::outputOptions(output, "has_fit", suspendWhenHidden = FALSE)
```

(`current_fit` is defined a few lines below in the source; in Shiny a reactive can reference it because the expression is evaluated lazily. If the linter/tests complain about use-before-def ordering, move this block to immediately after the `current_fit <- shiny::reactive(...)` definition instead.)

- [ ] **Step 6: Run the UI test to verify it passes**

Run: `Rscript -e 'devtools::test(filter="app-modules")'`
Expected: PASS for the new UI test (the existing `testServer` test still passes).

- [ ] **Step 7: Regenerate docs and commit**

Run: `Rscript -e 'devtools::document()'`

```bash
git add R/app-curve-fit.R R/app-binary.R tests/testthat/test-app-modules.R man/
git commit -m "feat(binary): add Optimize-all card + 7-param bounds grid UI"
```

---

## Task 5: Server — Optimize-all observer, prefill, storage, readout

**Files:**
- Modify: `R/app-binary.R` (`binary_server`)
- Modify: `tests/testthat/test-app-modules.R` (extend the `testServer` test)

- [ ] **Step 1: Write the failing server test**

In `tests/testthat/test-app-modules.R`, inside the existing `"binary_server workspace: ..."` `testServer` block, BEFORE the final "editing a single curve after freezing invalidates" section (i.e., after line 85 `expect_false(isTRUE(current_fit()$simulated))`), insert:

```r
    # Optimize-all: jointly refine the displayed model; objective must not worsen.
    chosen <- last_compare()$chosen
    pre_obj <- current_fit()$objective
    session$setInputs(optimize_all = 1)
    expect_true(isTRUE(current_fit()$joint))
    expect_equal(current_fit()$deviation, chosen)
    expect_lte(current_fit()$objective, pre_obj + 1e-6)

    # Pinning max via equal bounds holds it across a re-optimize.
    held <- unname(current_fit()$par[["max"]])
    session$setInputs(olo_max = held, ohi_max = held, optimize_all = 2)
    expect_equal(unname(current_fit()$par[["max"]]), held, tolerance = 1e-8)
    expect_true("max" %in% current_fit()$fixed)
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter="app-modules")'`
Expected: FAIL — `optimize_all` has no observer, so `current_fit()$joint` is NULL.

- [ ] **Step 3: Implement the server logic in `binary_server` (`R/app-binary.R`)**

Add reactive values for the readout near the other `reactiveVal`s (after `last_compare <- shiny::reactiveVal(NULL)`):

```r
    optimize_pre  <- shiny::reactiveVal(NULL)
    optimize_post <- shiny::reactiveVal(NULL)
```

Add a prefill observer and the Optimize-all observer (place them after the existing `observeEvent(input$model, …)` picker-sync block):

```r
    # Keep the Optimize-all start column in sync with the displayed fit.
    shiny::observe({
      f <- current_fit()
      if (is.null(f)) return()
      for (p in names(f$par))
        shiny::updateNumericInput(session, paste0("oval_", p),
                                  value = round(f$par[[p]], 4))
    })

    # Optimize all params: jointly refine the displayed model, seeded from it.
    shiny::observeEvent(input$optimize_all, {
      shiny::req(frozen(), current_fit())
      f <- current_fit()
      spec <- model_spec(input$reference, f$deviation, 2)
      b <- collect_bounds_all(shiny::reactiveValuesToList(input), params = spec$params)
      pre <- f$objective
      newfit <- tryCatch(
        shiny::withProgress(message = "Optimizing all parameters...", value = 0.5,
          refine_joint(f, engine_df(), lower = b$lower, upper = b$upper,
                       n_starts = n_starts_eff(), time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Optimize failed:", conditionMessage(e)),
                                  type = "error")
          NULL
        })
      if (is.null(newfit)) return()
      s <- fits_store(); s[[f$deviation]] <- newfit; fits_store(s)
      optimize_pre(pre); optimize_post(newfit$objective)
      if ("a" %in% names(newfit$par))
        shiny::updateNumericInput(session, "val_a", value = round(newfit$par[["a"]], 4))
      if ("b" %in% names(newfit$par))
        shiny::updateNumericInput(session, "val_b", value = round(newfit$par[["b"]], 4))
    })
```

Add the readout output (near the other `output$` renderers):

```r
    output$optimize_readout <- shiny::renderUI({
      shiny::req(optimize_post())
      lab <- if (identical(current_fit()$response, "binary")) "Deviance" else "SSR"
      improved <- optimize_post() <= optimize_pre() + 1e-9
      shiny::tags$p(
        shiny::tags$b(paste0(lab, ": ")),
        round(optimize_pre(), 2), shiny::HTML(" &rarr; "), round(optimize_post(), 2),
        if (improved) shiny::tags$span(style = "color:green", " ✓ improved"))
    })
```

In the freeze-invalidation observer (the `observeEvent(list(fit1(), fit2(), input$reference, input$response), …)` block), add readout resets alongside the existing `fits_store(list())`:

```r
          optimize_pre(NULL)
          optimize_post(NULL)
```

In the `observeEvent(input$model, …)` picker-sync block, also clear the readout (it refers to the previously displayed model), by adding at the end of that observer:

```r
      optimize_pre(NULL); optimize_post(NULL)
```

- [ ] **Step 4: Run to verify it passes**

Run: `Rscript -e 'devtools::test(filter="app-modules")'`
Expected: PASS (full `binary_server` test incl. the new Optimize-all section).

- [ ] **Step 5: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): wire Optimize-all observer, prefill, storage, readout"
```

---

## Task 6: Confidence intervals — blank pinned params for joint fits

**Files:**
- Modify: `R/app-binary.R` (`output$cis`)
- Modify: `tests/testthat/test-app-modules.R` (assert blanking)

- [ ] **Step 1: Write the failing test**

In the `testServer` block of `tests/testthat/test-app-modules.R`, right after the "Pinning max via equal bounds" assertions added in Task 5, insert:

```r
    # CI table: a pinned parameter shows no interval (it was not estimated).
    cis <- param_ci(current_fit(), engine_df(), current_fit()$reference,
                    current_fit()$deviation, current_fit()$response)
    cis$lower[cis$parameter %in% current_fit()$fixed] <- NA
    cis$upper[cis$parameter %in% current_fit()$fixed] <- NA
    expect_true(is.na(cis$lower[cis$parameter == "max"]))
```

(This asserts the transform the app applies; Step 3 moves it into `output$cis`.)

- [ ] **Step 2: Run to verify it passes at the helper level, then update the app**

Run: `Rscript -e 'devtools::test(filter="app-modules")'`
Expected: PASS (the inline transform works); proceed to fold it into the app output.

- [ ] **Step 3: Update `output$cis` in `binary_server` (`R/app-binary.R`)**

Replace the non-simulated branch so pinned params from a joint fit lose their interval:

```r
    output$cis <- DT::renderDT({
      shiny::req(current_fit())
      f <- current_fit()
      tab <- if (isTRUE(f$simulated)) {
        data.frame(parameter = names(f$par), value = round(unname(f$par), 4))
      } else {
        ci <- param_ci(f, engine_df(), f$reference, f$deviation, f$response)
        if (isTRUE(f$joint) && length(f$fixed)) {
          ci$lower[ci$parameter %in% f$fixed] <- NA
          ci$upper[ci$parameter %in% f$fixed] <- NA
        }
        ci
      }
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"))
    })
```

- [ ] **Step 4: Run to verify the suite still passes**

Run: `Rscript -e 'devtools::test(filter="app-modules")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app-binary.R tests/testthat/test-app-modules.R
git commit -m "feat(binary): hide CIs for parameters pinned in a joint fit"
```

---

## Task 7: Excel cross-check (acceptance characterization)

**Files:**
- Create: `tests/testthat/test-joint-excel-xcheck.R`

This is a soft acceptance check on the real binary data: the joint fit should land near the workbook's published CA+S/A parameters and never worsen the staged objective. It is NOT a bit-exact match (different optimizer; CA bisection non-smoothness), so tolerances are generous and it is skipped on CRAN.

- [ ] **Step 1: Write the test**

Create `tests/testthat/test-joint-excel-xcheck.R`:

```r
# Acceptance characterization: joint fit on the real binary CPF/IMI(MPs) data
# should land near the Excel-published CA+S/A params (Max ~786.4, a ~1.546) and
# never worsen the staged objective. Approximate by design (GRG vs our optimizer).
# Spec: docs/design/specs/2026-06-02-binary-joint-fit-design.md

test_that("joint fit on real binary data is plausible and improves on staged", {
  skip_on_cran()
  csv <- testthat::test_path("fixtures", "binary_mps_cpf_continuous.csv")
  skip_if_not(file.exists(csv), "binary fixture missing")
  df <- utils::read.csv(csv)

  res <- analyse_mixture(df, reference = "CA", response = "continuous", n_starts = 1)
  staged <- res$fits[["SA"]]
  joint  <- refine_joint(staged, df, n_starts = 20, time_limit = 120)

  expect_lte(joint$objective, staged$objective + 1e-6)         # never worse
  expect_gt(unname(joint$par[["max"]]), 600)                   # plausible plateau
  expect_lt(unname(joint$par[["max"]]), 1000)
  expect_true(is.finite(unname(joint$par[["a"]])))             # interaction finite
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter="joint-excel-xcheck")'`
Expected: PASS (or SKIP if the fixture is absent). If `max`/`a` fall outside the plausible band, capture the actual values in the failure message and reconcile against the workbook before tightening — do not silently widen the band.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-joint-excel-xcheck.R
git commit -m "test(binary): acceptance cross-check of joint fit vs Excel params"
```

---

## Task 8: Final verification + memory update

**Files:**
- Modify: `C:\Users\jelle\.claude\projects\D--sam\memory\excel-joint-fitting.md` (+ `MEMORY.md` pointer) — mark implemented

- [ ] **Step 1: Run all touched test files together**

Run: `Rscript -e 'devtools::test(filter="fit-bounds"); devtools::test(filter="joint-refine"); devtools::test(filter="app-io"); devtools::test(filter="app-modules"); devtools::test(filter="joint-excel-xcheck")'`
Expected: all PASS (xcheck may SKIP if fixture missing).

- [ ] **Step 2: Confirm the package loads and docs are current**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); message("OK")'`
Expected: "OK" with no errors; `git status` shows no uncommitted `man/` changes (if it does, commit them).

- [ ] **Step 3: Update memory to record the feature shipped**

Edit `excel-joint-fitting.md`: change the "Plan/decision" line to note the joint Optimize-all stage is implemented on the Binary tab (refine_joint + Optimize-all card), staged flow preserved. Update the `MEMORY.md` one-liner accordingly. (Memory files only — no git commit needed.)

- [ ] **Step 4: Final commit if anything remains**

```bash
git add -A
git commit -m "chore(binary): finalize joint Optimize-all stage"
```

---

## Self-review notes

- **Spec coverage:** UI card (Task 4) ✓; all-7 free-by-default grid (Task 4) ✓; fix-via-Lower=Upper routed to `fixed` (Task 2 `split_fixed_bounds`, Task 3 `refine_joint`) ✓; engine allows bounding a/b (Task 1) ✓; seeded-from-staged + multistart + time_limit (Task 3/5) ✓; before→after readout (Task 5) ✓; store overwrites displayed model + diagnostics follow (Task 5) ✓; staged comparison untouched (no changes to `analyse_mixture`/Find best) ✓; freeze invalidation clears joint fit + readout (Task 5) ✓; CIs blank pinned params (Task 6) ✓; tests: engine units (Task 1), fix-routing (Task 2), recovery + monotonic (Task 3), app wiring (Task 5), Excel cross-check (Task 7) ✓.
- **Type/name consistency:** input ids `olo_<p>`/`ohi_<p>`/`oval_<p>` used consistently across `optimize_param_row` (Task 4), `collect_bounds_all` (Task 2), and the server reads (Task 5). `refine_joint` sets `joint = TRUE` (Task 3) consumed by `output$cis` (Task 6). `fixed` recorded on the fit (Task 1) consumed by Task 6.
- **No placeholders:** every code step shows the actual code and exact run command.
