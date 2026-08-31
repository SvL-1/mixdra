# Design: "About this model" panel on the Single Chemical tab

**Date:** 2026-06-01
**Status:** Approved (brainstorm)
**Scope:** Single Chemical tab of the mixdra Shiny app only.

## Problem

The Single Chemical tab shows a "Parameters" table (`max`, `slope`, `ec50`, `SSR`,
`n`) with no indication of what model those parameters belong to. A reader cannot
tell from the table alone that this is a three-parameter log-logistic curve, what
each parameter means, or that `SSR`/`n` are fit diagnostics rather than curve
parameters. For the domain expert (Sam) the vocabulary is familiar, but the table
is still ambiguous — e.g. `max` reads like a function name, and the direction of
`slope` is not obvious without the equation.

## Goal

Make the tab self-documenting at near-zero cost: show the fitted model's equation
and a one-line gloss per parameter, without cluttering the tab for an expert user
and without introducing a network dependency.

## The model (ground truth)

From `R/single.R`, `analyse_single()` fits a three-parameter log-logistic curve:

```
Y = max / (1 + (C / ec50)^slope)
```

- `max`  — control/baseline response at `C = 0` (upper plateau).
- `slope` — steepness of the decline (Hill slope); `slope > 0` ⇒ response decreases with dose.
- `ec50` — concentration that halves the response (classic EC50).

`SSR` (residual sum of squares) and `n` (number of data points) are reported
alongside as fit diagnostics; they are not parameters of the curve.

## Design

Add a **collapsed** `bslib::accordion` panel titled **"About this model"** to
`single_ui()` in `R/app-single.R`, positioned just above the Parameters card.
It opens on demand, so it stays out of the way for expert users.

### Content (static HTML, independent of any fit)

1. **Equation** — `Y = max / (1 + (C / EC50)^slope)`, rendered as plain HTML
   (`<sup>` for the exponent), with a one-line caption: "three-parameter
   log-logistic; the response falls from `max` at zero dose."
2. **Parameter glossary** — an `htmltools` definition list (`<dl>`):
   - **max** — control/baseline response at `C = 0`.
   - **slope** — steepness of the decline (`> 0` ⇒ decreasing).
   - **EC50** — concentration that halves the response.
3. **Diagnostics note** — "**SSR** = residual sum of squares (goodness of fit);
   **n** = number of data points. These describe the fit, not the curve."

### Rendering choice

Plain HTML / unicode, **not** MathJax. The equation is simple enough to render
with `<sup>`, and avoiding `shiny::withMathJax()` keeps the panel working
**offline** (MathJax pulls from a CDN) and adds no dependency.

### Component boundary

Factor the markup into a small internal helper:

```r
model_help_single <- function() { ... }   # returns an htmltools/bslib tag
```

`single_ui()` calls `model_help_single()` so it stays readable, and the helper
is unit-testable in isolation. No server-side change — the panel is pure markup
and does not depend on whether a fit exists.

## Testing

One `testthat` case (in the existing app-modules test file or a small new file):

- `model_help_single()` returns a Shiny/htmltools tag, and its rendered text
  (via `htmltools::renderTags()$html` or `as.character()`) contains `"EC50"`
  and `"log-logistic"`.

The existing fast suite must stay green; no existing behaviour changes.

## Out of scope

- **Binary tab.** The Binary Mixture tab has a larger parameter set (per-chemical
  curves + deviation parameters). The same pattern can be applied later as a
  follow-up; this change is single-tab only.
- **Parameter table rows.** Left exactly as-is (`max`, `slope`, `ec50`, `SSR`,
  `n`). The glossary carries the clarity; no relabelling.
