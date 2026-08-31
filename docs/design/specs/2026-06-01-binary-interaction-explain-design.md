# Binary tab — make the interaction fit legible: design

**Date:** 2026-06-01
**Status:** ✅ Approved design (pre-implementation). To be implemented on branch
`binary-explain-fit` off `main`.
Builds on the staged Binary tab (`R/app-binary.R`,
`docs/design/specs/2026-06-01-binary-staged-ux-design.md`) and the engine's
deviation models (`R/fit.R`, `R/mix-response.R`, `select_parsimonious()` in
`R/compare.R`).

---

## 1. Purpose

After the staged rework, the Binary tab's **Freeze → interaction fit** step is
opaque. Two problems:

1. **It feels like magic.** Clicking *Freeze* runs a fit and four models appear
   with a "chosen" one, but nothing explains *what* is being fitted (only the
   interaction terms `a`/`b`, on the now-fixed curves) or what the `a`/`b`
   parameters and the four deviation models *mean*.
2. **The controls are disconnected.** `alpha`, `n_starts`, the Thorough
   checkbox, and `time_limit` sit in a global sidebar "Advanced" accordion, far
   from the one moment they act (the `analyse_mixture()` call triggered by
   *Freeze*).

This redesign makes the step legible: relocate the controls to the step,
explain the step inline, and show the relevant `a`/`b` formula for the model
being viewed.

## 2. Decisions (locked)

| Decision | Choice |
|---|---|
| Scope | **Binary tab only** this round. The Single tab's dropped "About this model" panel is a separate, tracked follow-up (see §8). |
| Controls | **Keep all four** (`alpha`, `n_starts`, Thorough, `time_limit`); **relocate** them from the sidebar to the Freeze step. No removal, no behaviour change — input ids stay identical so the server logic is untouched. |
| Explanation style | **Inline at the relevant location** (not a separate/collapsible "About" panel): a short sentence at the Freeze step, and a per-model formula block in Stage 2. |
| Formula display | Updates with the **displayed model** (the Stage 2 picker) and the **chosen reference** (CA/IA). |
| Implementation shape | Formula markup comes from a **pure helper** `interaction_help(reference, deviation)`, unit-testable on its own (mirrors the old `model_help_single()`). |

## 3. Layout

**Stage 1 — at the Freeze checkpoint** (below the two curve panels):

```
Stage 1 · Single curves
  [ chem 1 curve panel ]      [ chem 2 curve panel ]

  Freezing locks both curves above and fits only the interaction terms
  (a, b) to the mixture rows — then compares four models (no interaction →
  S/A → dose-ratio → dose-level) and flags the most parsimonious.

  Fit options:  α [0.05]   n_starts [1]   ☐ Thorough fit   time_limit (s) [30]
  [ Freeze curves → fit interactions ]
  <freeze_note>
```

**Stage 2 — the formula for the displayed model** (next to the picker, above
the comparison table):

```
Stage 2 · Interaction models
  Model to display: [ DR ▾ ]

  Dose-ratio dependent (DR):   F = (a + Σ bᵢ·zᵢ) · ∏z
    a — overall strength/direction of the interaction
    b — how the interaction shifts with the mixture ratio
    (zᵢ = each chemical's toxic-unit share; F adjusts the CA baseline)

  [ comparison table: model / parent / chi / df / p ]
```

Switching the picker (reference / SA / DR / DL) updates the formula block to
match; switching the **reference** (CA/IA) in the header updates the DL form
(and the "adjusts the CA/IA baseline" wording).

## 4. The explanation content

A short, always-visible intro sentence at the Freeze step (static markup):

> "Freezing locks both curves above and fits only the interaction terms (a, b)
> to the mixture rows, then compares four models (no interaction → S/A →
> dose-ratio → dose-level) and flags the most parsimonious."

`interaction_help(reference, deviation)` returns the per-model block. Content
per deviation (the interaction factor `F` acts on the toxic-unit shares
`z = TU / ΣTU`, `TU = C / EC50`; `∏z` is the product across chemicals):

- **reference** — *"No interaction term: the mixture follows the {CA|IA}
  reference exactly."*
- **S/A** (`F = a·∏z`) — *"`a` sets the overall strength and direction of the
  interaction."*
- **DR** (`F = (a + Σ bᵢ·zᵢ)·∏z`) — *"`a` = overall interaction; `b` = how it
  shifts with the mixture ratio (which chemical dominates)."* For a binary
  mixture this has one free `b` (Jonker Eq. 8: `b₂` is redundant).
- **DL** — CA: `a·(1 − b·ΣTU)·∏z`; IA: `a·(1 − b·P)·∏z` — *"`a` = overall
  interaction; `b` = how it shifts with the dose level (`P` = the IA-predicted
  effect at the mixture point)."*

**Sign convention** (verified against `R/mix-response.R` + `R/plot.R`: the CA
interaction factor is `exp(F)`, "< 1 synergy, > 1 antagonism", and the IA shift
is monotone in `F`): **`a > 0` → antagonism** (mixture *less* toxic than the
reference predicts), **`a < 0` → synergism** (mixture *more* toxic). Stated once
in the block, e.g. *"a > 0 antagonism, a < 0 synergism."* The implementation
must keep this wording consistent with the engine; if a future engine change
flips the factor, this copy changes with it.

Markup is plain HTML (no MathJax), matching the offline-safe style of the former
`model_help_single()` and the existing `single_model_equation()`. Superscripts
via `tags$sup`, Greek/symbols via HTML entities (`∏` `Σ` `α`).

## 5. Architecture

- **New pure helper `interaction_help(reference, deviation)` in `R/app-binary.R`.**
  Args: `reference` ∈ {CA, IA}, `deviation` ∈ {reference, SA, DR, DL}. Returns a
  `shiny::tagList` of the formula + `a`/`b` meaning for that combination. No
  Shiny reactivity, no fit object — pure, unit-testable.
- **`binary_ui` changes:**
  - Remove the sidebar `bslib::accordion("Advanced", …)` block entirely.
  - The sidebar keeps response, reference, template, file, and errors.
    `thorough_note` moves to the Fit options group.
  - Add, in the Stage 1 card after the curve columns and before the Freeze
    button: a static intro `shiny::p(...)` (the §4 sentence) and a **"Fit
    options"** group holding the four controls — same input ids (`alpha`,
    `n_starts`, `thorough`, `time_limit`) and same defaults as today.
  - Add, in the Stage 2 card next to the `model` picker: a
    `shiny::uiOutput(ns("interaction_help"))` block.
- **`binary_server` changes:**
  - Add `output$interaction_help <- renderUI({ interaction_help(input$reference,
    <displayed deviation>) })`, where the displayed deviation is the same value
    `shown_fit()` resolves (picker value, else `res_r()$chosen`). It must render
    only once frozen (`req(frozen())`), consistent with the rest of Stage 2.
  - `output$thorough_note` is unchanged (only its UI position moves).
  - No change to `res_r`, the freeze/invalidation logic, or how the four inputs
    are read — they keep their ids, so `analyse_mixture(...)` is called exactly
    as now.

## 6. Data flow

```
header reference (CA/IA) ─┐
Stage 2 picker (input$model, else res_r()$chosen) ─┴─► interaction_help(ref, dev)
                                                          → renderUI block (Stage 2)

Fit options (alpha/n_starts/thorough/time_limit, unchanged ids)
   └─ read by res_r() on Freeze → analyse_mixture(...)   (unchanged)
```

## 7. Testing

- **`interaction_help()` unit tests** (no Shiny server): for each
  `reference × deviation` pair, the returned HTML contains the right formula
  fragment and `a`/`b` description — e.g. `SA → "a·∏z"`, `DR → "Σ"` and "mixture
  ratio", `DL` CA vs IA differ (`ΣTU` vs `P`), `reference → "No interaction"`.
  One test asserts the sign-convention sentence is present.
- **`testServer` (binary):** with a frozen fixture fit, `output$interaction_help`
  follows the picker — selecting `DR` yields the DR formula; the displayed block
  matches `shown_fit()$deviation`. (Access the rendered text via the existing
  `shown_fit()`/`input$model` mechanics already used in the staged-flow test.)
- **Regression:** the existing staged-flow binary test stays green — control ids
  are unchanged, so freeze → four-model fit → picker override is unaffected.
- **UI smoke:** `binary_ui("binary")` builds and its HTML contains the Fit
  options labels and the Freeze-step intro sentence.

## 8. Out of scope / tracked follow-ups

- **Single tab "About this model" panel** (`model_help_single()` + accordion +
  its test, committed in `ace9b8a`, dropped during the curve-fit-panel refactor):
  not restored here per the scoping decision. The Single tab still shows the
  model equation, the SSR formula, and the Continuous/Quantal helptext, so it is
  not bare. Restore in a later round (ideally in this same inline style).
- **Cutting any control** (the YAGNI option of dropping `n_starts`/Thorough/
  `time_limit` as vestigial in the staged design): explicitly *not* done — the
  decision is to keep all four and relocate.
- The **scroll fix** (`page_navbar(fillable = FALSE)` in `R/app-run.R`) is a
  separate, already-made working-tree change; not part of this spec.
