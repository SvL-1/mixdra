# Campaign workflow — design

**Date:** 2026-08-20
**Status:** Approved — ready for implementation plan.
**Issue:** [#10 "Workflow suggestion"](https://github.com/SvL-1/mixdra/issues/10).
**Supersedes** the standalone-tab structure introduced by
`docs/superpowers/specs/2026-06-03-ternary-app-tab-design.md` (that tab is kept,
but is re-parented under the campaign and stops fitting its own Stages 1–2).

## 1. Goal

Restructure the Shiny app around a **campaign**: one uploaded dataset covering
2–3 stressors, whose single-stressor curves are fitted **once** and then held
fixed by every downstream stage. Today the same curve can be fitted up to four
times (standalone Single tab, inside each binary tab, inside the ternary tab),
producing several EC50s for one chemical with no indication which is
authoritative. The campaign makes one set of curves the single source of truth.

## 2. Origin — Sam's decisions (issue #10, 2026-08-20)

1. Do not refit singles in every tab; fit once and reuse downstream — **and
   likewise reuse the fitted binaries in the ternary fit**.
2. **One upload**, structured like `FBSA CPF IMI ternary -simplified.xls`:
   single rows carry zeros in the other two concentrations, binary rows carry
   zero in the third.
3. **Sub-tabs** for the three binary pairs, not three top-level tabs.

He also asks for this before further small-issue work, because the new layout
makes the app easier to test.

## 3. What already exists (do not rebuild)

- **The engine is already campaign-shaped.** `fit_ternary_asa()` takes one frame
  containing every stratum, calls `classify_rows()`, and stages off it.
- `marginal_df3(df, chem)` — slices one stressor's single series out of a
  ternary frame (keeps the shared control row).
- `assemble_curve_params()` / `assemble_curve_params3()` — build the frozen
  base-parameter vectors.
- `curve_fit_ui()` / `curve_fit_server()` — the reusable single-curve panel,
  used as-is by the new Singles page.
- `model_spec()` and `selection_chain_order()` are parameterised by `n_chem`,
  and `R/app-binary.R:415` already derives `n_chem` from the data — the n=2/n=3
  generalisation is mostly already paid for.
- The bundled `inst/extdata/ternary_ca_fbsa_cpf_imi_continuous.csv` is already a
  campaign file: 19 control, 75 single, 225 binary, 100 ternary rows.

### The reuse of binaries is exact, not an approximation

In `ca_asa_tri()` (`R/models-ternary-asa.R:43-68`) each binary branch activates
exactly one interaction term — `exp(A1·z1·z2)` when `C3 == 0`, `exp(A2·z1·z3)`
when `C2 == 0`, `exp(A3·z2·z3)` when `C1 == 0`. The binary S/A model's core is
`a * prod(z)` under `exp()` (`R/mix-response.R:48,63`) — the same quantity.

Therefore, **given one shared base**, the pooled Stage-2 fit over all binary rows
decomposes exactly into the three per-pair fits, and each pair's S/A `a` **is**
its ternary A-term. The campaign store is what makes the shared-base
precondition true. This equivalence is the linchpin of Sam's request and is
guarded by a test (§11).

## 4. Decisions (from brainstorming, 2026-08-20)

1. **The campaign replaces the mixture tabs.** Navbar becomes Introduction /
   Single Stressor / Campaign. The Single Stressor tab survives as a
   **non-campaign scratchpad** (a lone heat-stress experiment is a legitimate
   use and Sam asked for "Stressor" wording precisely for it).
2. **A campaign is 2 *or* 3 stressors.** A two-stressor campaign is one pair
   with no ternary stage. This covers the 11 standalone binary-pair datasets
   through the same code path, so they cannot disagree with the campaign
   numbers.
3. **DR/DL pairs still hand their S/A `a` to the ternary.** The ternary formula
   has no DR/DL equivalent, and the staged chain fits S/A for every pair anyway.
   Where a pair's selected winner is not S/A, the pair sub-tab and the ternary
   stage both say so, so the approximation is visible rather than silent.
4. **CA/IA is a campaign-wide choice; the ternary stage is CA-only.** Choosing
   IA leaves the ternary sub-tab in an explicit "not yet supported for IA"
   state, replacing today's disabled radio with no explanation
   (`R/app-ternary.R:18`).
5. **Staleness is tracked by version stamping**, not by auto-recompute — see §6.

## 5. State — the campaign store

`meta` (today a flat bag of name/unit strings) is replaced by a `campaign`
`reactiveValues`:

| Field | Meaning |
|---|---|
| `raw` | the validated upload |
| `n_chem` | 2 or 3, detected from the file |
| `response` | `"continuous"` / `"quantal"` |
| `reference` | `"CA"` / `"IA"`, campaign-wide |
| `names`, `units` | per stressor, from the Introduction tab |
| `singles` | per-stressor fit objects |
| `pairs` | keyed `"12"` / `"13"` / `"23"` |
| `ternary` | the ASA result |

`base` stays **derived** — `assemble_curve_params*()` over `singles`, cheap
enough to be a plain reactive rather than stored state.

**The flat `chem1..3` / `unit1..3` fields are retained** alongside the new
campaign fields. `axis_label()` (`R/app-single.R:16-23`) reads them by name, and
`curve_fit_server(chem_field = …)` passes them through. Keeping them is what
allows the curve-fit panel to be reused **genuinely unchanged** (§9); `names` and
`units` in the table above are the campaign-level view of the same values, not a
replacement for them.

## 6. Staleness — version stamping

Expensive fits in this app are deliberately imperative: `fits_store` is a
`reactiveVal` filled by a button, not a derived reactive, because refitting four
models is slow. That is the right call, but it means a downstream stage can
silently belong to a superseded curve.

- `base_version` increments whenever any single fit changes.
- Every pair result records the `base_version` it was fitted at; the ternary
  result records both `base_version` and the pair versions.
- When a recorded stamp no longer matches, the stage shows a **stale banner with
  a refit button and greys its numbers**. It never silently presents stale
  numbers as current, and it never re-runs an expensive fit unasked.

A pure reactive chain was rejected: it would refit 3 pairs × 4 models on any
upstream twitch.

## 7. I/O layer (`R/app-io.R`)

- **`pair_df(df, i, j)` — new.** Mirrors `marginal_df3()`: keep rows where the
  third concentration is 0, drop that column, rename `Ci`→`C1`, `Cj`→`C2`.
- **`campaign_n_chem(df)` — new.** 2 or 3 from the concentration columns present.
- **`upload_schema("campaign", response)`** — `C1`, `C2`, optional `C3`, plus the
  response columns.
- **`validate_upload()`** gains a campaign branch: at least 4 distinct
  concentrations per single series; a pair with no rows **disables that sub-tab**
  rather than erroring (unlike the ternary tab's hard "no ternary rows" error,
  which stays).
- **`classify_rows()`** is reused from `R/fit-ternary-asa.R`. Its two-stressor
  semantics must be verified and extended if needed: with only `C1`/`C2`
  present, a row with both > 0 is `"binary"`, both == 0 is `"control"`, exactly
  one > 0 is `"single"`, and the `"ternary"` class never occurs.
- The existing ternary template is promoted to the campaign template rather than
  rewritten.

## 8. Engine change (small, isolated)

`fit_ternary_asa(..., pairwise = NULL)` — when supplied, **skip Stage 2** and
hold `A1`/`A2`/`A3` fixed at the given values. It mirrors the existing `base=`
argument exactly (`R/fit-ternary-asa.R:46-72`), including its validation shape,
and `analyse_ternary()` forwards it. Mapping: pair `"12"`→`A1`, `"13"`→`A2`,
`"23"`→`A3`. Back-compatible: `NULL` keeps today's behaviour.

No other engine change. `analyse_mixture()` is **not** given a `base=` argument
— the binary path already achieves the same via
`fit_model(start = curve_params(), fixed = names(curve_params()))`
(`R/app-binary.R:411-413`).

## 9. Modules

| File | Change |
|---|---|
| `R/app-campaign.R` | **new** — owns the upload, validation, `n_chem`, the campaign-wide response/reference choice, the store, and the sub-navigation |
| `R/app-singles.R` | **new** — 2–3 `curve_fit` panels on one page over `marginal_df*` slices; writes fits into the store; reuses `curve_fit_server` unchanged |
| `R/app-binary.R` | **decomposed** — upload and the two Stage-1 curve panels move out; the remainder becomes `pair_workspace_ui/server(id, pair_df, base, reference, response)`: staged selection loop, results table, joint refine, plots, inspect. 629 lines → target well under 400 |
| `R/app-ternary.R` | consumes frozen `base` **and** `pairwise` instead of fitting Stages 1–2 |
| `R/app-single.R` | untouched — the scratchpad |
| `R/app-intro.R` | feeds names/units into the campaign store |
| `R/app-run.R` | navbar: Introduction / Single Stressor / Campaign |

## 10. Gating

Sub-nav items **disable with an inline reason** rather than vanishing ("Fit all
single-stressor curves first" — two or three depending on the campaign, so the
message is generated from `n_chem`). The ternary sub-tab is hidden entirely when
`n_chem == 2`. Pairs whose winner is DR/DL carry a visible note that the ternary
uses their S/A term.

## 11. Testing

Pure helpers (`pair_df`, `campaign_n_chem`, the campaign `validate_upload`
branch) are unit-tested alongside the existing `app-io` tests. Two oracles carry
the weight:

1. **Regression** — on the bundled fixture, a pair fitted through the campaign
   path reproduces today's binary-tab numbers given the same base.
2. **The linchpin (§3)** — per-pair `a` fitted from `pair_df()` equals
   `fit_ternary_asa()`'s pooled Stage-2 `A1`/`A2`/`A3` at the same base, to
   optimizer tolerance. If this ever fails, Sam's reuse is no longer exact.

Module smoke tests follow `tests/testthat/test-app-modules.R`. Run filtered
(`devtools::test(filter = …)`), never the full suite.

## 12. Phasing

One branch, four reviewable phases:

1. I/O helpers + campaign store + one upload + Singles page — **showable to Sam**
2. Pair-workspace decomposition + the three sub-tabs
3. Ternary on frozen base + pairwise (engine `pairwise=`)
4. Remove the old top-level Binary/Ternary **nav entries**; README + vignette.
   The modules themselves are re-parented under the campaign, not deleted —
   only their standalone entry points and their own upload/Stage-1 paths go.

## 13. Explicitly out of scope

- Quantal or IA support in the ternary stage (unchanged from today).
- Any change to the fitting mathematics. The campaign changes *which* parameters
  are held fixed and *where they come from*, never how a fit is computed.
- The open small issues (#7, #8) and the `EC50 ≤ top tested dose` bound offered
  in #8 — separate work, deliberately not bundled.
