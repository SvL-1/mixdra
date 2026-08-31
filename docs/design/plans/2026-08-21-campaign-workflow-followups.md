# Campaign workflow — follow-ups

**Date:** 2026-08-21
**Branch:** `campaign-workflow`
**Source:** the final whole-branch review and its scoped re-review.

None of these blocks the branch. Each was found by review, judged a follow-up
rather than a blocker, and is recorded here so it is not lost when the review
workspace is deleted. Ordered roughly by value.

## 1. `curve_fit_server` keeps the previous upload's fit

`R/app-curve-fit.R:113` — `current_fit` is never cleared when `fit_df()` changes.
After a new upload the Singles panels draw the **old file's fitted curve over the
new file's points**, and `output$diagnostics` (`:173-179`) prints the old fit's
`ssr` beside `nrow(fit_df())` from the new data — a mixed-provenance readout.

This is the same class as C3 from the final review, and the only surviving
instance of it. It is display-only: nothing reads `store$singles`, and it cannot
reach `store$base` or any campaign number. It predates this branch and affects
the standalone Single Stressor tab equally.

Fix is one line: `observeEvent(fit_df(), current_fit(NULL), ignoreInit = TRUE)`.

## 2. Roxygen drift from the deleted stale banner

`campaign_bump_base()`'s roxygen (`R/app-campaign.R:17-25`) and its generated
`.Rd` still say downstream stages "show a stale banner". The banner was removed;
§6 of the design doc has been amended to record that, but the roxygen has not.
Fix the roxygen and re-`document()` — fixing the `.Rd` alone regresses.

## 3. The quantal→engine mapping now exists in two places

`R/app-singles.R:29` and `R/app-campaign.R:176` each map `"quantal"` to the
engine's `"binary"`, and both fall through silently to `"continuous"` for any
unexpected value. One `campaign_engine_response()` helper would remove the drift
risk. Pairs naturally with item 7.

## 4. C3's own fix is not tested in isolation

The end-to-end test would still pass with only C1's clear in place, because
`base()` also goes `NULL` on a new upload. A direct `pair_workspace_server` test
that changes `fit_df` while holding `base` constant would pin the
self-protection the fix actually added.

## 5. Ternary gating message precedence

The "no ternary rows" branch pre-empts the IA and quantal explanations, so an IA
campaign that also lacks a ternary arm is told the less actionable of the two
facts. Reorder so the reference/response limitation is reported first.

## 6. The sub-navigation is rebuilt on every store write

`R/app-campaign.R` renders the whole navset from `campaign_stage_nav()`, which
depends transitively on most of the store. The reported symptom — losing the
selected tab mid-analysis — is fixed via `id`/`selected`, but the rebuild still
tears down DT scroll position and accordion state. Moving the gating decisions
into per-panel `uiOutput`s would narrow the dependency set. Deliberately deferred
from the fix wave: it restructures the module's core rendering path and would
have landed unreviewed.

## 7. Naming and placement

- `R/app-binary.R` now contains only `interaction_help` and `pair_workspace_*`;
  the name no longer describes it. `R/app-singles.R` (campaign) and
  `R/app-single.R` (scratchpad) differ by one character for unrelated modules.
- The `campaign_*` API is split across three files with no rule:
  `campaign_base`/`campaign_fit_base` in `app-singles.R`; `campaign_pairs`,
  `pair_key`, `pair_has_rows`, `campaign_pairwise`, `campaign_pairs_fitted` in
  `app-campaign.R`; `pair_base` in `app-io.R`.
- `%||%` is defined in `R/app-campaign.R` but used across three modules.

## 8. Dead code left by the restructure

- `assemble_curve_params()` / `assemble_curve_params3()` (`R/app-io.R`) have no
  production caller since the joint fit landed — they are now a test oracle only.
  Delete with their tests, or comment them as a deliberate oracle.
- The `"binary"` and `"ternary"` branches of `upload_schema` / `template_df` /
  `validate_upload` are unreachable from the app (only `"single"` and
  `"campaign"` are used), as is
  `inst/extdata/binary_ca_cpf_imi_fbsa_continuous.csv`.
  `test-app-modules.R` still asserts that fixture "validates as binary
  continuous", testing a stage the app no longer offers.
- `store$ternary` is written nowhere — the ternary result is recomputed by
  `asa_res()`. Delete the field and its row in §5 of the design doc.

## 9. `campaign_stage_nav()` has thin test coverage

It implements the whole of §10 (two-stressor campaigns hiding Ternary, pairs with
no mixture rows, the IA/quantal explanations, the gating messages). The fix wave
added coverage for the no-ternary-rows and tab-persistence cases; the others are
still uncovered.

## 10. Smaller items

- A pair with no mixture rows deadlocks the ternary: `campaign_pairwise()`
  requires all three terms, so the user is told to "fit the interaction models on
  all three pair tabs" when the data makes one impossible. The gate is
  scientifically right; the message is not.
- `analyse_ternary()` is exported and gained `pairwise` as its 5th positional
  argument, before `lower`/`upper`. No in-repo caller is affected, but it is a
  positional break for any external script.
- `read_upload()` is unguarded, so a non-CSV upload throws out of `parsed()`.
  Pre-existing convention, but the campaign is now the app's main entry point.
- On a failed upload the campaign keeps analysing the previous file while showing
  errors: the store keeps the old `raw`, and the sidebar row-count summary still
  describes it.
- The joint fit's `tryCatch` catches `error=` but not `warning=`, so an `optim`
  convergence warning reaches the console on a button press.
- Stale test comment at `test-app-modules.R:781-784` points at "Task 3" for
  version-bump coverage that now lives in the newer test.
