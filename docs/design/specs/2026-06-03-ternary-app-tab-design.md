# Ternary app tab (T3) — design

**Date:** 2026-06-03
**Status:** Approved — ready for implementation plan.
**Subproject:** Ternary (T3, the deferred "app + I/O integration" piece). Engine
(T1) and plotting (T2) are DONE and green; see
`docs/design/specs/2026-05-30-ternary-subproject-design.md` and
`docs/design/specs/2026-06-01-ternary-plotting-design.md`.

## 1. Goal

Add a `Ternary Mixture` tab to the mixdra Shiny app that reuses the **binary
tab's UX skeleton** (sidebar + three stacked stage-cards + a table-as-hub +
click-to-inspect) while honoring the fact that ternary's science is different:
there is no SA/DR/DL model-selection winner. The ternary analysis is the staged
**Advanced-S/A** workflow whose scientific punchline is the **overall-A4 vs
per-ratio-A4 contrast** (the "averaging-out" story).

## 2. What already exists (do not rebuild)

- `analyse_ternary(df, reference, response, …)` → staged ASA fit returning
  `base`, `pairwise` (A1/A2/A3), `A4_overall`, `individual` (per-ratio A4 data
  frame), and `fits`.
- `ternary_effect_table(res, df)` → per-ratio near-EC50 readout: `pred_CA`,
  `pred_SA`, `pred_ASA`, `a4_effect`.
- `plot_isoplane(res, df, n)`, `plot_sigma_tu(res, n)`.
- `classify_rows()`, `ternary_ratio_key()`.
- `curve_fit_ui()` / `curve_fit_server()` — the reusable single-curve panel
  (used as-is for each of the three chemicals).

## 3. Decisions (from brainstorming)

1. **Stage 1 is interactive, binary-style** — three curve panels; the user fits
   each single curve and they are held fixed for the staged fit. (Requires the
   engine `base=` change in §4.)
2. **Stage 2/3 hub = the per-ratio A4 table** — `Overall` row + one row per
   ternary ratio; click a row to inspect it in Stage 3.
3. **Sidebar shows the radios with unsupported options disabled** — Continuous +
   CA enabled; Quantal and IA disabled with a "coming soon" note.
4. **Bundle the ternary fixture as the default example** so the tab is usable on
   open.

## 4. Engine change (small, isolated)

`fit_ternary_asa()` and `analyse_ternary()` gain an optional `base =` argument: a
named numeric vector `max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3`.

- When `base` is supplied, **skip the internal Stage-1 singles fit** and use the
  supplied curves as `base_par`, held fixed in every downstream stage exactly as
  the auto-fit values are today (Stage 2 A1/A2/A3, Stage 3a overall A4, Stage 3b
  per-ratio A4).
- When `base` is `NULL` (default), behaviour is unchanged — existing callers and
  tests are unaffected (back-compat).

This is the ternary analogue of how the binary tab holds the single-chemical
curves fixed via `curve_params()`. It is the only engine change.

## 5. I/O layer (`R/app-io.R`)

- `upload_schema("ternary", response)` — continuous → `c("C1","C2","C3","Res")`;
  quantal → `c("C1","C2","C3","Affected","Exposed")` (defined for completeness;
  unused while quantal is disabled).
- `template_df("ternary", response)` — a small dataset spanning **all tiers**
  (control + three singles + three binaries + at least one ternary ratio's
  dose-series) so the staged fit has data for every stage.
- `validate_upload(df, "ternary", response)` — extend the `stage` `match.arg`;
  add `C3` to the concentration range checks; add a non-fatal warning if the
  upload contains no ternary rows (all three chemicals > 0), since the per-ratio
  A4 step needs them.
- `to_engine_df(df, "ternary")` — identity (already `C1/C2/C3`).
- `marginal_df3(df, chem)` — chem *i*'s single-chemical series: keep rows where
  the **other two** concentrations are 0, drop those columns, rename chem *i*'s
  column to `C1`. (Generalises the existing two-chemical `marginal_df()`; keep
  `marginal_df()` for the binary tab.)
- `assemble_curve_params3(fit1, fit2, fit3)` — shared `max` (mean of the three
  fitted maxes, matching the engine's seeding) + per-chemical `slope1/2/3` and
  `ec50_1/2/3`. **Note the underscore names** (`ec50_1`, not binary's `ec501`):
  these must match what the ternary engine's `base` argument expects.

## 6. Bundled example

Copy the validated fixture
`tests/testthat/fixtures/ternary_fbsa_cpf_imi_continuous.csv` (the one T1 anchors
its validation to; confirm the exact path at plan time) to
`inst/extdata/ternary_ca_fbsa_cpf_imi_continuous.csv`. The server auto-loads it
via `system.file("extdata", …)` until the user uploads, mirroring binary.

## 7. UI module (`R/app-ternary.R`)

`ternary_ui(id)` + `ternary_server(id, meta)`, mirroring `binary_ui`/
`binary_server` structure.

### Sidebar
- Response radio: Continuous (enabled) / Quantal (**disabled** + "coming soon").
- Reference radio: CA (enabled) / IA (**disabled** + note).
- Download template, Upload CSV, example help text, error panel.
- Advanced accordion: `n_starts`, `time_limit`, `thorough` (multi-start).
  **No `alpha`** — there is no LR test / model selection.

### Stage 1 · Single curves
Three `curve_fit_ui` panels (`chem1`, `chem2`, `chem3`), each fed by
`marginal_df3()`, with a reveal note shown until all three are fitted. Gates
Stage 2/3 on `frozen = all three fits present`.

### Stage 2 · Advanced S/A
- "Fit Advanced S/A" button (gated on `frozen`) runs
  `analyse_ternary(engine_df, reference = "CA", response = "continuous",
  base = assemble_curve_params3(...), n_starts, time_limit)` with a progress bar.
- A base + pairwise readout: `max, slope1-3, ec50_1-3, A1, A2, A3`.
- The **hub table**: an `Overall` row (showing `A4_overall`) + one row per
  ternary ratio, columns: *ratio (C1:C2:C3 label), n, individual A4, pred_CA,
  pred_SA, pred_ASA, a4_effect* (the per-ratio columns come from
  `ternary_effect_table()`). Single-row selection; the `Overall` row highlighted.
- Help text explaining the averaging-out interpretation (per-ratio A4 differing
  in sign while overall A4 ≈ 0 means real interaction is being averaged away).

### Stage 3 · Inspect
- `plot_isoplane(res, engine_df)` — 3-D simplex isoplanes + per-ratio EC50
  markers.
- `plot_sigma_tu(res)` — ΣTU-vs-z curves.
- An effect-size readout for the **selected** ratio (pred_CA/SA/ASA + a4_effect
  at its near-EC50 point). Nice-to-have: highlight the selected ratio's EC50
  marker on the isoplane.
- Before the staged fit has run, Stage 3 shows nothing (or the raw simplex point
  cloud); the plots require a `res` object.

## 8. Wiring (`R/app-run.R`)

- `app_ui()`: add `bslib::nav_panel("Ternary Mixture", ternary_ui("ternary"))`
  after the Binary panel.
- `app_server()`: add `ternary_server("ternary", meta)`.
- `meta$chem3` is set wherever `chem1`/`chem2` are (axis labels for the third
  curve panel).

## 9. Data flow

```
upload / bundled example
  → parse → validate_upload(…, "ternary")
  → engine_df (C1,C2,C3,Res)
  → marginal_df3 ×3 → curve_fit_server ×3 → fit1/fit2/fit3
  → assemble_curve_params3 → frozen base
  → [Fit Advanced S/A] → analyse_ternary(engine_df, base = frozen, …) → res
  → ternary_effect_table(res, engine_df) → hub table (Overall + per-ratio rows)
  → row selection → Stage-3 effect readout + isoplane / ΣTU plots
```

Any change to the curves / response / reference clears the stored `res` (the hub
table blanks), exactly as the binary tab clears its interaction store.

## 10. Testing

- **Engine** (`test-fit-ternary-asa.R`): `base=` is used (skips Stage-1 fit,
  downstream A1-A4 computed against the supplied curves) and `base=NULL` is
  unchanged (back-compat).
- **I/O** (`test-app-io.R` or similar, no Shiny): `upload_schema("ternary")`,
  `template_df("ternary")`, `validate_upload(…, "ternary")`, `marginal_df3`,
  `assemble_curve_params3`.
- **Module** (`testServer`, like `test-app-modules.R`): feed the bundled example,
  fit three curves, run Advanced S/A, assert the hub table has an Overall row +
  the expected per-ratio rows and that selection drives the effect readout.
- Keep the existing suite green. Run **filtered** tests only
  (`devtools::test(filter = …)`), per the project's test-running convention.

## 11. Explicitly out of scope

- IA Advanced-S/A and quantal ternary (engine doesn't support them yet — radios
  disabled).
- A binary-style "Optimize all params (joint)" card — the staged Advanced-S/A
  method **is** the intended ternary fit; no joint refine.
- Any change to the binary or single tabs beyond adding `meta$chem3`.
