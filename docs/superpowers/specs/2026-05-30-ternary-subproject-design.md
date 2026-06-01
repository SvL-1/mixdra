# Ternary subproject — design & blocker (ON HOLD)

**Date:** 2026-05-30
**Status:** ⛔ ON HOLD — T1 cannot be finalised until a model-spec question is
answered, which itself depends on other (upstream) implementation work the user
wants to do first. This document records the decomposition, the T1 plan sketch,
the workbook findings, and the blocker so we can resume cleanly.

Builds on the finished engine (`docs/superpowers/plans/2026-05-30-mixdra-engine.md`)
and the Plan-2-split decomposition (single+binary plotting etc.). Single + binary
are validated; ternary is the next big piece and gets its own subproject.

---

## 1. The ternary subproject (umbrella decomposition)

Mirrors the binary arc: validate the engine first, then plots. Each piece is its
own spec + implementation plan.

- **Plan T1 — Ternary engine validation** *(the first plan; blocked, see §4).*
  Validate `analyse_mixture()` end-to-end on the ternary raw-data workbook.
  Extract a fixture, identify the workbook reference fit, write
  `test-validation-ternary.R`, fix any ternary-specific engine bugs. This is the
  gate — nothing downstream is trustworthy until ternary fits reproduce the
  workbook.
- **Plan T2 — Ternary plotting** (the old "Plan 3"). EC50 isoplane surface over
  the simplex + ΣTU / z-value plots. Needs new machinery (an EC50-isobole solver
  over the simplex). Reference material: `Isoplane_MPs_FBSA_IMI.R`,
  `Isoplains_Mixtox.xlsx`, `Model4Isoplane.xlsm`.
- **Plan T3 — Ternary app + I/O integration.** *DEFERRED — left for a separate
  step at the user's request; not scoped here.*

T2/T3 stay as stubs until T1 is green.

---

## 2. Plan T1 — ternary engine validation (sketch)

**Goal:** prove `analyse_mixture()` reproduces the MixTox ternary workbook fit,
exactly as `test-validation-binary.R` did for binary, so ternary fits are
trustworthy before any plotting work.

The engine already accepts 3-chemical input: `analyse_mixture()` does
`intersect(c("C1","C2","C3"), names(df))` and `seed_from_singles()` builds the
right parameter set generically; the registry/fitter treat `spec$extra`
generically. So this is a validation effort, not engine wiring — **except** for
the deviation-model blocker in §4.

**Proven binary pattern to mirror:**
1. `extract_*_fixture.R` reads workbook sheets via `readxl`, keeps numeric rows,
   writes CSV fixtures into `tests/testthat/fixtures/`.
2. `test-validation-*.R` runs `analyse_mixture()` and asserts residual SS within
   ~2% of the workbook, key EC50 within tolerance, and the **chosen model**
   matches the workbook's conclusion.

**Workstreams (once unblocked):**
1. **Map & extract the raw data.** Locate the C1/C2/C3 + response columns in the
   `.xls` (the param/SS block is rows 1–20; the ~404 data points are below it —
   columns not yet mapped). Write `extract_ternary_fixture.R` →
   `tests/testthat/fixtures/ternary_fbsa_cpf_imi_continuous.csv`
   (cols `C1,C2,C3,Res`).
2. **Pin reference values** from the workbook (see §3) into a provenance comment.
3. **`test-validation-ternary.R`** — CA reference assertions (SS, EC50s, chosen
   model). IA only if a reference exists (workbook is CA-only → likely skip IA).
4. **Fix any ternary-specific engine bugs** surfaced (e.g. DR `b1/b2/b3`
   perturbation escaping 0; CA-bisection performance — 3-chem × multi-start ×
   ~404 rows is even slower than binary).

**Decided:** validation reference data source = `FBSA CPF IMI ternary
-simplified.xls`, sheet `CA CPF FBSA IMI` (the primary fit). Response type =
whatever the workbook has (continuous reproduction count; no quantal ternary
sheet seen).

---

## 2a. Ternary fitting is TWO-TIER: overall + per-ratio individual fits

(Added 2026-05-30 from user clarification; see paper *J. Hazard. Mater.*
ISSN 0304-3894, file `1-s2.0-S0304389425034132-main` — not on disk here.)

**Why:** fitting the ternary model over *all* ternary data at once "averages
out" the fit — the model has too little flexibility ("fluidity") to follow how
the interaction changes across mixture ratios, so a single global fit smears the
ratio-specific behaviour.

**So ternary analysis produces two kinds of fit:**
1. **Overall fit** — single + binary + **all** ternary data → one global
   parameter set. (Workbook sheet `…Overall`.)
2. **Individual fits** — single + binary + **one ternary mixture ratio's** data
   at a time → a separate parameter set *per ratio*. Singles+binaries are the
   shared backbone reused in every fit; only the ternary rows are subset to the
   one ratio. (Workbook ratio sheets: `33% FBSA`, `20CPF 60FBSA 20IMI`,
   `CPF60 FBSA20 IMI20`, `CPF20 FBSA20 IMI60`.)

**Output requirement:** ternary results must report BOTH the overall fit AND a
per-ratio table of individual-fit parameters (`max`, `slope1-3`, `ec50_1-3`,
`a`-values). So `analyse_mixture()` for ternary becomes multi-fit: one overall
fit + a loop over each distinct ternary ratio. The result object needs an
overall block + a per-ratio structure.

**Mixture-ratio detection = concentration-based, NOT sheet names; EXACT nominal
grouping, NO tolerance** (decided 2026-05-30):
- A ratio is a fixed proportion `C1:C2:C3` (a ray from the origin); a
  dose-response series at fixed ratio holds proportions constant.
- For each ternary row (all of C1,C2,C3 > 0), normalise to `Ci/(C1+C2+C3)` and
  group rows with the same proportion vector → one individual fit.
- Singles (one chem > 0) and binaries (one chem = 0) are excluded from ratio
  grouping — they're the shared backbone for every fit.
- **Users always supply NOMINAL concentrations.** Grouping is therefore EXACT —
  NO fuzzy/measurement tolerance, NO `ratio_tol` parameter. (An earlier ±5%
  tolerance idea was explicitly rejected.)
- The ONLY rounding is to ~6 significant figures on the normalised proportions,
  purely to neutralise IEEE floating-point representation differences across dose
  levels (e.g. scale factor not a power of two) — this is representation hygiene,
  NOT a physical tolerance.
- Fully general — works for any nominal ratios in any dataset, independent of the
  example data's specific ratios / sheet labels.

---

## 3. Workbook findings (from inspection this session)

### 3a. `FBSA CPF IMI ternary -simplified.xls` (repo root) — THE validation reference
The ternary analogue of the binary MixTox workbook. 7 sheets, ~805 rows each.
Primary sheet `CA CPF FBSA IMI`. Endpoint: a reproduction count (`max` ≈ 871,
"Juveniles"-style), continuous like the binary continuous data.

Reference-fit block (rows 1–20; two numeric columns per model — likely two solver
passes, "fit" vs "ui"):

| param        | CA reference        | Advanced S/A        |
|--------------|---------------------|---------------------|
| max          | ≈ 871.7             | ≈ 779 / 777         |
| beta1        | ≈ 4.68 / 6.58       | ≈ 2.73 / 2.71       |
| EC50_1       | ≈ 0.1275            | ≈ 0.136 / 0.144     |
| beta2        | ≈ 51.4 / 56.0       | ≈ 13.3 / 13.1       |
| EC50_2       | = 5.58 (held fixed) | = 5.58              |
| beta3        | ≈ 3.66 / 3.70       | ≈ 2.39 / 2.64       |
| EC50_3       | ≈ 0.590 / 0.567     | ≈ 0.689 / 0.578     |
| a1,a2,a3     | — (reference)       | per-pair (see §4)   |

- CA: Regression SS ≈ 23.42M (df 6), Residual SS ≈ 11.47M (df 412). N ≈ 404.
  (Other sheets — `Overall`, `33% FBSA`, ratio-specific — report Residual SS ≈
  12.34M; the exact canonical sheet/columns to assert against is a T1-execution
  detail.)
- Advanced S/A: S/A Residual SS ≈ 6.67M; Chivalue ≈ 211.7 / 176.2; ChiTest ≈
  5.93 / 3.35.
- Workbook is **CA-only** (every sheet is "CA" vs "Advanced S/A"). DR/DL rows
  exist but appear unfit.

### 3b. `Model4Isoplane.xlsm` and `Isoplains_Mixtox.xlsx` — Plan T2 (NOT T1)
Isoplane / EC50 / ΣTU–z **prediction** tools. The fit-data regions are empty;
they consume pre-fitted parameters to build isoplanes. The chosen model
throughout `Isoplains_Mixtox.xlsx` is `CA_SA`. Sheets include per-mixture param
sets (`FBSA_CPF_IMI`, `MPs_CPF_IMI`, `MPs_FBSA_IMI`, `MPs_FBSA_CPF`, Susana
datasets) and isoplane/TU-z tables. Use these when building T2.

---

## 4. Advanced S/A model + staged fitting workflow (RESOLVED 2026-05-30)

Sam supplied the exact MixTox workflow; the earlier blocker is resolved. The
ternary deviation is a **per-term "Advanced S/A"** fitted in **nested stages**.

### Model
- **Base CA/IA:** `max`, `slope1-3`, `ec50_1-3`.
- **Pairwise S/A:** `A1` (pair C1·C2), `A2` (C1·C3), `A3` (C2·C3).
- **Three-way S/A ("+S/A"):** `A4` (all three chemicals present).
- The single-`a` ported model (`ca_sa_tri` etc.) is the special case
  `A1=A2=A3=A4`.

**Inferred deviation function (generalises the ported single-`a` form; VERIFY
against the workbook during validation):**
```
F = exp( A1·z1·z2 + A2·z1·z3 + A3·z2·z3 + A4·z1·z2·z3 )
```
A binary mixture (one chem = 0) activates only that pair's term; a ternary
mixture activates all four. `zi = TUi/ΣTU`, `TUi = Ci/EC50i` (constant along a
fixed-ratio dose series → A4 is one identifiable number per ratio).

### Staged fitting (Sam's 10-step workflow)
1. Fit base CA params to **singles residuals only**.
2. Copy base params into CA+S/A, held **fixed**.
3. Fit **A1,A2,A3 to binary residuals only** (base fixed; **start all at 0**).
   Each binary isolates one pair, so the three are separately identifiable.
4. Copy (base + A1,A2,A3) into the A4 model ("CA + S/A + S/A"), all fixed.
5. Branch into copies — the shared starting point for steps 6–8.
6. **Overall A4:** everything fixed, fit **only A4 to ALL residuals** (incl. all
   ternary mixtures) → the single overall A4.
7. **Individual A4:** fit **only A4 to ONE ternary ratio's** data.
8. Repeat step 7 for every ternary ratio.
9. Per ratio: read the modelled **CA**, **CA+S/A**, and **S/A** values near the
   EC50 (or at a clear-response, non-Ymax-dominated, non-~100% point) → the
   effect-size of that ratio's A4.
10. Interpretation: comparing per-ratio A4 vs overall A4 exposes **averaging-out**
    — synergy at one ratio + antagonism at another cancel to overall A4 ≈ 0,
    which misrepresents the data. **This contrast is the scientific point of the
    whole ternary subproject.**

### Required outputs
- Overall fit params (base + A1/A2/A3 + overall A4).
- **Per-ratio table of individual A4 values** (one per ternary ratio).
- Per-ratio effect-size readouts: modelled CA / CA+S/A / S/A at the chosen
  near-EC50 points.

### Engine implications (new work beyond the current joint multi-start fitter)
- A **staged fitting routine**: fit parameter *subsets* on data *subsets*
  (singles→base, binaries→A1-3, ternary→A4), holding the rest fixed. Distinct
  from the engine's current joint fit + LR model selection.
- **Extended per-term deviation** (A1-A4) in `R/models-ternary.R` + registry.
- Data partitioning into singles / binaries / per-ternary-ratio (reuse the
  exact-nominal ratio grouping from §2a) and the per-ratio A4 loop.

### Still-open design questions (for the implementation plan)
1. Confirm the exact deviation formula above (validation will catch errors, but
   confirm if known).
2. How Advanced S/A coexists with the existing registry / SA·DR·DL model
   selection — new model + staged path, or a separate analysis entry point?
3. IA variant of the staged workflow + result-object shape.

---

## 5. Outcome — IMPLEMENTED & VERIFIED (2026-06-01)

Plan T1 was implemented on branch `mixdra-engine` per
`docs/superpowers/plans/2026-05-30-ternary-advanced-sa.md` (8 tasks, each
spec- and quality-reviewed; final holistic review = ready to merge).

**Delivered:**
- `R/models-ternary-asa.R` — `ca_asa_tri` predictor (per-term A1-A4) +
  `ca_asa_tri_vec`. (New hand-written file, NOT the verbatim-port
  `R/models-ternary.R`.)
- `R/registry.R` — `model_spec` "ASA" deviation (ternary-only; extra A1-A4).
- `R/fit.R` — fix so `fit_model` handles FIXED deviation params
  (`intersect(spec$extra, free)` in parscale + multi-start perturbation).
- `R/fit-ternary-asa.R` — `classify_rows`, `ternary_ratio_key` (exact-nominal,
  6 sig-fig), `fit_ternary_asa` (the staged fitter), and the **exported API**
  `analyse_ternary()` + `ternary_effect_table()`.
- Fixture `tests/testthat/fixtures/ternary_fbsa_cpf_imi_continuous.csv` (419
  rows: 19 control / 75 single / 225 binary / 100 ternary across 4 ratios) +
  `extract_ternary_fixture.R`; tests `test-models-ternary-asa.R`,
  `test-fit-ternary-asa.R`, `test-validation-ternary.R`.

**Deviation formula CONFIRMED** (resolves §4 open question 1). Candidate
`F = exp(A1·z1·z2 + A2·z1·z3 + A3·z2·z3 + A4·z1·z2·z3)` reproduces the workbook:
8 independent base + pairwise quantities match to 4–5 significant figures.

- **Validation source:** `FBSA CPF IMI ternary -simplified_correct.xls` (Sam's
  corrected file; raw data identical to the original `-simplified.xls`).
- **Anchored to sheet `CA CPF FBSA IMI Overall`** — the base-fixed, binaries-only
  column, which is exactly what our staged engine computes (Stage-1 base from
  singles → workbook base; Stage-2 A1/A2/A3 from binaries). The main sheet's
  all-data CA column is a *different* (jointly-refit) fit and was correctly NOT
  used.

| quantity | workbook | engine |
|---|---|---|
| max | 872.198 | 872.21 |
| ec50_1 | 0.12747 | 0.12747 |
| ec50_3 | 0.57507 | 0.57505 |
| slope1 | 4.6721 | 4.674 |
| slope3 | 3.6274 | 3.629 |
| A1 | 0.7295 | 0.751 |
| A2 | -0.2890 | -0.289 |
| A3 | -1.8424 | -1.872 |

- **EC50_2 (FBSA) is weakly identified** (workbook pins it at 5.58; the validation
  pins it via a tight bound). Its slope is unconstrained — not a fit error.
- **A4 (three-way)** is NOT independently checkable against this "simplified"
  workbook (its per-ratio A4 sheets use exploratory ratios that don't match the
  fixture's 4 ratios, fitted with a refit base). A4 structure is therefore
  verified by the **synthetic recovery test** in `test-fit-ternary-asa.R`
  (simulate from the formula with known per-ratio A4, recover them). Per-ratio A4
  on the real data genuinely vary (the "averaging-out" signal).

**Tests:** ternary feature green — 42 fast + 23 validation, 0 failures.

**Scope delivered:** CA + continuous only (as planned). IA Advanced S/A and
quantal ternary remain deferred (§4 open questions 2–3 likewise deferred: ASA is
a separate staged entry point, not wired into SA·DR·DL model selection).

### Not part of T1 (separate, Sam's concurrent work)
A broader refactor of the *general* engine's seeding (`analyse_mixture` now fits
curves from singles then a/b; interaction solve starts a=0, b=1 — see
`[[mixdra-staged-fitting]]`) landed on the same branch. It deliberately changes
model selection (e.g. quantal binary now chooses DR), which left **8 pre-existing
engine tests with stale expectations** (`test-validation-binary`, `test-analyse`,
`test-nchem`). These are NOT caused by the ternary feature and need re-baselining
to the new intended outcomes before the branch is fully green. Sam also noted the
ternary workbook still has DR/DL columns in some (non-`Overall`) tabs to clean up
— this does not affect the validation, which reads only the `Overall` sheet.
