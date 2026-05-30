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

## 4. ⛔ THE BLOCKER — ternary deviation model mismatch

The workbook's ternary deviation model is **"Advanced S/A" with per-interaction
parameters `a1, a2, a3`** (and a fourth `a4` in the `Model4Isoplane` variant).
Our engine's ported ternary SA (`R/models-ternary.R::ca_sa_tri` / `ia_sa_tri`)
uses a **single shared `a`** across all interactions:
`F4 <- exp(a * z1 * z2)`, `exp(a * z1 * z3)`, `exp(a * z2 * z3)`,
`exp(a * z1 * z2 * z3)`.

**Working interpretation (UNCONFIRMED):** in "Advanced S/A" each chemical *pair*
gets its own interaction magnitude (`a1`=pair 1‑2, `a2`=1‑3, `a3`=2‑3) and `a4`
is the three‑way (all‑present) term; our model reuses one `a` everywhere.

**Consequence:** the engine's current ternary SA **cannot reproduce the
workbook's deviation fit as-is** — it would need extending to the per-pair model.
The CA **reference** fit (max, betas, EC50s) *does* match our engine and is
validatable today.

### Open questions (need the user's domain knowledge)
1. Is the per-pair `a1/a2/a3` (+`a4` three-way) interpretation correct? What is
   the exact "Advanced S/A" formula?
2. Should ternary be considered "validated" on the **CA reference alone** for
   now, with the Advanced S/A deviation implemented + validated as a later plan?
   Or is reproducing the deviation fit essential to T1?
3. Does our ported single-`a` ternary SA correspond to any model the user
   actually uses, or is it effectively superseded by Advanced S/A for 3+ chems?

### Why on hold
Per the user (2026-05-30): other things must be implemented first before these
questions can be answered. T1 is therefore parked. **Decided so far:** data
source = `FBSA CPF IMI ternary -simplified.xls`; T3 deferred to a separate step.

### To resume
1. Get the Advanced S/A formula / scope decision (§4 Q1–Q2).
2. If "CA reference only": T1 = workstreams §2 (1–4) restricted to CA reference,
   no deviation assertions. Small, unblocked.
3. If "extend engine": add a per-pair Advanced S/A ternary model to
   `R/models-ternary.R` + `R/registry.R` first, then validate.
