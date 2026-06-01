# Mixture-model unification — design & findings

**Date:** 2026-06-01
**Status:** 🚧 IN PROGRESS on branch `mixdra-unify-models`.
Step 1 DONE — a single unified predictor (`mix_response()`) is proven equivalent
to all 16 ported model functions (576/576 equivalence assertions), and the work
uncovered + fixed a real bug in the ternary IA-DL model. Step 2 (rewire
`registry.R` to dispatch through `mix_response()` and delete the 16 originals)
NOT yet done. Builds on the validated engine
(`docs/superpowers/specs/2026-05-30-mixdra-design.md`).

---

## 1. Motivation

A control-engineer's reading of the engine: despite 16 separate predictor
functions (CA/IA × reference/SA/DR/DL × binary/ternary, ~1,600 lines across
`R/models-binary.R` + `R/models-ternary.R`), the whole family is **one
parametrised static nonlinearity**. The apparent complexity is *accidental*, not
essential — it comes from the functions having been ported verbatim from the
Excel/VBA workbook (one hand-written branch per concentration subset), which was
the right call for initial validation (`tools/port_models.R`: `maxdiff == 0`) but
leaves the science buried under copy-paste.

Reading all 16, only **two things actually vary**:

1. **The interaction scalar `F`** — the single term that distinguishes the
   deviation models (Jonker et al. 2005):
   | Deviation | `F` |
   |---|---|
   | reference | `1` (CA) / shift `0` (IA) |
   | S/A | `a · ∏ zᵢ` |
   | DR | `(a + Σ bᵢ zᵢ) · ∏ zᵢ` |
   | DL | `a · (1 − b · Σ TUᵢ) · ∏ zᵢ` (CA) / `a · (1 − b · P) · ∏ zᵢ` (IA) |

2. **The combination rule** — the only thing distinguishing CA from IA:
   - **CA**: root-find `Y` such that `Σ Cᵢ/ECᵢ(Y) = F`, where
     `ECᵢ(Y) = EC50ᵢ·((max−Y)/Y)^(1/slopeᵢ)` (bisection).
   - **IA**: `base = ∏ qᵢ` (or its inclusion–exclusion complement for increasing
     endpoints), then `max · Φ(Φ⁻¹(base) ± F)`.

**`n` (binary vs ternary) is not a difference at all** — it is just how many
terms are in the sum/product. The ternary functions enumerate every subset
(`C1>0 & C2>0 & C3==0`, …) by hand; a loop over `which(concs > 0)` replaces all
of it. The degenerate cases (all-zero → control; one chemical → the plain
log-logistic curve) and the increasing/decreasing slope-sign duality are
identical across all 16.

## 2. The unified predictor

`R/mix-response.R` — `mix_response(concs, max, slopes, ec50s, reference,
deviation, a, b)` — predicts one mixture measurement point for any `n`. ~80 lines
replace the ~1,600 lines of per-case functions, and it generalises to `n ≥ 4` for
free. Equivalence is pinned by `tests/testthat/test-mixture-predict.R`: a dense
concentration grid × both slope-sign regimes × all 16 (reference × deviation ×
n) combinations, comparing `mix_response()` to each original (CA at a loose
tolerance because of bisection; IA tightly because it is closed-form).

## 3. Results

- **576 / 576 equivalence assertions pass.** The single `mix_response()`
  reproduces all 16 ported models — 15 of them to numerical tolerance with no
  changes needed.
- **A real bug was found in ternary IA-DL** (see §4) — the 16th model. The
  equivalence test flagged it; the fix makes it 16/16.
- **Regression-clean.** No failures across `test-models-binary`,
  `test-models-ternary`, `test-validation-ternary`, `test-predict-mixture`,
  `test-recovery`, `test-analyse`. The 28 ternary-validation warnings are
  pre-existing optimiser noise (identical count before and after the change).

## 4. Bug found: ternary IA-DL deviation term

While proving equivalence, `mix_response()` disagreed with `ia_dl_tri` on 32 grid
points — and **only** that function. The original applied two *different*
deviation terms across concentration subsets:

| Active subset | original term |
|---|---|
| C1+C2 | `a · (1 − b · P50) · z1·z2` ✅ |
| C1+C3 | `(a + b · P50) · z1·z3` ❌ |
| C2+C3 | `(a + b · P50) · z2·z3` ❌ |
| C1+C2+C3 | `(a + b · P50) · z1·z2·z3` ❌ |

**Jonker et al. 2005, Eq. 13** specifies the IA dose-level deviation as
`G = a · (1 − b_DL · P) · ∏ zᵢ` for the general n-chemical case. The binary
`ia_dl_bi` and the C1+C2 ternary block use this correct form; the other three
ternary subsets use a structurally different `a + b·P50` (a `+` instead of `·`,
and no `1 −`), which also loses the `b_DL` switch interpretation of Table 1
(no deviation at the EC50 level when `b_DL = 1`/`2`). It is almost certainly a
copy-paste slip when the ternary blocks were replicated from the workbook.

**Fix:** the six occurrences in `ia_dl_tri` were corrected to `a·(1 − b·P50)`
(documented inline on the function). Practical impact: any ternary fit using the
IA reference with a DL deviation was wrong whenever chemical 3 was involved
(i.e. essentially all real ternary data). CA models, all binary models, and
ternary IA reference/SA/DR were unaffected. **Sam has been notified** to check
whether any published/in-progress results used ternary IA-DL specifically.

## 5. Caveats & remaining work

- **Binary DR asymmetry (intentional, per the paper).** Jonker Eq. 8 defines
  binary DR as `(a + b₁·z1)·z1·z2` with `b₂` redundant. The unified per-chemical
  `b`-vector reproduces this via the mapping `b → c(b, 0)`; the registry/fit layer
  must apply that mapping for binary DR.
- **The ternary-only ASA model** (`R/models-ternary-asa.R`, params `A1..A4`) is a
  genuinely different functional form, not a member of this family. It stays
  separate.
- **Step 2 (not yet done):** rewire `R/registry.R` so `model_spec(...)$fn`
  dispatches through `mix_response()` (assembling vector args from the per-name
  parameters), then delete the 16 functions in `models-binary.R` /
  `models-ternary.R`. The fitter, objectives, and LR tests call only through the
  registry, so they need no changes. This is where the ~1,600 → ~120 line
  reduction actually lands; it builds on the IA-DL fix, so it is gated on Sam's
  sign-off.
