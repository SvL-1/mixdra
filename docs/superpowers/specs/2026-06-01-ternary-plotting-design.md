# Ternary plotting (T2) — design

**Date:** 2026-06-01
**Status:** ✅ Approved — ready for implementation plan
**Subproject:** ternary T2 (the old "Plan 3"). Builds on the validated ternary
engine (T1, `docs/superpowers/specs/2026-05-30-ternary-subproject-design.md`)
and the single+binary plotting layer
(`docs/superpowers/specs/2026-05-30-mixdra-binary-plots-design.md`).

Delivers the EC50 isoplane + ΣTU/z computation and plots for the ternary
Advanced-S/A fit, as a tested public analytical API plus thin plotly renderers,
validated against the precomputed `Isoplains_Mixtox.xlsx` reference.

---

## 1. Scope

**In scope:**
- Model-driven EC50 isoplane + ΣTU/z computation, exposed as exported, documented,
  pure functions (Approach B — the values are genuine scientific outputs, like
  Skylar's exported Excel tables, not just plot internals).
- Two plotly renderers: a 3D isoplane and a ΣTU-vs-z plot.
- Validation of the computation against `Isoplains_Mixtox.xlsx`.

**Out of scope (deferred):**
- App wiring — hooking these into the Shiny app stays in the deferred T3.
- IA reference and quantal response (the ternary engine is CA + continuous only).
- Arbitrary effect levels other than EC50 (see §2: the clean `ΣTU = F4` identity
  is EC50-specific; generalising needs the implicit solver and is a clean later
  extension).

**Input:** an `analyse_ternary()` result (`res`) plus the source data frame. `res`
already carries everything needed: `base` (max, slope1-3, ec50_1-3), `pairwise`
(A1,A2,A3), `A4_overall`, and the per-ratio `individual` A4 table.

**Audience / rendering decision:** plots are for **interactive / in-app use**
(pure plotly). We do NOT reproduce Skylar's exact publication styling — in
particular the ggplot broken-axis (`scale_y_break`) z-plots are replaced by a
plain linear plotly y-axis, since interactive zoom serves the same purpose. This
keeps the project's all-plotly convention intact.

---

## 2. The math (settled)

The ternary Advanced-S/A CA predictor (`ca_asa_tri`, `R/models-ternary-asa.R`)
solves, for response `Y`:

```
(C1/ec1) + (C2/ec2) + (C3/ec3) = F4(z),
  ec_i = EC50_i * ((Max - Y)/Y)^(1/slope_i),
  F4(z) = exp(A1*z1*z2 + A2*z1*z3 + A3*z2*z3 + A4*z1*z2*z3),
  z_i = TU_i / ΣTU,  TU_i = C_i / EC50_i.
```

**At the EC50 surface** (`Y = Max/2`), `(Max - Y)/Y = 1`, so every
`((Max - Y)/Y)^(1/slope_i) = 1` and `ec_i = EC50_i`. The equation collapses to a
closed form — no bisection needed:

```
ΣTU = F4(z)                         on the EC50 isoplane
C_i = EC50_i * z_i * F4(z)          the isoplane point for direction z
```

Consequences:
- **Pure CA** (all A = 0): `F4 ≡ 1`, so `ΣTU ≡ 1` — exactly Skylar's reference
  line, and `C_i = EC50_i * z_i`.
- **CA+S/A** ("SA"): `F4` with A1,A2,A3 and **A4 = 0**.
- **CA+S/A+S/A** ("ASA"): `F4` with A1,A2,A3 and an A4.
- **EC50 markers** for a *tested* ratio = its isoplane point computed with that
  ratio's **individual A4**. ΣTU < 1 ⇒ synergy (less mixture needed to reach
  EC50); ΣTU > 1 ⇒ antagonism. Plotting tested-ratio markers against the overall
  curve is the visual form of the per-ratio-vs-overall A4 "averaging-out" story
  that is the scientific point of the whole ternary subproject.

---

## 3. Architecture (Approach B)

Three layers, mirroring the Plan-2 plotting split (compute is first-class and
tested; renderers are thin and guarded):

```
analyse_ternary() result (res) + df
        │
        ▼
R/ternary-isoplane.R   exported, pure, plotly-free  ← validated vs Excel
  ec50_isoplane(res, model, n)
  sigma_tu_curve(res, model)
  ec50_markers(res, df)
        │
        ▼
R/plot-data-ternary.R  pure builders (shape API output for plotting)
        │
        ▼
R/plot.R (extend)      thin plotly renderers (guarded on plotly)
  plot_isoplane(res, df, ...)
  plot_sigma_tu(res, df, ...)
```

### 3.1 Computation API — `R/ternary-isoplane.R` (NEW, exported)

- **`ec50_isoplane(res, model = c("SA", "ASA"), n = 30)`**
  Returns a data frame of isoplane points over a triangular `z`-grid on the
  2-simplex (`z_i ≥ 0`, `Σz = 1`, `n` points per edge):
  `z1, z2, z3, C1, C2, C3, sigma_tu, model`. `model = "SA"` uses `A4 = 0`;
  `"ASA"` uses `res$A4_overall`. Pure function of `res`.
- **`sigma_tu_curve(res, model = c("SA", "ASA"))`**
  Returns the ΣTU-vs-z curve data per chemical: `chem, z, sigma_tu, model`.
  **The exact z-path / axis convention is pinned by the `Isoplains_Mixtox.xlsx`
  `TU-zValues` sheet** — the validation fixture (§5) maps those columns and this
  function reproduces them. Working interpretation (to confirm against the sheet):
  for chemical *x*, vary `z_x` from 0→1 while the other two are held at their
  mutual experimental ratio; at `z_x = 0` it is the binary of the other two, at
  `z_x = 1` pure *x*.
- **`ec50_markers(res, df)`**
  Returns per-tested-ratio EC50 points using each ratio's **individual A4**:
  `ratio, C1, C2, C3, sigma_tu, z1, z2, z3`. Reuses `classify_rows()` and
  `ternary_ratio_key()` from `R/fit-ternary-asa.R`.

No plotly dependency; fully unit-testable; reusable by T3 to display numeric
tables.

### 3.2 Plot builders — `R/plot-data-ternary.R` (NEW, pure)

Shape the API output for plotting, e.g. combine the SA + ASA isoplane clouds and
the EC50 markers into one frame for `plot_isoplane`; reshape `sigma_tu_curve`
into per-chemical line series plus the ΣTU = 1 reference for `plot_sigma_tu`.
Pure, unit-tested.

### 3.3 Renderers — `R/plot.R` (EDIT)

- **`plot_isoplane(res, df, ...)`** — `scatter3d` of the CA+S/A and CA+S/A+S/A
  isoplane point-clouds plus the tested-ratio EC50 markers; axes = the three
  concentrations. Mirrors Skylar's `Isoplain*` plots, sourced from our fit.
- **`plot_sigma_tu(res, df, ...)`** — ΣTU-vs-z lines per chemical (SA solid, ASA
  dashed), a horizontal ΣTU = 1 reference line, and EC50 markers. Plain plotly
  linear y-axis.
- Both first guard `requireNamespace("plotly", quietly = TRUE)`, matching the
  existing four renderers; `plotly` stays in Suggests.

---

## 4. Validation against `Isoplains_Mixtox.xlsx`

Mirrors T1's fixture→assert pattern.

- **`tests/testthat/fixtures/extract_isoplane_fixture.R`** (NEW, tracked,
  alongside `extract_ternary_fixture.R` per the existing convention) — reads the
  `FBSA_CPF_IMI`, `TU-zValues`, and `TER_EC50` sheets (filtered to the FBSA/CPF/IMI mixture we already validated in
  T1), writes small CSV fixtures into `tests/testthat/fixtures/`. This step also
  **pins the z-path / axis convention** left open in §3.1 by revealing the
  sheet's actual columns.
- **`test-validation-isoplane.R`** (NEW) — runs `analyse_ternary()` on the
  existing ternary fixture, then asserts `ec50_isoplane` / `sigma_tu_curve` /
  `ec50_markers` match the Excel values within tolerance. Because our base +
  A1-A3 already match the workbook to 4-5 sig figs (T1), any mismatch here
  isolates a computation/convention bug, not a fit difference.
- **Caveat to verify during extraction:** the Excel tool's chosen model is
  `CA_SA`; confirm its A-parameters equal ours before asserting. If the Excel
  used a jointly-refit base (as the main ternary sheet did), assert against the
  matching column — exactly as T1 anchored to the `Overall` sheet.

---

## 5. Testing strategy

- **`test-ternary-isoplane.R`** (NEW) — unit tests for the three API functions:
  analytic identities (all A = 0 ⇒ `ΣTU ≡ 1` and `C_i = EC50_i * z_i`; SA reduces
  to `F4` with A4 = 0), grid shape, marker count = number of ternary ratios.
- **`test-plot-data-ternary.R`** (NEW) — unit tests for the builders.
- **`test-plot-ternary.R`** (NEW) — smoke tests: renderers return plotly objects
  with the expected number of traces; guarded to skip if plotly is absent.
- **`test-validation-isoplane.R`** (NEW) — the §4 validation; slow (it calls
  `analyse_ternary`). Use `n_starts = 1` + the cached fixture; run in background
  and Read a marker line, per the project's R-here gotchas.

---

## 6. File structure summary

```
R/ternary-isoplane.R          NEW  exported: ec50_isoplane, sigma_tu_curve, ec50_markers
R/plot-data-ternary.R         NEW  builders (pure)
R/plot.R                      EDIT add plot_isoplane, plot_sigma_tu
tests/testthat/fixtures/extract_isoplane_fixture.R  NEW  tracked extraction tool
tests/testthat/fixtures/isoplane_*.csv      NEW
tests/testthat/test-ternary-isoplane.R      NEW
tests/testthat/test-plot-data-ternary.R     NEW
tests/testthat/test-plot-ternary.R          NEW
tests/testthat/test-validation-isoplane.R   NEW
NAMESPACE / man/              EDIT via roxygen (document())
```

---

## 7. Notes / constraints

- `plotly` and `webshot2` are installed in the user library; renderers return
  plotly objects (render headlessly via `saveWidget` → `webshot` if a PNG is
  needed).
- `R CMD check` / `rcmdcheck` cannot run here (no Rtools) — flag as a release-time
  step in the plan, not a per-task gate.
- The validation test is long; background it and Monitor for a marker, do not poll.

---

## 8. Open items for the implementation plan

1. Confirm the `TU-zValues` / `TER_EC50` sheet columns and the exact z-path
   definition during fixture extraction (Task 1), then lock `sigma_tu_curve` to
   it.
2. Confirm the Excel isoplane's A-parameters match our fit before choosing the
   assertion columns (§4 caveat).
3. Triangular `z`-grid resolution default (`n = 30`) — adjust if the isoplane
   cloud is too sparse/dense for a usable 3D plot.
