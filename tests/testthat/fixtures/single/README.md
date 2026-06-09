# Single-chemical fixtures

CSV test fixtures for the **single-chemical** (one-compound) dose-response
engine and the app's "single" tab. Columns: `Conc, Res` — the dose (mg/kg) and
the reproduction count (the `upload_schema("single", "continuous")` schema).

| Fixture | Chemical | Conc range | Rows | Provenance |
|---|---|---|---|---|
| `single_cpf_continuous.csv` | CPF | 0 – 0.4 | 45 | The CPF single-chemical arm (20 controls + doses 0.025/0.05/0.1/0.2/0.4 ×5) of the **CPF·MPs·IMI** experiment — i.e. the CPF-only rows of `binary/cpf_mps_imi/` and `ternary/cpf_mps_imi/`. Exported standalone for single-tab tests; no separate source workbook or extractor. |
| `single_chlorfenapyr_continuous.csv` | chlorfenapyr | 1e-6 – 900 | 39 | Sam's chlorfenapyr trial attached to issue #6. Response drops cleanly from ~1000 to 0 (EC50 ≈ 1.1), but spans ~9 decades with a `1e-6` pseudo-control. Regression fixture for the "fitted curve drawn as a straight line on the log axis" bug: `dr_curve_data()` sampled the fitted curve on a *linear* grid over the full dose range, so only ~2 of 200 points fell in the low-dose region where the sigmoid bends. |
