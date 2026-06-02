# Single-chemical fixtures

CSV test fixtures for the **single-chemical** (one-compound) dose-response
engine and the app's "single" tab. Columns: `Conc, Res` — the dose (mg/kg) and
the reproduction count (the `upload_schema("single", "continuous")` schema).

| Fixture | Chemical | Conc range | Rows | Provenance |
|---|---|---|---|---|
| `single_cpf_continuous.csv` | CPF | 0 – 0.4 | 45 | The CPF single-chemical arm (20 controls + doses 0.025/0.05/0.1/0.2/0.4 ×5) of the **CPF·MPs·IMI** experiment — i.e. the CPF-only rows of `binary/cpf_mps_imi/` and `ternary/cpf_mps_imi/`. Exported standalone for single-tab tests; no separate source workbook or extractor. |
