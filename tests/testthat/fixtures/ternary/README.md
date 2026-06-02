# Ternary mixture fixtures

CSV test fixtures for the **ternary** (three-chemical) mixture engine, grouped
by experiment — one subfolder per 3-chemical experiment, **named to match the
`binary/` experiment folders** so the same experiment lines up across both
trees. Each subfolder holds the extracted CSV plus its source `.xls`
(gitignored, local-only — the committed artifact is the CSV).

Columns: `C1, C2, C3, Res` — the three chemical doses (mg/kg) and the
reproduction (RGR) count.

> ⚠️ **`C1`/`C2`/`C3` follow the source sheet's *column* order, not the
> filename/folder.** The chemical in each column is in the table below. E.g.
> `ternary_mps_cpf_imi` has **C1 = CPF, C2 = MPs, C3 = IMI**.

The 3 ternary experiments actually run with the 4 chemicals (CPF / MPs / IMI /
FBSA): FBSA+CPF+IMI, MPs+CPF+IMI, MPs+FBSA+CPF. The 4th triple (MPs+FBSA+IMI)
was not run. Each source workbook's primary `CA <chems>` sheet holds the full
pooled design (controls + every single arm + all binary/ternary rays); the
per-ratio tabs are subset views. Non-integer "refitted to singles" values are
kept verbatim. See the `binary-workbook-structure` memory.

Regenerate (reads the local source `.xls` in each subfolder):
`R -q -e 'source("tests/testthat/fixtures/ternary/extract_ternary_fixtures_all.R")'`

| Subfolder | Fixture | C1 | C2 | C3 | Rows | Source file / sheet |
|---|---|---|---|---|---|---|
| `fbsa_cpf_imi/` | `ternary_fbsa_cpf_imi_continuous.csv` | CPF | FBSA | IMI | 419 | `...DR FBSA CPF IMI_simplified.xls` / `CA CPF FBSA IMI` |
| `cpf_mps_imi/` | `ternary_mps_cpf_imi_continuous.csv` | CPF | MPs | IMI | 420 | `...DR MPs CPF IMI_simplified.xls` / `CA MPs CPF IMI` |
| `mps_cpf_fbsa/` | `ternary_mps_fbsa_cpf_continuous.csv` | MPs | FBSA | CPF | 418 | `...DR MPs FBSA CPF_simplified.xls` / `CA MPs FBSA CPF` |
