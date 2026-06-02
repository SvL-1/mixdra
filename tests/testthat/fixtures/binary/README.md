# Binary mixture fixtures

CSV test fixtures for the **binary** (two-chemical) mixture engine, grouped by
the **3-chemical experiment** they come from. The binary sheets are pairwise
*projections* of those ternary experiments — verified: the `(IMI)`-campaign
controls are byte-for-byte the ternary `MPs CPF IMI` controls (see the
`ternary/` fixtures and the `binary-workbook-structure` memory). Within one
experiment folder the files share the same controls + single-chemical arms, so
they are **not** statistically independent.

Columns:
- **continuous** (reproduction): `C1, C2, Res` — the two chemical doses (mg/kg)
  and the reproduction count.
- **quantal** (survival): `C1, C2, Affected, Exposed` — `Affected` is the
  **survivor** count; `pi = Affected / Exposed`.

> ⚠️ **`C1`/`C2` follow the source sheet's *column* order, not the filename.**
> The chemical in each column is in the tables below. E.g.
> `binary_ca_mps_cpf_imi` has **C1 = CPF, C2 = MPs**.

## Sources & regeneration

- Most fixtures: `binary mixtox calculations CA continuous refitted to singles and binary residuals_correct.xlsm`, kept here at `binary/` root (gitignored, local-only). It is a single workbook covering **all** experiments, so it sits at the root rather than in a per-experiment subfolder.
  Regenerate: `R -q -e 'source("tests/testthat/fixtures/binary/extract_binary_fixtures_all.R")'`
- `survival/binary_mps_cpf_quantal.csv` only: `MixTox Model_binary_MPs_CPF.xlsm` (older workbook, left at repo root since it is git-tracked).
  Regenerate: `R -q -e 'source("tests/testthat/fixtures/binary/extract_binary_fixture.R")'`

## `cpf_mps_imi/` — experiment CPF · MPs · IMI  (↔ `ternary_mps_cpf_imi`)

| Fixture | C1 | C2 | Rows | Source sheet |
|---|---|---|---|---|
| `binary_ca_mps_cpf_imi_continuous.csv` | CPF | MPs | 145 | `CA - MPs CPF (IMI)` |
| `binary_ca_mps_imi_cpf_continuous.csv` | MPs | IMI | 145 | `CA - MPs IMI (CPF)` |
| `binary_ca_imi_cpf_mps_continuous.csv` | CPF | IMI | 145 | `CA -IMI CPF (MPs)` |

## `fbsa_cpf_imi/` — experiment CPF · FBSA · IMI  (↔ `ternary_fbsa_cpf_imi`)

| Fixture | C1 | C2 | Rows | Source sheet |
|---|---|---|---|---|
| `binary_ca_cpf_fbsa_imi_continuous.csv` | CPF | FBSA | 144 | `CA - CPF FBSA (IMI)` |
| `binary_ca_fbsa_imi_cpf_continuous.csv` | FBSA | IMI | 144 | `CA - FBSA IMI (CPF)` |
| `binary_ca_cpf_imi_fbsa_continuous.csv` | CPF | IMI | 144 | `CA - CPF IMI (FBSA)` |
| `binary_ia_cpf_imi_continuous.csv` | CPF | IMI | 144 | `IA - continous data model (new)` (IA fit of the CPF×IMI pair) |

## `mps_fbsa_imi/` — experiment MPs · FBSA · IMI  (no ternary file — not sent)

| Fixture | C1 | C2 | Rows | Source sheet |
|---|---|---|---|---|
| `binary_ca_mps_fbsa_imi_continuous.csv` | MPs | FBSA | 145 | `CA - MPs FBSA (IMI)` |
| `binary_ca_mps_imi_fbsa_continuous.csv` | MPs | IMI | 144 | `CA - MPs IMI (FBSA)` |
| `binary_ca_fbsa_imi_mps_continuous.csv` | FBSA | IMI | 145 | `CA - FBSA IMI (MPs)` |

## `mps_cpf_fbsa/` — experiment MPs · CPF · FBSA  (↔ `ternary_mps_fbsa_cpf`; 3rd pair MPs×FBSA not in workbook)

| Fixture | C1 | C2 | Rows | Source sheet |
|---|---|---|---|---|
| `binary_ca_mps_cpf_fbsa_continuous.csv` | MPs | CPF | 145 | `CA - MPs CPF (FBSA)` |
| `binary_ca_cpf_fbsa_mps_continuous.csv` | FBSA | CPF | 143 | `CA - CPF FBSA (MPs)` |

## `survival/` — quantal survival endpoint (CPF×IMI), kept together

| Fixture | C1 | C2 | Rows | Source sheet |
|---|---|---|---|---|
| `binary_ca_cpf_imi_quantal.csv` | CPF | IMI | 144 | `CA - binary data model (new)` |
| `binary_ia_cpf_imi_quantal.csv` | CPF | IMI | 144 | `IA - binary data model (new)` |
| `binary_mps_cpf_quantal.csv` | MPs | IMI | 145 | older workbook, `CA - binary data model (new)` |
