# mixdra 0.1.0

First tagged release: the point at which the engine, the app and the validation
all cover single, binary and ternary mixtures.

## The engine

* `analyse_mixture()` fits the whole reference → S/A → {DR, DL} chain, runs the
  nested likelihood-ratio tests and selects the most parsimonious model in one
  call. `analyse_single()` and `analyse_ternary()` cover the 1- and 3-chemical
  cases.
* One n-agnostic predictor (`mix_response()`) serves single, binary and ternary
  mixtures for both reference models (CA, IA) and every deviation pattern (S/A,
  DR, DL, and Advanced S/A for ternary), replacing sixteen hand-written per-case
  model functions. The legacy functions are retained as a test oracle.
* Staged maximum-likelihood fitting: curve parameters are identified from the
  single-compound rows and then held fixed while the interaction parameters are
  fitted, so an interaction cannot be absorbed by refitted curves.
* A "joint refine" mode re-fits every parameter at once, mirroring the single-SSR
  fit the Excel workbook performs, for digit-for-digit comparison.
* Continuous (sum-of-squares) and quantal (binomial likelihood) responses.
* Base-R `optim` with L-BFGS-B, a Nelder-Mead fallback for the non-smooth CA
  bisection surface, multi-start seeding, hard bounds, and per-parameter pinning.

## The app

* A `bslib` dashboard (`run_app()`) with Introduction, Single Stressor and
  Campaign tabs. A campaign is one uploaded file covering two or three
  stressors, run in order: singles → one sub-tab per pair → ternary.
* Plots are labelled with the stressor names, units and endpoint entered on the
  Introduction tab.
* An example campaign (CPF + FBSA + IMI, continuous) is loaded until you upload
  your own file, so the app is explorable without any data to hand.

## Validation

* The engine is checked against the original Excel/VBA + Solver workbooks
  (fixtures extracted from the binary and ternary workbooks) and against
  synthetic data generated with known parameters, which it must recover on
  noise-free input.
* `simulate_single()` and `simulate_mixture()` are exported, so the recovery
  check is reproducible by users rather than internal to the test suite.

## Known limitations

* The ternary stage is CA + continuous only.
* More than three chemicals, and non-monotonic (J-shaped / hormesis) curves, are
  not supported.
* The API may still change; see the pre-release note in the README.
