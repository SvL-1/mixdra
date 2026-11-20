# mixdra (development version)

* The **Mixture toxicity** tab recognises a (de)potentiation design on upload
  (#13): when one stressor is applied at only one to three concentrations next
  to another stressor's full dose-response curve, it says so and points to the
  **(De)potentiation** tab. With one or two levels the upload was already
  rejected for having too few points; the error now also names where the data
  do belong. With three levels the upload passes and a warning appears under
  the upload summary.

# mixdra 0.2.0

Adds a second kind of experiment: (de)potentiation, where stressor 2 is applied
at only one to three fixed levels. The app now separates it from full mixture
toxicity designs.

* New **(De)potentiation** tab: when stressor 2 is applied at only one or a few
  fixed levels, it fits stressor 1's dose-response curve per series and tests,
  per level, whether the EC50 differs from the reference series with a
  one-degree-of-freedom likelihood-ratio test (residual SS for continuous,
  deviance for quantal data). As in the Excel workbook, a series without its
  own controls borrows the reference's, which are then counted in both curves;
  the app says so next to the results. The method is generic: any stressor
  pair, any number of levels.
* Two new exported functions behind that tab: `analyse_potentiation()` runs the
  analysis (curves, tests and the conclusion text) and `plot_potentiation()`
  plots stressor 1's observed points and fitted curve for every series.
* The **Multiple stressors** tab is now called **Mixture toxicity**. Only the
  user-facing labels changed; the uploaded-file format is unchanged.
* The **Campaign** tab is now called **Multiple stressors** (#12). Only the
  user-facing labels changed; the uploaded-file format is unchanged.
* Lower/upper constraints and fixed values entered on the Singles page are now
  carried into the joint single-stressor fit and into each pair's "Optimize all
  params (joint)" (#12). Previously the panels were exploratory only and every
  constraint was silently discarded. Each parameter row gained a **Fix**
  checkbox, so pinning no longer means typing the same number into Lower and
  Upper.
* `max` is a single campaign-wide quantity (one control group, one asymptote),
  so it now has one shared constraint row above the stressor panels instead of
  one per stressor.
* `fit_single()` / `analyse_single()` gained a `fixed` argument; pinning through
  equal bounds is not possible because L-BFGS-B cannot difference inside a
  zero-width box.
* The ternary hub table reports each mixture ratio in **toxic units**
  (`TU = C / EC50`) alongside the concentration ratio it was dosed at (#11), and
  warns when a stressor's fitted EC50 lies above every dose tested for it, since
  its TU share is then an extrapolation. The TU ratio is a relabelling only: it
  groups the data identically, so no fitted value changes.

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
