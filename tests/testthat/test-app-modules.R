skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

test_that("intro_server writes form inputs into the shared meta store", {
  meta <- shiny::reactiveValues()
  shiny::testServer(intro_server, args = list(meta = meta), {
    session$setInputs(expID = "E1", species = "Daphnia", endpoint = "reproduction",
                      chem1 = "CPF", chem2 = "MPs",
                      unit1 = "mg/kg", unit2 = "particles/kg")
    expect_equal(meta$expID, "E1")
    expect_equal(meta$chem1, "CPF")
    expect_equal(meta$chem2, "MPs")
    expect_equal(meta$unit1, "mg/kg")
    expect_equal(meta$unit2, "particles/kg")
    expect_equal(meta$endpoint, "reproduction")
  })
})

test_that("intro_ui builds a Shiny UI fragment", {
  ui <- intro_ui("intro")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})

test_that("axis_label uses the per-stressor name + unit, plain for the single tab", {
  meta <- list(chem1 = "CPF", unit1 = "mg/kg", unit2 = "",
               chem0 = "Heat", unit0 = "C")
  expect_equal(axis_label(meta, "chem1"), "CPF (mg/kg)")     # named + unit
  expect_equal(axis_label(meta, "chem2"), "Concentration")   # no name, no unit
  expect_equal(axis_label(meta, "chem0"), "Heat (C)")        # single tab: own name + unit
  expect_equal(axis_label(meta), "Concentration")            # no field: generic fallback
})

test_that("single_server feeds its stressor name + unit into the axis-label meta fields", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    session$setInputs(name = "Heat", unit = "C")
    expect_equal(meta$chem0, "Heat")
    expect_equal(meta$unit0, "C")
    expect_equal(axis_label(meta, "chem0"), "Heat (C)")
  })
})

test_that("single_server surfaces validation errors and withholds a fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(single_server, args = list(meta = meta), {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(Conc = c(0, 1, 2), Res = c(100, 50, 10)),
                     path, row.names = FALSE)   # only 3 distinct concentrations
    session$setInputs(response = "continuous",
                      file = list(datapath = path, name = "short.csv"))
    expect_gt(length(errs()), 0)
    expect_match(paste(errs(), collapse = " "), "distinct")
  })
})

test_that("single_ui builds a Shiny UI fragment", {
  expect_true(inherits(single_ui("single"), c("shiny.tag", "shiny.tag.list", "bslib_fragment")))
})

test_that("pair_workspace_server: Fit-interactions fills the staged block, Optimize-all appends an isolated joint block", {
  skip_on_cran()
  csv <- testthat::test_path("fixtures", "binary", "cpf_mps_imi",
                             "binary_ca_mps_cpf_imi_continuous.csv")
  skip_if_not(file.exists(csv), "binary fixture missing")
  df <- to_engine_df(read_upload(csv), "binary")

  # `base` mirrors the frozen curve-parameter vector binary_server builds from
  # the two Stage-1 curve fits (assemble_curve_params(fit1, fit2)); mutating
  # this reactiveVal mid-test stands in for a Stage-1 curve edit.
  base_rv <- shiny::reactiveVal(
    assemble_curve_params(list(par = c(max = 700, slope = 2, ec50 = 1)),
                          list(par = c(max = 600, slope = 1, ec50 = 5))))

  shiny::testServer(pair_workspace_server, args = list(
    fit_df    = shiny::reactive(df),
    base      = shiny::reactive(base_rv()),
    reference = shiny::reactive("CA"),
    response  = shiny::reactive("continuous")), {
    session$setInputs(thorough = FALSE, n_starts = 1, alpha = 0.05, time_limit = 30)
    # the comparison table renders from the very start, before any interaction
    # fit (all four rows present, empty) -- it is no longer gated on `base`.
    expect_false(is.null(output$results))

    # "Fit interaction models" runs the whole staged chain in one (blocking) go.
    session$setInputs(fit_interactions = 1)
    expect_setequal(names(fits_store()), c("reference", "SA", "DR", "DL"))
    chosen     <- last_compare()$chosen
    staged_cmp <- last_compare()$comparison

    # with no row clicked, the displayed model defaults to the winner
    expect_equal(selected_model(), chosen)
    expect_equal(current_fit()$deviation, chosen)

    # only the staged block exists so far: 4 rows, all variant "staged"
    expect_equal(length(table_rows()), 4L)
    expect_true(all(vapply(table_rows(), function(e) e$variant == "staged", logical(1))))

    # clicking a staged row (index in reference,SA,DR,DL order) selects that model
    ord <- intersect(c("reference", "SA", "DR", "DL"), names(fits_store()))
    session$setInputs(results_rows_selected = match("SA", ord))
    expect_equal(selected_model(), "SA")
    expect_equal(current_fit()$deviation, "SA")
    expect_false(isTRUE(display_fit()$joint))      # a staged row shows the staged fit

    # "Optimize all params" jointly refines EVERY model -> a separator + joint block
    staged_sa_obj <- fits_store()[["SA"]]$objective
    session$setInputs(optimize_all = 1)
    expect_setequal(names(refined_fits()), c("reference", "SA", "DR", "DL"))
    expect_true(all(vapply(refined_fits(), function(f) isTRUE(f$joint), logical(1))))
    expect_lte(refined_fits()[["SA"]]$objective, staged_sa_obj + 1e-6)

    # the table now has the staged block, a separator, and the joint block
    expect_equal(length(table_rows()), 4L + 1L + 4L)
    sep <- vapply(table_rows(), function(e) isTRUE(e$separator), logical(1))
    expect_equal(sum(sep), 1L)

    # the staged verdict is untouched (still built from fits_store only)
    expect_identical(last_compare()$chosen, chosen)
    expect_equal(last_compare()$comparison, staged_cmp)
    expect_equal(fits_store()[["SA"]]$objective, staged_sa_obj)

    # clicking the JOINT SA row shows the joint-refined fit in the diagnostics
    sa_joint <- which(vapply(table_rows(), function(e)
      identical(e$model, "SA") && identical(e$variant, "joint"), logical(1)))
    session$setInputs(results_rows_selected = sa_joint)
    expect_equal(selected_model(), "SA")
    expect_true(isTRUE(display_fit()$joint))

    # a curve edit clears BOTH the staged and joint fits -- the user re-runs the
    # buttons (curve columns still show the new values, read live in the renderer).
    base_rv(assemble_curve_params(list(par = c(max = 720, slope = 2, ec50 = 1)),
                                  list(par = c(max = 600, slope = 1, ec50 = 5))))
    session$flushReact()
    expect_null(last_compare())
    expect_equal(length(fits_store()), 0)
    expect_equal(length(refined_fits()), 0)
  })
})

test_that("a bundled example dataset ships and validates as binary continuous", {
  ex <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled example not installed")
  df <- read_upload(ex)
  expect_length(validate_upload(df, "binary", "continuous"), 0)
  expect_true(all(c("C1", "C2", "Res") %in% names(df)))
})

test_that("binary_server falls back to the bundled example before any upload", {
  skip_on_cran()
  ex <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled example not installed")
  meta <- shiny::reactiveValues(chem1 = "CPF", chem2 = "IMI")
  shiny::testServer(binary_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_length(errs(), 0)                 # no upload, yet valid
    expect_gt(nrow(engine_df()), 0)          # data is available
  })
})

test_that("binary_ui: always-visible table + Fit/Optimize buttons, then diagnostics", {
  html <- as.character(binary_ui("binary"))
  expect_false(grepl("Freeze curves", html, fixed = TRUE))
  # Stage 2 + 3 now live in the pair_workspace_ui module, nested under "work".
  expect_match(html, "binary-work-results", fixed = TRUE)
  expect_match(html, "binary-work-fit_interactions", fixed = TRUE)
  expect_match(html, "Fit interaction models", fixed = TRUE)
  expect_match(html, "binary-work-alpha", fixed = TRUE)
  expect_match(html, "threshold for the model comparison", fixed = TRUE)
  expect_match(html, "binary-work-optimize_all", fixed = TRUE)
  expect_match(html, "Optimize all params (joint)", fixed = TRUE)
  expect_match(html, "binary-work-refine_readout", fixed = TRUE)
  expect_match(html, "binary-work-refined_badge", fixed = TRUE)
  expect_match(html, "binary-work-interaction_help", fixed = TRUE)
  expect_match(html, "binary-work-objective", fixed = TRUE)
  # Stage 3: diagnostics
  expect_match(html, "binary-work-surface", fixed = TRUE)
  expect_match(html, "binary-work-isobole", fixed = TRUE)
  expect_match(html, "binary-work-op", fixed = TRUE)
  expect_match(html, "binary-work-cis", fixed = TRUE)
  expect_match(html, "binary-work-n_starts", fixed = TRUE)
  # the manual path + compare-all button + model picker are gone
  expect_false(grepl("binary-compare_all", html, fixed = TRUE))
  expect_false(grepl("binary-model", html, fixed = TRUE))
  expect_false(grepl("binary-autofit", html, fixed = TRUE))
  expect_false(grepl("binary-simulate", html, fixed = TRUE))
  expect_false(grepl("binary-val_a", html, fixed = TRUE))
  expect_false(grepl("binary-find_best", html, fixed = TRUE))
})

test_that("pair_workspace_ui builds a Shiny UI fragment", {
  ui <- pair_workspace_ui("work")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})

test_that("single_ui shows the model equation and Autofit/Simulate buttons", {
  html <- as.character(single_ui("single"))
  expect_match(html, "Y = max / (1 + (C / EC50)", fixed = TRUE)
  expect_match(html, "Autofit parameters", fixed = TRUE)
  expect_match(html, "Simulate", fixed = TRUE)
  expect_match(html, "the curve is fit to the proportion", fixed = TRUE)
})

test_that("interaction_param_row renders a value input with inert bound cells", {
  ns <- shiny::NS("binary")
  html <- as.character(interaction_param_row(ns, "a", "a", "overall strength"))
  expect_match(html, "binary-val_a", fixed = TRUE)   # editable Value input present
  expect_match(html, "—", fixed = TRUE)          # em-dash placeholder present
  expect_false(grepl("binary-lo_a", html, fixed = TRUE))  # no Lower input
  expect_false(grepl("binary-hi_a", html, fixed = TRUE))  # no Upper input
})

test_that("a bundled example dataset ships and validates as ternary continuous", {
  ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled ternary example not installed")
  df <- read_upload(ex)
  expect_length(validate_upload(df, "ternary", "continuous"), 0)
  expect_true(all(c("C1", "C2", "C3", "Res") %in% names(df)))
})

test_that("intro_server writes chem3 into the shared meta store", {
  meta <- shiny::reactiveValues()
  shiny::testServer(intro_server, args = list(meta = meta), {
    session$setInputs(chem1 = "CPF", chem2 = "FBSA", chem3 = "IMI")
    expect_equal(meta$chem3, "IMI")
  })
})

test_that("ternary_server: three curves -> Fit Advanced S/A -> hub + selection", {
  skip_on_cran()
  ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled ternary example not installed")
  meta <- shiny::reactiveValues(chem1 = "CPF", chem2 = "FBSA", chem3 = "IMI")
  shiny::testServer(ternary_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA",
                      thorough = FALSE, n_starts = 1, time_limit = 30)
    expect_length(errs(), 0)                  # bundled example is valid
    expect_gt(nrow(engine_df()), 0)

    # Simulate the three single curves (values near the validated fit).
    session$setInputs(`chem1-val_max` = 872, `chem1-val_slope` = 4.67,
                      `chem1-val_ec50` = 0.1275, `chem1-simulate` = 1)
    session$setInputs(`chem2-val_max` = 872, `chem2-val_slope` = 13,
                      `chem2-val_ec50` = 5.58, `chem2-simulate` = 1)
    session$setInputs(`chem3-val_max` = 872, `chem3-val_slope` = 3.63,
                      `chem3-val_ec50` = 0.575, `chem3-simulate` = 1)
    expect_true(isTRUE(frozen()))

    # Run the staged Advanced-S/A fit.
    session$setInputs(fit_asa = 1)
    res <- asa_res()
    expect_false(is.null(res))
    expect_setequal(names(res$pairwise), c("A1", "A2", "A3"))
    expect_gt(nrow(res$individual), 0)

    # Hub: Overall row + one row per ternary ratio; Overall selected by default.
    h <- hub_df()
    expect_equal(h$Ratio[1], "Overall")
    expect_equal(nrow(h), 1L + nrow(res$individual))
    expect_null(selected_ratio())             # row 1 (Overall) -> NULL ratio

    # Clicking a ratio row selects that ratio.
    session$setInputs(hub_rows_selected = 2)
    expect_equal(selected_ratio(), h$.key[2])

    # A curve edit clears the staged result (user re-runs Fit Advanced S/A).
    session$setInputs(`chem1-val_max` = 880, `chem1-simulate` = 2)
    expect_null(asa_res())
  })
})

test_that("ternary_ui builds the sidebar, three curve panels, and the three stages", {
  html <- as.character(ternary_ui("ternary"))
  # sidebar controls + disabled-radio script
  expect_match(html, "ternary-response", fixed = TRUE)
  expect_match(html, "ternary-reference", fixed = TRUE)
  expect_match(html, "coming soon", fixed = TRUE)
  expect_match(html, "prop('disabled', true)", fixed = TRUE)
  # Stage 1: three curve panels
  expect_match(html, "ternary-chem1-autofit", fixed = TRUE)
  expect_match(html, "ternary-chem2-autofit", fixed = TRUE)
  expect_match(html, "ternary-chem3-autofit", fixed = TRUE)
  # Stage 2: fit button + hub table
  expect_match(html, "ternary-fit_asa", fixed = TRUE)
  expect_match(html, "ternary-hub", fixed = TRUE)
  # Stage 3: plots + effect readout
  expect_match(html, "ternary-isoplane", fixed = TRUE)
  expect_match(html, "ternary-sigma_tu", fixed = TRUE)
  expect_match(html, "ternary-effect", fixed = TRUE)
  # no alpha / joint controls (deliberately dropped vs binary)
  expect_false(grepl("ternary-alpha", html, fixed = TRUE))
  expect_false(grepl("ternary-optimize_all", html, fixed = TRUE))
})

test_that("campaign_server detects stressor count from the upload", {
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_equal(store$n_chem, 3L)          # bundled example is a 3-stressor campaign
    expect_equal(store$chems, c(1L, 2L, 3L))
    expect_equal(store$reference, "CA")
  })
})

test_that("changing the campaign reference invalidates every downstream fit", {
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    v0 <- store$base_version
    store$singles <- list(`1` = list(par = c(max = 1, slope = 1, ec50 = 1)))

    session$setInputs(reference = "IA")
    expect_true(store$base_version > v0)
    expect_equal(length(store$singles), 0)   # fits cleared, not silently reused
  })
})

test_that("campaign_bump_base increments the stamp downstream stages compare against", {
  store <- shiny::reactiveValues(base_version = 0L)
  campaign_bump_base(store)
  campaign_bump_base(store)
  # isolate(): reading a reactiveValues field outside a reactive consumer
  # errors on current shiny versions.
  expect_equal(shiny::isolate(store$base_version), 2L)
})

test_that("campaign_ui builds a Shiny UI fragment", {
  ui <- campaign_ui("camp")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})

test_that("campaign_fit_base reproduces the engine's own Stage-1 joint fit", {
  df <- to_engine_df(read_upload(system.file(
    "extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")), "campaign")

  got <- campaign_fit_base(df, c(1L, 2L, 3L), "CA", "continuous")
  ref <- analyse_ternary(df, reference = "CA", response = "continuous",
                         n_starts = 1)$base

  expect_equal(sort(names(got)), sort(names(ref)))
  expect_equal(unname(got[names(ref)]), unname(ref), tolerance = 1e-4)
  expect_equal(unname(got[["max"]]), 872.2, tolerance = 1e-2)   # workbook value
})

test_that("campaign_fit_base uses binary parameter names for two stressors", {
  df <- to_engine_df(read_upload(system.file(
    "extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")), "campaign")
  two <- df[df$C3 == 0, c("C1", "C2", "Res")]

  got <- campaign_fit_base(two, c(1L, 2L), "CA", "continuous")
  expect_equal(sort(names(got)),
               sort(c("max", "slope1", "slope2", "ec501", "ec502")))
})

test_that("campaign_fit_base ignores a present-but-never-dosed third column", {
  df <- to_engine_df(read_upload(system.file(
    "extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")), "campaign")
  two <- df[df$C3 == 0, ]                    # keeps the all-zero C3 column
  expect_equal(campaign_chems(two), c(1L, 2L))

  got <- campaign_fit_base(two, campaign_chems(two), "CA", "continuous")
  expect_false(any(grepl("3", names(got))))  # no slope3 / ec50_3 leaked in
})

test_that("campaign_fit_base handles a non-contiguous stressor set", {
  df <- to_engine_df(read_upload(system.file(
    "extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")), "campaign")

  # Stressors 1 and 3 dosed, stressor 2 never dosed but its column remains.
  d13 <- df[df$C2 == 0, ]
  expect_equal(campaign_chems(d13), c(1L, 3L))

  got <- campaign_fit_base(d13, c(1L, 3L), "CA", "continuous")
  expect_equal(sort(names(got)),
               sort(c("max", "slope1", "slope2", "ec501", "ec502")))

  # Slot 2 must carry stressor 3's curve, not stressor 2's. Stressor 3's own
  # EC50 from the full campaign is ~0.575; stressor 2's is ~34.7, so a
  # mis-mapped rename would be off by nearly two orders of magnitude.
  expect_lt(unname(got[["ec502"]]), 5)
})

test_that("campaign_base reads the stored joint fit, not the panel fits", {
  store <- shiny::reactiveValues(base = c(max = 1, slope1 = 2, ec501 = 3),
                                 singles = list())
  expect_equal(unname(shiny::isolate(campaign_base(store))[["max"]]), 1)

  store$base <- NULL
  expect_null(shiny::isolate(campaign_base(store)))
})

test_that("singles_server starts with no fits and does not invent one", {
  # The store-write and version-bump behaviour is covered by the
  # campaign_bump_base unit test in Task 3; driving nested curve_fit_server
  # modules through testServer is not worth the harness cost here. The Step 6
  # manual app check is what confirms a fit reaches the store.
  store <- shiny::reactiveValues(
    n_chem = 2L, chems = c(1L, 2L), base_version = 0L, singles = list(),
    raw = data.frame(C1 = c(0, 1, 2, 3, 0, 0, 0),
                     C2 = c(0, 0, 0, 0, 1, 2, 3),
                     Res = c(100, 80, 60, 40, 90, 70, 50)))
  shiny::testServer(singles_server, args = list(store = store), {
    expect_equal(length(store$singles), 0)
    expect_equal(store$base_version, 0L)
  })
})

test_that("singles_ui builds a Shiny UI fragment", {
  ui <- singles_ui("s")
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
              inherits(ui, "bslib_fragment"))
})

test_that("pair_workspace reproduces the binary tab's staged S/A fit", {
  skip_if_not_installed("shiny")
  df <- read.csv(system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv",
                             package = "mixdra"))
  p12 <- pair_df(df, 1, 2)

  f1 <- analyse_single(single_df(df, 1))
  f2 <- analyse_single(single_df(df, 2))
  base <- assemble_curve_params(f1, f2)

  fit <- fit_model(p12, "CA", "SA", "continuous",
                   start = c(base, a = 0), fixed = names(base),
                   n_starts = 1, time_limit = 30)

  expect_true(is.finite(fit$par[["a"]]))
  expect_named(fit$par[names(base)], names(base))
  expect_equal(unname(fit$par[names(base)]), unname(base), tolerance = 1e-8)
})
