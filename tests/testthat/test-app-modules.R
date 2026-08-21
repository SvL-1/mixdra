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

  # `base` mirrors the frozen curve-parameter vector a campaign pair builds
  # from the two Stage-1 curve fits (assemble_curve_params(fit1, fit2));
  # mutating this reactiveVal mid-test stands in for a Stage-1 curve edit.
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

test_that("ternary_server: consumes the store's frozen base + pairwise -> hub + selection", {
  skip_on_cran()
  ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
  skip_if_not(nzchar(ex) && file.exists(ex), "bundled ternary example not installed")
  df <- to_engine_df(read_upload(ex), "campaign")
  base <- campaign_fit_base(df, c(1L, 2L, 3L), "CA", "continuous")

  # Real per-pair S/A `a` values, fit exactly as pair_workspace_server does
  # (curves fixed at the pair's slice of `base`) -- this stands in for what
  # its on_fit callback would have published into store$pairs.
  pair_a <- function(i, j) {
    b <- pair_base(base, i, j)
    fit_model(pair_df(df, i, j), "CA", "SA", "continuous",
             start = b, fixed = names(b), n_starts = 1, time_limit = 30)$par[["a"]]
  }
  a12 <- pair_a(1, 2); a13 <- pair_a(1, 3); a23 <- pair_a(2, 3)

  store <- shiny::reactiveValues(
    raw = df, chems = c(1L, 2L, 3L), reference = "CA", response = "continuous",
    base = base,
    pairs = list(`12` = list(a = a12, chosen = "SA"),
                `13` = list(a = a13, chosen = "DR"),   # exercises the sa_note
                `23` = list(a = a23, chosen = "SA")))

  shiny::testServer(ternary_server, args = list(store = store), {
    session$setInputs(n_starts = 1)

    res <- asa_res()
    expect_false(is.null(res))
    expect_setequal(names(res$pairwise), c("A1", "A2", "A3"))
    # the frozen pairwise values are used AS-IS, not re-derived (A1 = pair
    # 1-2, A2 = pair 1-3, A3 = pair 2-3 -- ca_asa_tri()'s term order).
    expect_equal(unname(res$pairwise[["A1"]]), unname(a12), tolerance = 1e-8)
    expect_equal(unname(res$pairwise[["A2"]]), unname(a13), tolerance = 1e-8)
    expect_equal(unname(res$pairwise[["A3"]]), unname(a23), tolerance = 1e-8)
    expect_gt(nrow(res$individual), 0)

    # Hub: Overall row + one row per ternary ratio; Overall selected by default.
    h <- hub_df()
    expect_equal(h$Ratio[1], "Overall")
    expect_equal(nrow(h), 1L + nrow(res$individual))
    expect_null(selected_ratio())             # row 1 (Overall) -> NULL ratio

    # Clicking a ratio row selects that ratio.
    session$setInputs(hub_rows_selected = 2)
    expect_equal(selected_ratio(), h$.key[2])

    # The sa_note flags pair 1-3 (chosen = "DR", not S/A) and only that pair.
    # renderUI's output slot is list(html, deps); html is what the browser sees.
    note <- as.character(output$sa_note$html)
    expect_match(note, "13", fixed = TRUE)
    expect_false(grepl("12", note, fixed = TRUE))
    expect_false(grepl("23", note, fixed = TRUE))
  })
})

test_that("ternary_server: no fit until store$raw/base/pairwise are all present", {
  store <- shiny::reactiveValues(raw = NULL, chems = NULL, reference = NULL,
                                 response = NULL, base = NULL, pairs = list())
  shiny::testServer(ternary_server, args = list(store = store), {
    # asa_res() is a plain reactive gated by shiny::req(); with nothing in the
    # store yet it fails req() (a silent-error condition), not a plain NULL.
    expect_error(asa_res(), class = "shiny.silent.error")
  })
})

test_that("ternary_ui has no upload/curve-fit UI of its own and shows the three stages", {
  html <- as.character(ternary_ui("ternary"))
  # the sa_note output belongs to ternary_server, not the pair tab
  expect_match(html, "ternary-sa_note", fixed = TRUE)
  # the hub table + inspect plots remain
  expect_match(html, "ternary-hub", fixed = TRUE)
  expect_match(html, "ternary-isoplane", fixed = TRUE)
  expect_match(html, "ternary-sigma_tu", fixed = TRUE)
  expect_match(html, "ternary-effect", fixed = TRUE)
  expect_match(html, "ternary-n_starts", fixed = TRUE)
  # no local upload, response/reference toggle, or Stage-1 curve panels --
  # the campaign store provides all of that now
  expect_false(grepl("ternary-file", html, fixed = TRUE))
  expect_false(grepl("ternary-template", html, fixed = TRUE))
  expect_false(grepl("ternary-response", html, fixed = TRUE))
  expect_false(grepl("ternary-reference", html, fixed = TRUE))
  expect_false(grepl("ternary-chem1-autofit", html, fixed = TRUE))
  expect_false(grepl("ternary-chem2-autofit", html, fixed = TRUE))
  expect_false(grepl("ternary-chem3-autofit", html, fixed = TRUE))
  expect_false(grepl("ternary-fit_asa", html, fixed = TRUE))
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

test_that("campaign wiring carries upload, base and settings into a pair workspace", {
  skip_if_not_installed("shiny")
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")

    # The upload reaches the store, and the pair slice is derived from it.
    expect_gt(nrow(store$raw), 0)
    expect_equal(store$n_chem, 3L)
    p12 <- pair_df(store$raw, 1, 2)
    expect_true(all(c("C1", "C2") %in% names(p12)))
    expect_true(any(p12$C1 > 0 & p12$C2 > 0))

    # Settings propagate as the workspace will read them.
    expect_equal(store$reference, "CA")
    expect_equal(store$response, "continuous")

    # Base starts absent, and the pair tabs are gated on it.
    expect_null(shiny::isolate(campaign_base(store)))
    store$base <- c(max = 800, slope1 = 2, slope2 = 3, ec501 = 0.1, ec502 = 0.5)
    expect_false(is.null(shiny::isolate(campaign_base(store))))
  })
})

test_that("campaign wiring: each pair workspace fits its OWN pair's data with the correctly-sliced base", {
  # This is the permanent end-to-end coverage for the campaign_server ->
  # pair_workspace_server delegation (Task 5's lost glue, restored here per
  # Task 6). It must fail if EITHER the registration loop wires the wrong
  # pair's data/settings into a slot (the local()-capture bug), OR the base
  # reaching a pair is not correctly reduced to that pair's binary parameter
  # names (the ec50_1..3/ec501-502 naming-mismatch defect found in review:
  # without pair_base(), every 3-stressor pair's EC50s silently collapsed to
  # ~0 and the second slot's slope/EC50 came from the wrong stressor).
  skip_if_not_installed("shiny")
  skip_on_cran()
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_equal(store$chems, c(1L, 2L, 3L))   # bundled example: 3 stressors

    # Freeze the singles, as the Singles page would, so the pair tabs unlock.
    store$base <- campaign_fit_base(store$raw, store$chems, store$reference,
                                    store$response)
    session$flushReact()

    # --- Pair 1-3 -------------------------------------------------------
    b13 <- pair_base(store$base, 1, 3)
    session$setInputs(`pair13-thorough` = FALSE, `pair13-n_starts` = 1,
                      `pair13-alpha` = 0.05, `pair13-time_limit` = 30)
    session$setInputs(`pair13-fit_interactions` = 1)

    # Property 2 (base + settings reach the child): the fitted SA model for
    # pair 1-3, run through the campaign, holds the curve params fixed at
    # pair_base(store$base, 1, 3) exactly -- stressor 1's EC50 in slot 1,
    # stressor 3's in slot 2, NOT zeroed and NOT stressor 2's. This is the
    # single assertion that would have caught the pre-fix defect.
    sa13 <- pair_workspaces[["13"]]$last_compare()$fits$SA
    expect_equal(unname(sa13$par[["ec501"]]), unname(b13[["ec501"]]), tolerance = 1e-6)
    expect_equal(unname(sa13$par[["ec502"]]), unname(b13[["ec502"]]), tolerance = 1e-6)
    expect_equal(unname(sa13$par[["ec502"]]), unname(store$base[["ec50_3"]]), tolerance = 1e-6)
    expect_false(isTRUE(all.equal(unname(sa13$par[["ec502"]]),
                                  unname(store$base[["ec50_2"]]))))  # not stressor 2's

    # And it matches a fit computed directly on the pair's own slice with the
    # same base and settings -- reference/response/base all propagate, not
    # just "some" fit ran.
    p13 <- pair_df(store$raw, 1, 3)
    direct13 <- fit_model(p13, "CA", "SA", "continuous",
                          start = b13, fixed = names(b13),
                          n_starts = 1, time_limit = 30)
    expect_equal(unname(sa13$par), unname(direct13$par), tolerance = 1e-4)
    expect_equal(sa13$objective, direct13$objective, tolerance = 1e-4)

    # --- Pair 1-2 (property 1: its OWN data, not pair 1-3's) ------------
    b12 <- pair_base(store$base, 1, 2)
    session$setInputs(`pair12-thorough` = FALSE, `pair12-n_starts` = 1,
                      `pair12-alpha` = 0.05, `pair12-time_limit` = 30)
    session$setInputs(`pair12-fit_interactions` = 1)
    sa12 <- pair_workspaces[["12"]]$last_compare()$fits$SA
    expect_equal(unname(sa12$par[["ec501"]]), unname(b12[["ec501"]]), tolerance = 1e-6)
    expect_equal(unname(sa12$par[["ec502"]]), unname(b12[["ec502"]]), tolerance = 1e-6)

    # A closure bug that let every slot capture the last-iterated pair (2-3)
    # would make pair 1-2's frozen ec502 equal pair 1-3's; they must differ,
    # since slot 2 is a different stressor (2 vs 3) with a very different EC50.
    expect_false(isTRUE(all.equal(unname(sa12$par[["ec502"]]),
                                  unname(sa13$par[["ec502"]]))))
    expect_equal(unname(sa12$par[["ec502"]]), unname(store$base[["ec50_2"]]),
                tolerance = 1e-6)
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

test_that("campaign_pairs enumerates the pairs in A-B, A-C, B-C order", {
  expect_equal(campaign_pairs(c(1L, 2L, 3L)),
               list(c(1L, 2L), c(1L, 3L), c(2L, 3L)))
  expect_equal(campaign_pairs(c(1L, 2L)), list(c(1L, 2L)))
  expect_equal(pair_key(2, 3), "23")
})

test_that("campaign_pairwise maps pair S/A values onto A1/A2/A3", {
  store <- shiny::reactiveValues(
    chems = c(1L, 2L, 3L),
    pairs = list(`12` = list(a = 0.5), `13` = list(a = -0.2)))
  # bare reactiveValues reads error outside a reactive consumer in current
  # shiny -- isolate() at the test call site (not inside campaign_pairwise
  # itself, which must stay reactive for asa_res() to invalidate correctly).
  expect_null(shiny::isolate(campaign_pairwise(store)))

  shiny::isolate(store$pairs$`23` <- list(a = 1.1))
  pw <- shiny::isolate(campaign_pairwise(store))
  expect_equal(unname(pw[["A1"]]), 0.5)    # pair 1-2
  expect_equal(unname(pw[["A2"]]), -0.2)   # pair 1-3
  expect_equal(unname(pw[["A3"]]), 1.1)    # pair 2-3
})

test_that("campaign_pairwise goes stale (NULL) when a pair's base_version no longer matches the campaign's", {
  # Direct-store unit test: all three pairs present and internally consistent
  # (same base_version as the campaign) -> non-NULL. Bump ONLY the campaign's
  # stamp (as a single-stressor refit would, via campaign_bump_base()) without
  # touching the pairs -- exactly the silent-staleness scenario: the pairs'
  # `a` values are unchanged but were fitted against a now-superseded base.
  store <- shiny::reactiveValues(
    base_version = 2L,
    chems = c(1L, 2L, 3L),
    pairs = list(`12` = list(a = 0.5, chosen = "SA", base_version = 2L),
                `13` = list(a = -0.2, chosen = "SA", base_version = 2L),
                `23` = list(a = 1.1, chosen = "SA", base_version = 2L)))
  expect_false(is.null(shiny::isolate(campaign_pairwise(store))))

  campaign_bump_base(store)   # simulates a single-stressor curve being refit
  expect_equal(shiny::isolate(store$base_version), 3L)
  expect_null(shiny::isolate(campaign_pairwise(store)))   # stale -- must not mix generations
})

test_that("campaign_server: refitting the base after all three pairs are fit makes campaign_pairwise NULL again (no silent generation mixing)", {
  # End-to-end version of the guard above, through the real production wiring:
  # Singles' "Fit single-stressor curves" button bumps store$base_version but
  # does NOT clear store$pairs (only a file/response/reference change does) --
  # so without the base_version guard in campaign_pairwise(), the ternary
  # stage would silently combine the fresh base with stale pairwise terms.
  skip_if_not_installed("shiny")
  skip_on_cran()
  meta <- shiny::reactiveValues()
  shiny::testServer(campaign_server, args = list(meta = meta), {
    session$setInputs(response = "continuous", reference = "CA")
    expect_equal(store$chems, c(1L, 2L, 3L))

    session$setInputs(`singles-fit_singles` = 1)
    session$flushReact()
    expect_false(is.null(campaign_base(store)))
    v0 <- store$base_version

    # Fit all three pairs for real; each on_fit publishes its S/A `a` AND the
    # base_version it was fitted against (v0).
    session$setInputs(`pair12-n_starts` = 1, `pair12-alpha` = 0.05, `pair12-time_limit` = 30,
                      `pair13-n_starts` = 1, `pair13-alpha` = 0.05, `pair13-time_limit` = 30,
                      `pair23-n_starts` = 1, `pair23-alpha` = 0.05, `pair23-time_limit` = 30)
    session$setInputs(`pair12-fit_interactions` = 1)
    session$setInputs(`pair13-fit_interactions` = 1)
    session$setInputs(`pair23-fit_interactions` = 1)

    expect_equal(store$pairs[["12"]]$base_version, v0)
    expect_equal(store$pairs[["13"]]$base_version, v0)
    expect_equal(store$pairs[["23"]]$base_version, v0)
    expect_false(is.null(campaign_pairwise(store)))    # all three present + fresh -> usable

    # Re-fit the singles: the base changes and base_version bumps, but
    # store$pairs is left untouched (still stamped at v0) -- exactly the
    # scenario the guard exists for.
    session$setInputs(`singles-fit_singles` = 2)
    session$flushReact()
    expect_gt(store$base_version, v0)
    expect_equal(store$pairs[["12"]]$base_version, v0)   # still the OLD stamp
    expect_null(campaign_pairwise(store))                # now correctly gated
  })
})

test_that("pair_has_rows spots a pair with no mixture rows", {
  df <- campaign_fixture()
  expect_true(pair_has_rows(df, 1, 2))
  expect_true(pair_has_rows(df, 2, 3))
  df2 <- df[!(df$C2 > 0 & df$C3 > 0), ]
  expect_false(pair_has_rows(df2, 2, 3))
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
