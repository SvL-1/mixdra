# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows) and fit each chemical's curve. "Fit interaction models" fills the
# Stage-2 comparison table -- CA/IA reference + SA/DR/DL fit via the staged method
# (curves fixed from the single compounds, only a/b per model), compared by LR
# test with the parsimonious winner highlighted. Click a row to inspect that model
# (Stage 3 diagnostics). "Optimize all params (joint)" re-fits every parameter of
# every model (refine_joint, Excel-style) and appends them as a second "Optimized
# (joint)" block of rows for inspection -- stored separately, so the staged verdict
# (p-values + winner) never changes.

#' Per-model explanation of the interaction fit (pure, for the Binary tab)
#'
#' Returns the formula and `a`/`b` meaning for one reference x deviation
#' combination. Plain HTML (no MathJax), matching the offline-safe style of
#' [single_model_equation()]. No reactivity or fit object -- unit-testable.
#'
#' The interaction factor `F` adjusts the chosen reference baseline on the
#' toxic-unit shares `z = TU / SigmaTU` (`TU = C / EC50`; `prod z` is the
#' product across chemicals).
#'
#' @param reference Reference model: `"CA"` or `"IA"`.
#' @param deviation Deviation key: `"reference"`, `"SA"`, `"DR"`, or `"DL"`.
#' @return A [shiny::tagList()] of formula + parameter meanings.
#' @keywords internal
interaction_help <- function(reference, deviation) {
  ref <- if (identical(reference, "IA")) "IA" else "CA"

  intro <- shiny::tags$p(shiny::tags$small(shiny::HTML(paste0(
    "The interaction factor F adjusts the ", ref,
    " baseline on the toxic-unit shares ",
    "(z = TU / &Sigma;TU, TU = C / EC50; &prod;z = product across chemicals)."))))

  sign_note <- shiny::tags$p(shiny::tags$small(
    shiny::tags$b("Sign of a: "),
    "a > 0 → antagonism (mixture less toxic than the reference predicts); ",
    "a < 0 → synergism (more toxic)."))

  body <- switch(
    deviation,
    reference = shiny::tags$p(shiny::HTML(paste0(
      "<b>No interaction term.</b> The mixture follows the ", ref,
      " reference exactly."))),
    SA = shiny::tagList(
      shiny::tags$p(shiny::tags$b("Similar action (S/A): "),
                    shiny::tags$code(shiny::HTML("F = a&middot;&prod;z"))),
      shiny::tags$p(shiny::tags$small(
        shiny::tags$code("a"),
        " sets the overall strength and direction of the interaction."))),
    DR = shiny::tagList(
      shiny::tags$p(shiny::tags$b("Dose-ratio dependent (DR): "),
                    shiny::tags$code(shiny::HTML(
                      "F = (a + &Sigma;b<sub>i</sub>&middot;z<sub>i</sub>)&middot;&prod;z"))),
      shiny::tags$p(shiny::tags$small(
        shiny::tags$code("a"), " = overall interaction; ", shiny::tags$code("b"),
        " = how it shifts with the mixture ratio (which chemical dominates). ",
        "For a binary mixture there is one free b (Jonker Eq. 8)."))),
    DL = {
      dl_form <- if (ref == "IA")
        "F = a&middot;(1 &minus; b&middot;P)&middot;&prod;z"
      else
        "F = a&middot;(1 &minus; b&middot;&Sigma;TU)&middot;&prod;z"
      dl_b <- if (ref == "IA")
        "how it shifts with the dose level (P = the IA-predicted effect at the mixture point)."
      else
        "how it shifts with the dose level (ΣTU = the summed toxic units at the mixture point)."
      shiny::tagList(
        shiny::tags$p(shiny::tags$b("Dose-level dependent (DL): "),
                      shiny::tags$code(shiny::HTML(dl_form))),
        shiny::tags$p(shiny::tags$small(
          shiny::tags$code("a"), " = overall interaction; ",
          shiny::tags$code("b"), " = ", dl_b)))
    },
    shiny::tags$p()  # unknown deviation -> empty
  )

  shiny::tagList(intro, body,
                 if (!identical(deviation, "reference")) sign_note)
}

#' Binary Mixture stage UI
#' @param id Module id.
#' @keywords internal
binary_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example dataset (CPF + IMI, continuous) is loaded until you upload your own.")),
      shiny::uiOutput(ns("errors")),
      # Optimizer-tuning knobs apply to EVERY fit on this tab (the staged
      # comparison loop and the joint Optimize). They are advanced/rarely-changed,
      # so they live in a collapsed accordion at the bottom of the sidebar --
      # grouped with the other tab-wide settings, out of the Stage 1->2->3 flow.
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced fitting options",
          shiny::helpText("Apply to the staged comparison loop and Optimize all params (joint)."),
          shiny::numericInput(ns("alpha"), "alpha (α)", value = 0.05,
                              min = 0, max = 1, step = 0.01),
          shiny::helpText("Significance threshold for the model comparison: a more ",
                          "complex model is kept only if it improves the fit at p < α."),
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1),
          shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
          shiny::uiOutput(ns("thorough_note"))
        )
      )
    ),

    # Stage 1 -- two single-chemical curve panels. The interaction workspace
    # (Stage 2) appears automatically once both curves are fitted.
    bslib::card(
      bslib::card_header("Stage 1 · Single curves"),
      shiny::p("Fit each chemical's dose-response curve (Autofit or Simulate). ",
               "These curves are held fixed when the interaction models are ",
               "compared below. The comparison workspace appears once both are fitted."),
      # One row per chemical, stacked vertically (sets up the ternary case --
      # chemical 3 is simply another row). Each row is settings | plots.
      shiny::div(shiny::h5(shiny::textOutput(ns("chem1_title"))),
                 curve_fit_ui(ns("chem1"))),
      shiny::div(class = "mt-4",
                 shiny::h5(shiny::textOutput(ns("chem2_title"))),
                 curve_fit_ui(ns("chem2"))),
      shiny::uiOutput(ns("reveal_note"))
    ),

    # Stage 2 -- the comparison table IS the hub and is ALWAYS visible: it shows
    # all four model rows (reference/SA/DR/DL) from the start. The curve columns
    # fill in as soon as both Stage-1 curves are fitted; "Fit interaction models"
    # fills the a/b/objective/df/p columns. Click a row to inspect it (Stage 3);
    # "Optimize all params (joint)" appends a joint-refined block of rows.
    bslib::card(
      bslib::card_header("Stage 2 · Interaction models"),
      shiny::p("Fit both single-chemical curves in Stage 1, then ",
               shiny::tags$b("Fit interaction models"), " to fill the table ",
               "(reference → S/A → DR/DL, curves held fixed from the singles). ",
               "Click a row to inspect that model; the best is highlighted."),
      shiny::div(
        shiny::actionButton(ns("fit_interactions"), "Fit interaction models",
                            class = "btn-primary")),
      shiny::helpText(
        "The highlighted row is the selected (best) model. ",
        "The \"… (joint)\" column is filled by Optimize all params and is ",
        "display-only — it never changes which model is selected."),
      DT::DTOutput(ns("results")),
      shiny::uiOutput(ns("interaction_help")),
      shiny::uiOutput(ns("refined_badge")),
      shiny::uiOutput(ns("objective")),
      shiny::div(class = "mt-2",
        shiny::actionButton(ns("optimize_all"), "Optimize all params (joint)")),
      shiny::helpText("Re-fits every parameter (curves + interaction) of each ",
                      "model at once, seeded from its staged fit — the Excel-style ",
                      "joint fit. The results are appended as a second \"Optimized ",
                      "(joint)\" block of rows; click one to inspect its surface. ",
                      "The staged block above (p-values / winner) is unchanged."),
      shiny::uiOutput(ns("refine_readout"))
    ),

    # Stage 3 -- diagnostics. ALWAYS visible: the 3-D surface shows the raw
    # observed point cloud from the moment data is loaded, and gains the fitted
    # surface once a model is selected. The isobole / obs-pred / CI panels are
    # model-based, so they fill in once a model is fitted.
    bslib::card(
      bslib::card_header("Stage 3 · Inspect the selected model"),
      bslib::layout_columns(
        bslib::card(bslib::card_header("3-D response surface"),
                    plotly::plotlyOutput(ns("surface"), height = "520px")),
        bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                    plotly::plotlyOutput(ns("isobole"), height = "520px"))
      ),
      bslib::card(bslib::card_header("Observed vs predicted"),
                  plotly::plotlyOutput(ns("op"))),
      bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                  DT::DTOutput(ns("cis")))
    )
  )
}

#' Binary Mixture stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
binary_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() paste0("binary_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("binary", input$response), file, row.names = FALSE)
    )

    output$thorough_note <- shiny::renderUI({
      if (isTRUE(input$thorough))
        shiny::div(class = "text-warning",
                   shiny::tags$small("Multi-start fitting may take several minutes."))
    })

    output$chem1_title <- shiny::renderText({
      nm <- meta$chem1; if (!is.null(nm) && nzchar(nm)) nm else "Chemical 1"
    })
    output$chem2_title <- shiny::renderText({
      nm <- meta$chem2; if (!is.null(nm) && nzchar(nm)) nm else "Chemical 2"
    })

    # Data source: the user's upload, or -- before any upload -- a bundled
    # example dataset (CPF + IMI, continuous) so the tab is usable on open.
    # `system.file` resolves under inst/ in dev and the install tree.
    upload_path <- shiny::reactive({
      if (!is.null(input$file)) return(input$file$datapath)
      ex <- system.file("extdata", "binary_ca_cpf_imi_fbsa_continuous.csv", package = "mixdra")
      if (nzchar(ex)) ex else NULL
    })
    parsed <- shiny::reactive({
      shiny::req(upload_path()); read_upload(upload_path())
    })
    errs <- shiny::reactive({
      shiny::req(upload_path()); validate_upload(parsed(), "binary", input$response)
    })
    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    engine_df <- shiny::reactive({
      shiny::req(upload_path(), length(errs()) == 0)
      to_engine_df(parsed(), "binary")
    })
    m1 <- shiny::reactive(marginal_df(engine_df(), 1))
    m2 <- shiny::reactive(marginal_df(engine_df(), 2))

    # Stage 1: two embedded single-chemical fitters, one per marginal series.
    # The chemical panels are the single source of truth for the curves; the
    # staged compare-all loop reads them (curve_params()) and never writes back.
    fit1 <- curve_fit_server("chem1", fit_df = m1, meta = meta,
                             chem_field = "chem1")
    fit2 <- curve_fit_server("chem2", fit_df = m2, meta = meta,
                             chem_field = "chem2")

    # Stages 2-3 are gated on `frozen`: both single curves fitted. There is no
    # manual freeze step -- the workspace simply appears once both fits exist.
    frozen <- shiny::reactive(!is.null(fit1()) && !is.null(fit2()))

    # Frozen curve-parameter vector (shared max = average of the two fits).
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2())
      assemble_curve_params(fit1(), fit2())
    })

    # Per-model staged interaction fits, keyed by model name. The auto-fill loop
    # writes all four (reference/SA/DR/DL). `current_fit` is whatever is stored for
    # the selected (clicked) row, NULL if none yet. A curve/reference/response
    # change re-runs the loop, rebuilding the store at the new curves.
    fits_store   <- shiny::reactiveVal(list())
    last_compare <- shiny::reactiveVal(NULL)

    # Number of chemicals present (binary tab -> 2; future-proofs the loop order).
    n_chem <- shiny::reactive(
      length(intersect(c("C1", "C2", "C3"), names(engine_df()))))

    # Joint-refined fits, keyed by model name -- a POST-SELECTION polish. These
    # appear as their OWN rows (a second "Optimized (joint)" block) below the
    # staged rows; the staged verdict (table p-values + winner) is built only
    # from fits_store and is never affected.
    refined_fits <- shiny::reactiveVal(list())

    # Staged model order (reference/SA/DR/DL for binary) and the clicked row.
    # The table shows the staged block first, then -- once Optimize-all has run --
    # a separator and the joint block. `table_model()` below holds the per-row
    # model/variant mapping the selection uses; `sel_row` persists the clicked
    # index across re-renders.
    results_order <- shiny::reactive(selection_chain_order(n_chem()))
    sel_row <- shiny::reactiveVal(integer(0))
    shiny::observeEvent(input$results_rows_selected, {
      sel_row(input$results_rows_selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    # The ordered list of table rows and their model/variant mapping. The staged
    # block is always present (one row per model in results_order()); once any
    # joint refine exists, a separator row and the joint block are appended. Each
    # entry: list(model, variant "staged"/"joint", separator TRUE/FALSE, fit).
    table_rows <- shiny::reactive({
      ord <- results_order()
      s   <- fits_store()
      rf  <- refined_fits()
      rows <- lapply(ord, function(m)
        list(model = m, variant = "staged", separator = FALSE, fit = s[[m]]))
      jmods <- ord[vapply(ord, function(m) !is.null(rf[[m]]), logical(1))]
      if (length(jmods)) {
        rows <- c(rows, list(list(model = NA_character_, variant = "joint",
                                  separator = TRUE, fit = NULL)))
        rows <- c(rows, lapply(jmods, function(m)
          list(model = m, variant = "joint", separator = FALSE, fit = rf[[m]])))
      }
      rows
    })

    # The selected row's entry (model + variant). A separator click, or no click,
    # falls back to the staged winner.
    selected_entry <- shiny::reactive({
      rows <- table_rows()
      if (length(rows) == 0) return(NULL)
      r <- sel_row()
      if (length(r) >= 1 && r[[1]] >= 1 && r[[1]] <= length(rows)) {
        e <- rows[[r[[1]]]]
        if (!isTRUE(e$separator)) return(e)
      }
      ord <- results_order()
      ch  <- if (!is.null(last_compare())) last_compare()$chosen else NULL
      m   <- if (!is.null(ch) && ch %in% ord) ch else ord[[1]]
      list(model = m, variant = "staged", separator = FALSE, fit = fits_store()[[m]])
    })

    selected_model <- shiny::reactive({
      e <- selected_entry(); if (is.null(e)) NULL else e$model
    })

    # The staged fit of the selected model (used to seed/guard Optimize-all).
    current_fit <- shiny::reactive({
      m <- selected_model()
      if (is.null(m)) return(NULL)
      fits_store()[[m]]
    })

    # The fit shown in the diagnostics / CIs: the staged or joint fit for the
    # selected row, exactly as clicked (so the user can A/B the two surfaces).
    display_fit <- shiny::reactive({
      e <- selected_entry(); if (is.null(e)) return(NULL)
      e$fit
    })

    engine_response <- shiny::reactive(
      if (input$response == "quantal") "binary" else "continuous")
    n_starts_eff <- shiny::reactive(
      if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts)

    # The reference (CA/IA) model has NO interaction parameters -- its surface is
    # fully determined by the two frozen marginal curves. So we can build it the
    # moment both curves are fitted, with no "Fit interaction models" step: it is
    # the additivity baseline. The 3-D surface shows it by default until a fitted
    # interaction model is selected; the gap between it and the points IS the
    # interaction. Cheap (no free params to estimate).
    reference_fit_live <- shiny::reactive({
      shiny::req(frozen())
      tryCatch(
        fit_model(engine_df(), input$reference, "reference", engine_response(),
                  start = curve_params(), fixed = names(curve_params()),
                  n_starts = 1, time_limit = input$time_limit),
        error = function(e) NULL)
    })

    output$reveal_note <- shiny::renderUI({
      if (is.null(fit1()) || is.null(fit2()))
        shiny::div(class = "text-muted",
                   shiny::tags$small(
                     "Fit both single curves (Autofit or Simulate) to reveal ",
                     "the interaction workspace."))
    })

    # Any change to the curves / reference / response makes the stored interaction
    # fits stale, so clear them (the table's interaction columns blank out; its
    # curve columns still show the new curve values, read live from curve_params()
    # in the renderer). The user re-runs "Fit interaction models" to refit.
    # `curve_params()` req()s both single fits, so this is inert until frozen.
    shiny::observeEvent(
      list(curve_params(), input$reference, input$response),
      {
        fits_store(list())
        last_compare(NULL)
        refined_fits(list())
      },
      ignoreInit = TRUE)

    # "Fit interaction models": fit the staged chain (reference -> SA -> DR -> DL)
    # in one go, with a progress bar. Each model is its own staged fit -- curves
    # FIXED at the single-compound values (curve_params()), only the interaction
    # params (a, b) fitted to the mixture rows -- the same fit analyse_mixture()
    # does. This keeps the CA/IA reference built only from the single compounds,
    # so genuine interactions surface as a/b rather than being absorbed by
    # re-fitted curves (Sam's required method; do NOT switch to a joint fit).
    # The whole table renders once at the end (no live streaming) -- simple and
    # robust.
    shiny::observeEvent(input$fit_interactions, {
      if (!isTRUE(frozen())) {
        shiny::showNotification(
          "Fit both single-chemical curves first (Autofit / Simulate in Stage 1).",
          type = "warning")
        return()
      }
      refined_fits(list())
      models <- selection_chain_order(n_chem())
      s <- list()
      shiny::withProgress(message = "Fitting interaction models", value = 0, {
        for (i in seq_along(models)) {
          dev <- models[[i]]
          shiny::incProgress(0, detail = paste0("model ", dev, " (", i, "/",
                                                length(models), ")"))
          fit <- tryCatch(
            fit_model(engine_df(), input$reference, dev, engine_response(),
                      start = curve_params(), fixed = names(curve_params()),
                      n_starts = n_starts_eff(), time_limit = input$time_limit),
            error = function(e) {
              shiny::showNotification(
                paste0("Fit failed (", dev, "): ", conditionMessage(e)), type = "error")
              NULL
            })
          if (!is.null(fit)) s[[dev]] <- fit
          shiny::incProgress(1 / length(models))
        }
      })
      fits_store(s)
      cmp <- compare_fits(s, nrow(engine_df()), engine_response(), input$alpha)
      last_compare(list(fits = s, comparison = cmp$comparison,
                        chosen = cmp$chosen, reference = input$reference,
                        response = engine_response()))
      sel_row(match(cmp$chosen, results_order()))   # default-select the winner row
    })

    # Optimize all params (joint): for EVERY staged model, re-fit every parameter
    # (curves + interaction) at once, seeded from its staged fit (Excel-style).
    # Stored in refined_fits ONLY -- never in fits_store -- and shown as a second
    # "Optimized (joint)" block of rows, so the staged verdict (table p-values +
    # highlighted winner) is never affected. Click a joint row to see its surface.
    shiny::observeEvent(input$optimize_all, {
      shiny::req(frozen())
      s <- fits_store()
      if (length(s) == 0) {
        shiny::showNotification(
          "Run \"Fit interaction models\" first -- the joint optimization refines those fits.",
          type = "warning")
        return()
      }
      models <- intersect(results_order(), names(s))
      r <- list()
      shiny::withProgress(message = "Optimizing all parameters (joint)", value = 0, {
        for (i in seq_along(models)) {
          m <- models[[i]]
          shiny::incProgress(0, detail = paste0("model ", m, " (", i, "/",
                                                length(models), ")"))
          newfit <- tryCatch(
            refine_joint(s[[m]], engine_df(),
                         n_starts = n_starts_eff(), time_limit = input$time_limit),
            error = function(e) {
              shiny::showNotification(
                paste0("Optimize failed (", m, "): ", conditionMessage(e)), type = "error")
              NULL
            })
          if (!is.null(newfit)) r[[m]] <- newfit
          shiny::incProgress(1 / length(models))
        }
      })
      refined_fits(r)
    })

    # Per-model explanation: tracks the reference and the selected model.
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen(), selected_model())
      interaction_help(input$reference, selected_model())
    })

    # Fit-objective readout for the displayed model.
    output$objective <- shiny::renderUI({
      shiny::req(display_fit())
      f <- display_fit()
      lab <- if (identical(f$response, "binary")) "Deviance" else "SSR"
      shiny::tags$p(
        shiny::tags$b(paste0(lab, ": ")), round(f$objective, 2),
        "   |   ", shiny::tags$b("n: "), f$n,
        if (isTRUE(f$simulated)) shiny::tags$em(" (simulated)")
        else if (isTRUE(f$joint)) shiny::tags$em(" (joint-refined)"))
    })

    # Staged -> joint objective for the selected row, when a joint row (or a model
    # with a joint refine) is selected.
    output$refine_readout <- shiny::renderUI({
      e <- selected_entry()
      shiny::req(e, identical(e$variant, "joint"))
      staged <- fits_store()[[e$model]]; joint <- refined_fits()[[e$model]]
      shiny::req(staged, joint)
      lab <- if (identical(joint$response, "binary")) "Deviance" else "SSR"
      improved <- joint$objective <= staged$objective + 1e-9
      shiny::tags$p(
        shiny::tags$b(paste0(e$model, " joint refine ", lab, " (vs staged): ")),
        round(staged$objective, 2), shiny::HTML(" &rarr; "), round(joint$objective, 2),
        if (improved) shiny::tags$span(style = "color:green", " ✓"))
    })

    # Badge: the diagnostics are showing a joint-refined (Optimized) row.
    output$refined_badge <- shiny::renderUI({
      e <- selected_entry()
      if (!is.null(e) && identical(e$variant, "joint"))
        shiny::div(class = "text-info", shiny::tags$small(
          "Showing the joint-refined (all-parameters) fit for this model. ",
          "The staged comparison block above is unchanged."))
    })

    output$surface <- plotly::renderPlotly({
      shiny::req(engine_df())
      # raw point cloud before any curves; once both curves are frozen, the
      # additivity (reference) surface; once an interaction model is selected,
      # that model's surface (display_fit()).
      f <- display_fit()
      if (is.null(f) && isTRUE(frozen())) f <- reference_fit_live()
      plot_surface(f, engine_df())
    })
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit())
      ref <- if (!is.null(last_compare())) last_compare()$fits$reference else NULL
      plot_isobole(display_fit(), engine_df(), reference_fit = ref)
    })
    output$op <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit()); plot_obs_pred(display_fit(), engine_df())
    })

    # The comparison table. The staged block (one row per model) is always shown
    # -- empty at first, curve columns filling from curve_params() once both
    # Stage-1 curves are fitted, and a/b/objective/df/p filling when "Fit
    # interaction models" completes. Once "Optimize all params" has run, a
    # separator row ("Optimized (joint)") and one joint row per model are
    # appended; each joint row carries that model's fully re-fitted parameters so
    # the user can click it and compare its surface against the staged version.
    # Only the staged block feeds the verdict (p / highlighted winner). Selection
    # is read via isolate() so clicking a row does NOT re-render the table (DT
    # highlights client-side); the current selection is reapplied on the
    # re-renders that DO happen (after a fit / optimize completes).
    output$results <- DT::renderDT({
      rows <- table_rows()
      cmp  <- if (!is.null(last_compare())) last_compare()$comparison else NULL
      chosen <- if (!is.null(last_compare())) last_compare()$chosen else NULL
      cp   <- tryCatch(curve_params(), error = function(e) NULL)  # NULL until frozen
      olab <- if (identical(engine_response(), "binary")) "Deviance" else "SSR"

      vcurve <- function(e, nm) {                   # curve param (max/slope/ec50)
        if (isTRUE(e$separator)) return(NA_real_)
        f <- e$fit
        if (!is.null(f) && nm %in% names(f$par)) round(unname(f$par[[nm]]), 4)
        else if (identical(e$variant, "staged") && !is.null(cp) && nm %in% names(cp))
          round(unname(cp[[nm]]), 4)
        else NA_real_
      }
      vip <- function(e, nm) {                      # interaction param a/b
        if (isTRUE(e$separator) || is.null(e$fit) || !(nm %in% names(e$fit$par)))
          return(NA_real_)
        round(unname(e$fit$par[[nm]]), 4)
      }
      vobj <- function(e) if (!isTRUE(e$separator) && !is.null(e$fit))
        round(e$fit$objective, 4) else NA_real_
      vdf  <- function(e) if (!isTRUE(e$separator) && !is.null(e$fit))
        as.numeric(e$fit$df) else NA_real_
      vp   <- function(e) {                         # staged verdict p only
        if (isTRUE(e$separator) || !identical(e$variant, "staged") || is.null(cmp))
          return(NA_real_)
        i <- which(cmp$model == e$model); if (length(i)) round(cmp$p[i[[1]]], 4) else NA_real_
      }
      vlabel <- function(e) if (isTRUE(e$separator)) "─ Optimized (joint) ─" else e$model
      vhl <- function(e) {
        if (isTRUE(e$separator)) return("sep")
        if (identical(e$variant, "staged") && !is.null(chosen) &&
            identical(e$model, chosen)) return("win")
        ""
      }

      disp <- data.frame(check.names = FALSE,
        `Interaction model` = vapply(rows, vlabel, character(1)),
        max    = vapply(rows, vcurve, numeric(1), "max"),
        slope1 = vapply(rows, vcurve, numeric(1), "slope1"),
        slope2 = vapply(rows, vcurve, numeric(1), "slope2"),
        ec501  = vapply(rows, vcurve, numeric(1), "ec501"),
        ec502  = vapply(rows, vcurve, numeric(1), "ec502"),
        a      = vapply(rows, vip, numeric(1), "a"),
        b      = vapply(rows, vip, numeric(1), "b"))
      disp[[olab]] <- vapply(rows, vobj, numeric(1))
      disp[["df"]] <- vapply(rows, vdf, numeric(1))
      disp[["p"]]  <- vapply(rows, vp, numeric(1))
      disp[[".hl"]] <- vapply(rows, vhl, character(1))   # hidden style flag

      hl_idx <- ncol(disp) - 1L                          # 0-based index of .hl
      dt <- DT::datatable(
        disp, rownames = FALSE,
        selection = list(mode = "single", target = "row",
                         selected = shiny::isolate(sel_row())),
        options = list(dom = "t", ordering = FALSE, scrollX = TRUE,
                       columnDefs = list(list(visible = FALSE, targets = hl_idx))))
      DT::formatStyle(dt, "Interaction model", valueColumns = ".hl", target = "row",
                      backgroundColor = DT::styleEqual(c("win", "sep"),
                                                       c("#d8f0d8", "#eeeeee")),
                      fontWeight = DT::styleEqual("win", "bold"),
                      fontStyle  = DT::styleEqual("sep", "italic"))
    })
    # CIs for a fitted model; a simulated (hand-entered) set has no CIs, so show
    # its entered values instead (honest -- those are the numbers you set).
    # For a joint fit, parameters that were pinned (Lower == Upper) were not
    # estimated, so their CI is meaningless -- blank those rows.
    output$cis <- DT::renderDT({
      shiny::req(display_fit())
      f <- display_fit()
      tab <- if (isTRUE(f$simulated)) {
        data.frame(parameter = names(f$par), value = round(unname(f$par), 4))
      } else {
        ci <- param_ci(f, engine_df(), f$reference, f$deviation, f$response)
        if (isTRUE(f$joint)) ci <- blank_pinned_ci(ci, f$fixed)
        ci
      }
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"))
    })

    list(current_fit = current_fit, last_compare = last_compare,
         refined_fits = refined_fits, display_fit = display_fit,
         selected_model = selected_model)  # return for testability
  })
}
