# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows), fit each chemical's curve, then "Fit & compare all models" fits
# CA/IA reference + SA/DR/DL via the staged method (curves fixed from the single
# compounds, only a/b fitted per model -- one model per tick for the live table),
# compares them by LR test, and selects the parsimonious winner.

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
      # Optimizer-tuning knobs apply to EVERY fit on this tab (Autofit,
      # Simulate, and the compare-all loop). They are advanced/rarely-changed,
      # so they live in a collapsed accordion at the bottom of the sidebar --
      # grouped with the other tab-wide settings, out of the Stage 1->2->3 flow.
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced fitting options",
          shiny::helpText("Apply to Autofit, Simulate, and the compare-all loop."),
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

    # Stage 2 -- the hero: fit & compare ALL interaction models (staged), live.
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Compare interaction models"),
        shiny::p("Fits every interaction model with the curves held fixed from ",
                 "the single compounds (only the interaction terms are fitted), ",
                 "reference → S/A → DR/DL, and compares them. Rows fill in as each ",
                 "model finishes; the best model is highlighted."),
        bslib::layout_columns(
          col_widths = c(7, 5),
          shiny::div(class = "mt-4",
                     shiny::actionButton(ns("compare_all"),
                                         "Fit & compare all models",
                                         class = "btn-primary")),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05,
                              min = 0, max = 1, step = 0.01)),
        shiny::helpText(
          "alpha (α) is the significance threshold for the model comparison: ",
          "a more complex model is kept only if it improves the fit at p < α ",
          "(default 0.05). The table shows every fitted model and highlights the ",
          "selected (best) one."),
        DT::DTOutput(ns("results"))
      )
    ),

    # Stage 3 -- inspect the chosen (or any) model: diagnostics, CIs, and an
    # optional by-hand explore panel.
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 3 · Inspect a model"),
        shiny::selectInput(
          ns("model"), "Model",
          choices = c("No interaction (reference)" = "reference",
                      "Similar action (S/A)"       = "SA",
                      "Dose-ratio (DR)"            = "DR",
                      "Dose-level (DL)"            = "DL")),
        shiny::uiOutput(ns("interaction_help")),
        shiny::uiOutput(ns("objective")),
        bslib::layout_columns(
          bslib::card(bslib::card_header("3-D response surface"),
                      plotly::plotlyOutput(ns("surface"), height = "520px")),
          bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                      plotly::plotlyOutput(ns("isobole"), height = "520px"))
        ),
        bslib::card(bslib::card_header("Observed vs predicted"),
                    plotly::plotlyOutput(ns("op"))),
        bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                    DT::DTOutput(ns("cis"))),

        # Explore by hand: fit/enter the displayed model's a/b at the current
        # curves (staged) -- optional, off the main compare-all path.
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel(
            "Explore by hand",
            shiny::helpText("Fit or enter this model's interaction params at the ",
                            "current curves, without re-running the comparison."),
            shiny::conditionalPanel(
              condition = "input.model != 'reference'", ns = ns,
              bslib::card(
                bslib::card_header("Parameters"),
                shiny::fluidRow(
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
                  shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
                  shiny::column(2, shiny::tags$small(shiny::tags$b("Value")))
                ),
                interaction_param_row(ns, "a", "a",
                  "Overall strength & direction (a > 0 antagonism, a < 0 synergism)."),
                shiny::conditionalPanel(
                  condition = "input.model == 'DR' || input.model == 'DL'", ns = ns,
                  interaction_param_row(ns, "b", "b",
                    "How the interaction shifts with the mixture ratio / dose level."))
              )
            ),
            shiny::div(
              shiny::actionButton(ns("autofit"), "Autofit (a, b)", class = "btn-primary"),
              shiny::actionButton(ns("simulate"), "Simulate")
            )
          )
        )
      )
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
    # example dataset (continuous IA binary mixture) so the tab is usable on
    # open. `system.file` resolves under inst/ in dev and the install tree.
    upload_path <- shiny::reactive({
      if (!is.null(input$file)) return(input$file$datapath)
      ex <- system.file("extdata", "binary_ia_cpf_imi_continuous.csv", package = "mixdra")
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
    output$frozen <- shiny::reactive(isTRUE(frozen()))
    shiny::outputOptions(output, "frozen", suspendWhenHidden = FALSE)

    # Frozen curve-parameter vector (shared max = average of the two fits).
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2())
      assemble_curve_params(fit1(), fit2())
    })

    # Per-model interaction fits, keyed by model name. Autofit/Simulate write one
    # entry; Compare all writes all four. `current_fit` is whatever is stored for
    # the selected model (NULL if none yet). When a curve changes, the stored fits
    # are re-evaluated at the new curve parameters (not cleared) so Stage 3 stays live.
    fits_store   <- shiny::reactiveVal(list())
    last_compare <- shiny::reactiveVal(NULL)

    # Number of chemicals present (binary tab -> 2; future-proofs the loop order).
    n_chem <- shiny::reactive(
      length(intersect(c("C1", "C2", "C3"), names(engine_df()))))

    # Live staged compare-all loop. `loop_queue` holds the models still to fit;
    # `stepper_on` arms the timer-driven stepper. We fit ONE model per tick and
    # let Shiny flush (paint the new table row) between ticks via invalidateLater.
    loop_queue <- shiny::reactiveVal(NULL)
    stepper_on <- shiny::reactiveVal(FALSE)

    # Joint-refined fits, keyed by model name -- a POST-SELECTION polish. The
    # comparison table NEVER reads this; only the selected model's diagnostics do.
    # This is what keeps the staged verdict (table p-values + winner) untouched.
    refined_fits <- shiny::reactiveVal(list())
    refine_pre   <- shiny::reactiveVal(NULL)   # objective before the last joint refine
    refine_post  <- shiny::reactiveVal(NULL)   # objective after

    current_fit  <- shiny::reactive({
      m <- input$model
      if (is.null(m)) return(NULL)
      fits_store()[[m]]
    })

    # The fit shown in the diagnostics / CIs for the selected model: the
    # joint-refined fit when one exists, else the staged fit. Never read by the
    # comparison table (which stays purely staged).
    display_fit <- shiny::reactive({
      m <- input$model
      if (is.null(m)) return(NULL)
      r <- refined_fits()[[m]]
      if (!is.null(r)) r else fits_store()[[m]]
    })

    engine_response <- shiny::reactive(
      if (input$response == "quantal") "binary" else "continuous")
    n_starts_eff <- shiny::reactive(
      if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts)

    output$reveal_note <- shiny::renderUI({
      if (is.null(fit1()) || is.null(fit2()))
        shiny::div(class = "text-muted",
                   shiny::tags$small(
                     "Fit both single curves (Autofit or Simulate) to reveal ",
                     "the interaction workspace."))
    })

    # A structural change (reference model or response type) makes the stored
    # fits meaningless, so clear everything.
    shiny::observeEvent(
      list(input$reference, input$response),
      {
        fits_store(list())
        last_compare(NULL)
        refined_fits(list()); refine_pre(NULL); refine_post(NULL)
      },
      ignoreInit = TRUE)

    # A changed curve does NOT blank the interaction: re-evaluate each stored fit
    # at the new curve parameters (keeping its a/b), so Stage 3 stays in sync with
    # Stage 1 live. This is what makes editing chemical 1 show up in the
    # diagnostics. The model comparison, however, was computed at the old curves,
    # so it is invalidated (re-run Compare all to compare at the new curves).
    shiny::observeEvent(
      list(fit1(), fit2()),
      {
        s <- fits_store()
        if (length(s)) {
          cp <- curve_params()
          for (m in names(s)) {
            old <- s[[m]]
            ab  <- old$par[intersect(c("a", "b"), names(old$par))]
            newf <- tryCatch(
              eval_mixture(engine_df(), old$reference, m, old$response, cp,
                           interaction = ab),
              error = function(e) NULL)
            if (is.null(newf)) next
            # carry the metadata the diagnostics / comparison / CI layer relies on
            newf$df        <- old$df
            newf$joint     <- old$joint
            newf$fixed     <- old$fixed
            newf$simulated <- old$simulated
            s[[m]] <- newf
          }
          fits_store(s)
        }
        last_compare(NULL)
        refined_fits(list()); refine_pre(NULL); refine_post(NULL)
      },
      ignoreInit = TRUE)

    # Kick: clear the store and arm the selection chain (reference -> SA -> DR -> DL).
    shiny::observeEvent(input$compare_all, {
      shiny::req(frozen(), curve_params())
      fits_store(list())
      last_compare(NULL)
      refined_fits(list()); refine_pre(NULL); refine_post(NULL)
      loop_queue(selection_chain_order(n_chem()))
      stepper_on(TRUE)
    })

    # Stepper: re-runs on each timer tick while armed. Reads the queue with
    # isolate() so only the timer (not its own writes) re-triggers it, which is
    # what forces a client paint between models.
    shiny::observe({
      if (!isTRUE(stepper_on())) return()
      shiny::invalidateLater(0)                 # schedule the next tick
      q <- shiny::isolate(loop_queue())
      if (is.null(q) || length(q) == 0) { stepper_on(FALSE); return() }
      shiny::isolate({
        dev <- q[[1]]
        s   <- fits_store()
        # Staged per-model fit: curves are FIXED at the single-compound values
        # (curve_params()) and only the interaction params (a, b) are fitted to
        # the mixture rows -- the same fit analyse_mixture() does, run one model
        # per tick for the live table. This keeps the CA/IA reference built only
        # from the single compounds, so genuine interactions surface as a/b
        # rather than being absorbed by re-fitted curves (Sam's required method;
        # do NOT switch this loop to a joint curve+interaction fit).
        fit <- tryCatch(
          fit_model(engine_df(), input$reference, dev, engine_response(),
                    start = curve_params(), fixed = names(curve_params()),
                    n_starts = n_starts_eff(), time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(
              paste0("Fit failed (", dev, "): ", conditionMessage(e)), type = "error")
            NULL
          })
        if (!is.null(fit)) { s[[dev]] <- fit; fits_store(s) }
        rest <- q[-1]
        loop_queue(rest)
        if (length(rest) == 0) {
          stepper_on(FALSE)
          cmp <- compare_fits(fits_store(), nrow(engine_df()),
                              engine_response(), input$alpha)
          last_compare(list(fits = fits_store(), comparison = cmp$comparison,
                            chosen = cmp$chosen, reference = input$reference,
                            response = engine_response()))
          shiny::updateSelectInput(session, "model", selected = cmp$chosen)
          ch <- fits_store()[[cmp$chosen]]
          shiny::updateNumericInput(session, "val_a",
            value = if ("a" %in% names(ch$par)) round(ch$par[["a"]], 4) else NA)
          shiny::updateNumericInput(session, "val_b",
            value = if ("b" %in% names(ch$par)) round(ch$par[["b"]], 4) else NA)
        }
      })
    })

    # Read the entered a/b as a named numeric (blank -> NA).
    read_ab <- function() {
      raw <- list(a = input$val_a, b = input$val_b)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

    # Optimize all params (joint): re-fit EVERY parameter of the selected model
    # at once, seeded from its staged fit (Excel-style). Stored in refined_fits
    # ONLY -- never in fits_store -- so the staged comparison verdict (table
    # p-values + highlighted winner) is never affected. Post-selection polish.
    shiny::observeEvent(input$optimize_all, {
      shiny::req(frozen(), current_fit())
      f   <- current_fit()
      pre <- f$objective
      newfit <- tryCatch(
        shiny::withProgress(message = "Optimizing all parameters...", value = 0.5,
          refine_joint(f, engine_df(),
                       n_starts = n_starts_eff(), time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Optimize failed:", conditionMessage(e)),
                                  type = "error")
          NULL
        })
      if (is.null(newfit)) return()
      r <- refined_fits(); r[[input$model]] <- newfit; refined_fits(r)
      refine_pre(pre); refine_post(newfit$objective)
    })

    # Autofit: fit only the selected model's interaction params, curves fixed.
    shiny::observeEvent(input$autofit, {
      shiny::req(frozen(), curve_params())
      dev <- input$model
      fit <- tryCatch(
        shiny::withProgress(message = "Fitting interaction...", value = 0.5,
          fit_model(engine_df(), input$reference, dev, engine_response(),
                    start = curve_params(), fixed = names(curve_params()),
                    n_starts = n_starts_eff(), time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      if ("a" %in% names(fit$par))
        shiny::updateNumericInput(session, "val_a", value = round(fit$par[["a"]], 4))
      if ("b" %in% names(fit$par))
        shiny::updateNumericInput(session, "val_b", value = round(fit$par[["b"]], 4))
      s <- fits_store(); s[[dev]] <- fit; fits_store(s)
      r <- refined_fits(); r[[dev]] <- NULL; refined_fits(r)  # staged edit supersedes a prior joint refine
      refine_pre(NULL); refine_post(NULL)
    })

    # Simulate: evaluate the selected model with the entered a/b (no refit).
    shiny::observeEvent(input$simulate, {
      shiny::req(frozen(), curve_params())
      dev <- input$model
      if (dev == "reference") {
        shiny::showNotification(
          "The reference model has no interaction parameters to simulate.", type = "message")
        return()
      }
      ab <- read_ab()
      need <- if (dev == "SA") "a" else c("a", "b")
      if (any(is.na(ab[need]))) {
        shiny::showNotification("Enter a (and b) to simulate.", type = "warning")
        return()
      }
      fit <- tryCatch(
        eval_mixture(engine_df(), input$reference, dev, engine_response(),
                     curve_params(), interaction = ab[need]),
        error = function(e) {
          shiny::showNotification(paste("Simulate failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      s <- fits_store(); s[[dev]] <- fit; fits_store(s)
      r <- refined_fits(); r[[dev]] <- NULL; refined_fits(r)  # staged edit supersedes a prior joint refine
      refine_pre(NULL); refine_post(NULL)
    })

    # Switching the picker syncs the grid to that model's stored a/b (blank if none).
    shiny::observeEvent(input$model, {
      f <- fits_store()[[input$model]]
      shiny::updateNumericInput(session, "val_a",
        value = if (!is.null(f) && "a" %in% names(f$par)) round(f$par[["a"]], 4) else NA)
      shiny::updateNumericInput(session, "val_b",
        value = if (!is.null(f) && "b" %in% names(f$par)) round(f$par[["b"]], 4) else NA)
      refine_pre(NULL); refine_post(NULL)
    }, ignoreInit = TRUE)

    # Per-model explanation: tracks the reference and the selected model.
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen())
      interaction_help(input$reference, input$model)
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

    # Before -> after objective for the last joint refine of the selected model.
    output$refine_readout <- shiny::renderUI({
      shiny::req(!is.null(refine_post()), display_fit())
      lab <- if (identical(display_fit()$response, "binary")) "Deviance" else "SSR"
      improved <- refine_post() <= refine_pre() + 1e-9
      shiny::tags$p(
        shiny::tags$b(paste0("Joint refine ", lab, ": ")),
        round(refine_pre(), 2), shiny::HTML(" &rarr; "), round(refine_post(), 2),
        if (improved) shiny::tags$span(style = "color:green", " ✓"))
    })

    # Badge: the diagnostics are showing the joint-refined fit for this model.
    output$refined_badge <- shiny::renderUI({
      m <- input$model
      if (!is.null(m) && !is.null(refined_fits()[[m]]))
        shiny::div(class = "text-info", shiny::tags$small(
          "Showing the joint-refined (all-parameters) fit for this model. ",
          "The comparison table above is unchanged (staged)."))
    })

    output$surface <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit()); plot_surface(display_fit(), engine_df())
    })
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit())
      ref <- if (!is.null(last_compare())) last_compare()$fits$reference else NULL
      plot_isobole(display_fit(), engine_df(), reference_fit = ref)
    })
    output$op <- plotly::renderPlotly({
      shiny::req(frozen(), display_fit()); plot_obs_pred(display_fit(), engine_df())
    })

    # One row per interaction model the user has fitted (model name in the first
    # column): params + objective/df + LR p-value across the columns. A single
    # Autofit/Simulate adds its own row; Compare all fills all four. The p-value
    # column and the highlighted best-model row only appear once Compare all has
    # run (it is the sole writer of `last_compare()`, which carries the parent
    # comparison + chosen model). Sorting is disabled -- row order is fixed.
    output$results <- DT::renderDT({
      fits <- fits_store()
      shiny::req(length(fits) > 0)
      ord    <- intersect(c("reference", "SA", "DR", "DL"), names(fits))
      cmp    <- if (!is.null(last_compare())) last_compare()$comparison else NULL
      chosen <- if (!is.null(last_compare())) last_compare()$chosen else NULL
      mat  <- round(result_table(list(fits = fits[ord], comparison = cmp)), 4)
      disp <- as.data.frame(t(mat), check.names = FALSE)        # models -> rows
      disp <- cbind(`Interaction model` = rownames(disp), disp, stringsAsFactors = FALSE)
      rownames(disp) <- NULL
      dt <- DT::datatable(disp, rownames = FALSE,
                          options = list(dom = "t", ordering = FALSE, scrollX = TRUE))
      if (!is.null(chosen) && chosen %in% disp[["Interaction model"]])
        dt <- DT::formatStyle(dt, "Interaction model", target = "row",
                              fontWeight = DT::styleEqual(chosen, "bold"),
                              backgroundColor = DT::styleEqual(chosen, "#d8f0d8"))
      dt
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
         refined_fits = refined_fits, display_fit = display_fit)  # return for testability
  })
}
