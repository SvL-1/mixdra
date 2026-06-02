# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows), fit CA/IA reference + SA/DR/DL deviations via analyse_mixture(),
# compare them, and visualise the chosen (or any) fit.

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
      shiny::uiOutput(ns("errors"))
    ),

    # Stage 1 -- two single-chemical curve panels + the freeze checkpoint.
    bslib::card(
      bslib::card_header("Stage 1 · Single curves"),
      shiny::p("Review and adjust each chemical's dose-response curve, then freeze ",
               "to fit the interaction. The mixture model uses one shared ",
               shiny::tags$code("max"), " (the average of the two fits)."),
      bslib::layout_columns(
        shiny::div(shiny::h5(shiny::textOutput(ns("chem1_title"))),
                   curve_fit_ui(ns("chem1"))),
        shiny::div(shiny::h5(shiny::textOutput(ns("chem2_title"))),
                   curve_fit_ui(ns("chem2")))
      ),
      shiny::actionButton(ns("freeze"), "Freeze curves", class = "btn-primary"),
      shiny::uiOutput(ns("freeze_note"))
    ),

    # Stage 2 -- per-model interaction workspace (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Interaction model"),

        # Model picker + the all-models comparison action (alpha drives the
        # parsimony test) on one row.
        bslib::layout_columns(
          col_widths = c(5, 3, 4),
          shiny::selectInput(
            ns("model"), "Model",
            choices = c("No interaction (reference)" = "reference",
                        "Similar action (S/A)"       = "SA",
                        "Dose-ratio (DR)"            = "DR",
                        "Dose-level (DL)"            = "DL")),
          shiny::div(class = "mt-4",
                     shiny::actionButton(ns("find_best"), "Find best model")),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01)
        ),
        shiny::helpText(
          "alpha (α) is the significance threshold for the model comparison: ",
          "Find best model keeps a more complex model only if it improves the fit ",
          "at p < α (default 0.05)."),

        shiny::uiOutput(ns("interaction_help")),

        # a/b parameter grid (hidden for the reference model; b shown for DR/DL).
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
        ),
        shiny::uiOutput(ns("objective")),
        DT::DTOutput(ns("comparison")),

        shiny::tags$b("Fit options"),
        bslib::layout_columns(
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1)
        ),
        shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
        shiny::uiOutput(ns("thorough_note")),
        shiny::conditionalPanel(
          condition = "output.has_fit", ns = ns,
          bslib::card(
            bslib::card_header("Optimize all parameters (joint)"),
            shiny::p("Refines the displayed model by fitting every parameter at ",
                     "once, seeded from the current fit. Fix any parameter by ",
                     "setting its Lower = Upper."),
            shiny::fluidRow(
              shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
              shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
              shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
              shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
              shiny::column(2, shiny::tags$small(shiny::tags$b("Value (start)")))
            ),
            optimize_param_row(ns, "max", "max", "Control response (upper plateau)."),
            optimize_param_row(ns, "slope1", "slope1", "Chemical 1 curve steepness."),
            optimize_param_row(ns, "slope2", "slope2", "Chemical 2 curve steepness."),
            optimize_param_row(ns, "ec501", "EC50 1", "Chemical 1 half-effect conc."),
            optimize_param_row(ns, "ec502", "EC50 2", "Chemical 2 half-effect conc."),
            optimize_param_row(ns, "a", "a", "Overall interaction strength/direction."),
            shiny::conditionalPanel(
              condition = "input.model == 'DR' || input.model == 'DL'", ns = ns,
              optimize_param_row(ns, "b", "b",
                                 "Interaction shift with ratio / dose level.")),
            shiny::actionButton(ns("optimize_all"), "Optimize all params",
                                class = "btn-primary"),
            shiny::uiOutput(ns("optimize_readout"))
          )
        )
      )
    ),

    # Stage 3 -- diagnostics for the displayed model (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 3 · Diagnostics"),
        bslib::layout_columns(
          bslib::card(bslib::card_header("3-D response surface"),
                      plotly::plotlyOutput(ns("surface"))),
          bslib::card(bslib::card_header("2-D isobole (vs reference)"),
                      plotly::plotlyOutput(ns("isobole")))
        ),
        bslib::card(bslib::card_header("Observed vs predicted"),
                    plotly::plotlyOutput(ns("op"))),
        bslib::layout_columns(
          bslib::card(bslib::card_header("Results table"), DT::DTOutput(ns("results"))),
          bslib::card(bslib::card_header("Confidence intervals (displayed model)"),
                      DT::DTOutput(ns("cis")))
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

    parsed <- shiny::reactive({
      shiny::req(input$file); read_upload(input$file$datapath)
    })
    errs <- shiny::reactive({
      shiny::req(input$file); validate_upload(parsed(), "binary", input$response)
    })
    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    engine_df <- shiny::reactive({
      shiny::req(input$file, length(errs()) == 0)
      to_engine_df(parsed(), "binary")
    })
    m1 <- shiny::reactive(marginal_df(engine_df(), 1))
    m2 <- shiny::reactive(marginal_df(engine_df(), 2))

    # Stage 1: two embedded single-chemical fitters, one per marginal series.
    fit1 <- curve_fit_server("chem1", fit_df = m1, meta = meta, chem_field = "chem1")
    fit2 <- curve_fit_server("chem2", fit_df = m2, meta = meta, chem_field = "chem2")

    # Checkpoint state. `frozen` gates Stages 2-3 (exposed to the UI as an output).
    frozen <- shiny::reactiveVal(FALSE)
    output$frozen <- shiny::reactive(isTRUE(frozen()))
    shiny::outputOptions(output, "frozen", suspendWhenHidden = FALSE)

    # Frozen curve-parameter vector (shared max = average of the two fits).
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2())
      assemble_curve_params(fit1(), fit2())
    })

    # Per-model interaction fits, keyed by model name. Autofit/Simulate write one
    # entry; Find best writes all four. `current_fit` is whatever is stored for the
    # selected model (NULL if none yet). Keeping a store -- rather than a single
    # reactiveVal that Find best would stomp when it switches the picker -- means
    # the programmatic picker switch never re-triggers or clears a fit.
    fits_store   <- shiny::reactiveVal(list())
    last_compare <- shiny::reactiveVal(NULL)
    optimize_pre  <- shiny::reactiveVal(NULL)
    optimize_post <- shiny::reactiveVal(NULL)
    current_fit  <- shiny::reactive({
      m <- input$model
      if (is.null(m)) return(NULL)
      fits_store()[[m]]
    })

    # Reveal the Optimize-all card only once a model is displayed.
    output$has_fit <- shiny::reactive(!is.null(current_fit()))
    shiny::outputOptions(output, "has_fit", suspendWhenHidden = FALSE)

    engine_response <- shiny::reactive(
      if (input$response == "quantal") "binary" else "continuous")
    n_starts_eff <- shiny::reactive(
      if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts)

    output$freeze_note <- shiny::renderUI({
      if (is.null(fit1()) || is.null(fit2()))
        shiny::div(class = "text-muted",
                   shiny::tags$small(
                     "Fit both single curves (Autofit or Simulate) before freezing."))
    })

    # Freeze only locks the curves and reveals Stage 2 -- it does not fit.
    shiny::observeEvent(input$freeze, {
      if (is.null(fit1()) || is.null(fit2())) {
        shiny::showNotification("Fit both single curves before freezing.", type = "warning")
        return()
      }
      frozen(TRUE)
    })

    # A changed curve or header model invalidates the freeze and clears all fits.
    shiny::observeEvent(
      list(fit1(), fit2(), input$reference, input$response),
      {
        if (isTRUE(frozen())) {
          frozen(FALSE)
          fits_store(list())
          last_compare(NULL)
          optimize_pre(NULL)
          optimize_post(NULL)
          for (p in c("max", "slope1", "slope2", "ec501", "ec502", "a", "b")) {
            shiny::updateNumericInput(session, paste0("olo_", p), value = NA)
            shiny::updateNumericInput(session, paste0("ohi_", p), value = NA)
          }
        }
      },
      ignoreInit = TRUE)

    # Read the entered a/b as a named numeric (blank -> NA).
    read_ab <- function() {
      raw <- list(a = input$val_a, b = input$val_b)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

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
    })

    # Find best model: fit all four + select; store every fit, land on the chosen.
    shiny::observeEvent(input$find_best, {
      shiny::req(frozen(), curve_params())
      res <- tryCatch(
        shiny::withProgress(message = "Comparing interaction models...", value = 0.5,
          analyse_mixture(engine_df(), reference = input$reference,
                          response = engine_response(), start = curve_params(),
                          alpha = input$alpha, n_starts = n_starts_eff(),
                          time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(res)) return()
      last_compare(res)
      s <- fits_store()
      for (m in names(res$fits)) s[[m]] <- res$fits[[m]]
      fits_store(s)
      shiny::updateSelectInput(session, "model", selected = res$chosen)
      ch <- res$fits[[res$chosen]]
      shiny::updateNumericInput(session, "val_a",
        value = if ("a" %in% names(ch$par)) round(ch$par[["a"]], 4) else NA)
      shiny::updateNumericInput(session, "val_b",
        value = if ("b" %in% names(ch$par)) round(ch$par[["b"]], 4) else NA)
    })

    # Switching the picker syncs the grid to that model's stored a/b (blank if none).
    shiny::observeEvent(input$model, {
      f <- fits_store()[[input$model]]
      shiny::updateNumericInput(session, "val_a",
        value = if (!is.null(f) && "a" %in% names(f$par)) round(f$par[["a"]], 4) else NA)
      shiny::updateNumericInput(session, "val_b",
        value = if (!is.null(f) && "b" %in% names(f$par)) round(f$par[["b"]], 4) else NA)
      optimize_pre(NULL); optimize_post(NULL)
      for (p in c("max", "slope1", "slope2", "ec501", "ec502", "a", "b")) {
        shiny::updateNumericInput(session, paste0("olo_", p), value = NA)
        shiny::updateNumericInput(session, paste0("ohi_", p), value = NA)
      }
    }, ignoreInit = TRUE)

    # Keep the Optimize-all start column in sync with the displayed fit.
    shiny::observe({
      f <- current_fit()
      if (is.null(f)) return()
      for (p in names(f$par))
        shiny::updateNumericInput(session, paste0("oval_", p),
                                  value = round(f$par[[p]], 4))
    })

    # Optimize all params: jointly refine the displayed model, seeded from it.
    shiny::observeEvent(input$optimize_all, {
      shiny::req(frozen(), current_fit())
      f <- current_fit()
      spec <- model_spec(input$reference, f$deviation, 2)
      b <- collect_bounds_all(shiny::reactiveValuesToList(input), params = spec$params)
      pre <- f$objective
      newfit <- tryCatch(
        shiny::withProgress(message = "Optimizing all parameters...", value = 0.5,
          refine_joint(f, engine_df(), lower = b$lower, upper = b$upper,
                       n_starts = n_starts_eff(), time_limit = input$time_limit)),
        error = function(e) {
          shiny::showNotification(paste("Optimize failed:", conditionMessage(e)),
                                  type = "error")
          NULL
        })
      if (is.null(newfit)) return()
      s <- fits_store(); s[[f$deviation]] <- newfit; fits_store(s)
      optimize_pre(pre); optimize_post(newfit$objective)
      if ("a" %in% names(newfit$par))
        shiny::updateNumericInput(session, "val_a", value = round(newfit$par[["a"]], 4))
      if ("b" %in% names(newfit$par))
        shiny::updateNumericInput(session, "val_b", value = round(newfit$par[["b"]], 4))
    })

    # Per-model explanation: tracks the reference and the selected model.
    output$interaction_help <- shiny::renderUI({
      shiny::req(frozen())
      interaction_help(input$reference, input$model)
    })

    output$optimize_readout <- shiny::renderUI({
      shiny::req(!is.null(optimize_post()))
      lab <- if (identical(current_fit()$response, "binary")) "Deviance" else "SSR"
      improved <- optimize_post() <= optimize_pre() + 1e-9
      shiny::tags$p(
        shiny::tags$b(paste0(lab, ": ")),
        round(optimize_pre(), 2), shiny::HTML(" &rarr; "), round(optimize_post(), 2),
        if (improved) shiny::tags$span(style = "color:green", " ✓ improved"))
    })

    # Fit-objective readout for the displayed model.
    output$objective <- shiny::renderUI({
      shiny::req(current_fit())
      f <- current_fit()
      lab <- if (identical(f$response, "binary")) "Deviance" else "SSR"
      shiny::tags$p(
        shiny::tags$b(paste0(lab, ": ")), round(f$objective, 2),
        "   |   ", shiny::tags$b("n: "), f$n,
        if (isTRUE(f$simulated)) shiny::tags$em(" (simulated)"))
    })

    output$surface <- plotly::renderPlotly({
      shiny::req(frozen(), current_fit()); plot_surface(current_fit(), engine_df())
    })
    output$isobole <- plotly::renderPlotly({
      shiny::req(frozen(), current_fit())
      ref <- if (!is.null(last_compare())) last_compare()$fits$reference else NULL
      plot_isobole(current_fit(), engine_df(), reference_fit = ref)
    })
    output$op <- plotly::renderPlotly({
      shiny::req(frozen(), current_fit()); plot_obs_pred(current_fit(), engine_df())
    })

    # Table-2 style all-models matrix -- only meaningful after Find best.
    output$results <- DT::renderDT({
      shiny::req(last_compare())
      tab <- round(result_table(last_compare()), 4)
      DT::datatable(as.data.frame(tab), options = list(dom = "t"))
    })
    output$comparison <- DT::renderDT({
      shiny::req(last_compare())
      DT::datatable(last_compare()$comparison, rownames = FALSE, options = list(dom = "t"))
    })
    # CIs for a fitted model; a simulated (hand-entered) set has no CIs, so show
    # its entered values instead (honest -- those are the numbers you set).
    # For a joint fit, parameters that were pinned (Lower == Upper) were not
    # estimated, so their CI is meaningless -- blank those rows.
    output$cis <- DT::renderDT({
      shiny::req(current_fit())
      f <- current_fit()
      tab <- if (isTRUE(f$simulated)) {
        data.frame(parameter = names(f$par), value = round(unname(f$par), 4))
      } else {
        ci <- param_ci(f, engine_df(), f$reference, f$deviation, f$response)
        if (isTRUE(f$joint)) ci <- blank_pinned_ci(ci, f$fixed)
        ci
      }
      DT::datatable(tab, rownames = FALSE, options = list(dom = "t"))
    })

    list(current_fit = current_fit, last_compare = last_compare)  # return for testability
  })
}
