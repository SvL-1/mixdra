# Binary Mixture stage: upload the full binary dataset (single-chemical series +
# mixture rows), fit CA/IA reference + SA/DR/DL deviations via analyse_mixture(),
# compare them, and visualise the chosen (or any) fit.

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
      shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
      shiny::uiOutput(ns("thorough_note")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced",
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::numericInput(ns("time_limit"), "time_limit (s/model)", value = 30, min = 1)
        )
      ),
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
      shiny::actionButton(ns("freeze"), "Freeze curves → fit interactions",
                          class = "btn-primary"),
      shiny::uiOutput(ns("freeze_note"))
    ),

    # Stage 2 -- interaction-model comparison + model picker (revealed once frozen).
    shiny::conditionalPanel(
      condition = "output.frozen", ns = ns,
      bslib::card(
        bslib::card_header("Stage 2 · Interaction models"),
        shiny::selectInput(ns("model"), "Model to display", choices = NULL),
        shiny::helpText("The most parsimonious model is pre-selected; ",
                        "override above to inspect another."),
        DT::DTOutput(ns("comparison"))
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
    ns <- session$ns

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

    fit_df <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      to_engine_df(parsed(), "binary")
    })

    res_r <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      df <- to_engine_df(parsed(), "binary")
      engine_response <- if (input$response == "quantal") "binary" else "continuous"
      n_starts <- if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts
      b <- collect_bounds(shiny::reactiveValuesToList(input))
      shiny::withProgress(message = "Fitting models...", value = 0.5, {
        tryCatch(
          analyse_mixture(df, reference = input$reference, response = engine_response,
                          alpha = input$alpha, lower = b$lower, upper = b$upper,
                          n_starts = n_starts, time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
            NULL
          })
      })
    })

    # Refresh the (static) model picker after each fit, defaulting to the chosen model.
    shiny::observeEvent(res_r(), {
      shiny::updateSelectInput(session, "model",
                               choices = names(res_r()$fits), selected = res_r()$chosen)
    })

    # Fall back to the chosen model until the picker's input has populated (also
    # makes the reactive testable under shiny::testServer, where updateSelectInput
    # does not round-trip an input value).
    shown_fit <- shiny::reactive({
      shiny::req(res_r())
      m <- input$model
      if (is.null(m) || !m %in% names(res_r()$fits)) m <- res_r()$chosen
      res_r()$fits[[m]]
    })

    output$dr1 <- plotly::renderPlotly({
      plotly::layout(plot_dose_response(shown_fit(), fit_df(), chem = 1),
                     xaxis = list(title = axis_label(meta, "chem1")))
    })
    output$dr2 <- plotly::renderPlotly({
      plotly::layout(plot_dose_response(shown_fit(), fit_df(), chem = 2),
                     xaxis = list(title = axis_label(meta, "chem2")))
    })
    output$surface <- plotly::renderPlotly({
      plot_surface(shown_fit(), fit_df())
    })
    output$isobole <- plotly::renderPlotly({
      plot_isobole(shown_fit(), fit_df(), reference_fit = res_r()$fits$reference)
    })
    output$op <- plotly::renderPlotly({
      plot_obs_pred(shown_fit(), fit_df())
    })

    output$results <- DT::renderDT({
      shiny::req(res_r())
      tab <- round(result_table(res_r()), 4)
      DT::datatable(as.data.frame(tab), options = list(dom = "t"))
    })
    output$comparison <- DT::renderDT({
      shiny::req(res_r())
      DT::datatable(res_r()$comparison, rownames = FALSE, options = list(dom = "t"))
    })
    output$cis <- DT::renderDT({
      shiny::req(shown_fit())
      f <- shown_fit()
      DT::datatable(param_ci(f, fit_df(), f$reference, f$deviation, f$response),
                    rownames = FALSE, options = list(dom = "t"))
    })

    res_r   # return for testability
  })
}
