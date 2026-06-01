# Single Chemical stage: upload one chemical's dose-response data, fit a
# three-parameter log-logistic curve via analyse_single(), and show the curve,
# observed-vs-predicted, and a parameter table.

#' Concentration-axis label from the shared meta store
#'
#' With `chem_field` (e.g. "chem1") uses that chemical's name; without it uses a
#' generic "Concentration". Appends the unit when present. Used by the single
#' stage (generic) and the binary stage (per chemical).
#' @keywords internal
axis_label <- function(meta, chem_field = NULL) {
  nm <- if (!is.null(chem_field)) meta[[chem_field]] else NULL
  unit <- meta$unit
  base <- if (!is.null(nm) && nzchar(nm)) nm else "Concentration"
  if (!is.null(unit) && nzchar(unit)) paste0(base, " (", unit, ")") else base
}

#' Static "About this model" help markup for the Single Chemical tab
#'
#' Explains the fitted three-parameter log-logistic curve and what each reported
#' value means. Pure markup, independent of any fit, so it is unit-testable on
#' its own. Rendered as plain HTML (no MathJax) so it works offline.
#' @return A `shiny` tag list.
#' @keywords internal
model_help_single <- function() {
  shiny::tagList(
    shiny::tags$p(
      "The Single Chemical tab fits a three-parameter ",
      shiny::tags$b("log-logistic"), " dose-response curve:"
    ),
    shiny::tags$p(shiny::tags$code(
      "Y = max / (1 + (C / EC50)", shiny::tags$sup("slope"), ")"
    )),
    shiny::tags$p(shiny::tags$small(
      "The response falls from ", shiny::tags$code("max"), " at zero dose."
    )),
    shiny::tags$dl(
      shiny::tags$dt("max"),
      shiny::tags$dd("Control / baseline response at C = 0 (the upper plateau)."),
      shiny::tags$dt("slope"),
      shiny::tags$dd("Steepness of the decline (> 0 means the response decreases with dose)."),
      shiny::tags$dt("EC50"),
      shiny::tags$dd("Concentration that halves the response.")
    ),
    shiny::tags$p(shiny::tags$small(
      shiny::tags$b("SSR"), " = residual sum of squares (goodness of fit); ",
      shiny::tags$b("n"), " = number of data points. ",
      "These describe the fit, not the curve."
    ))
  )
}

#' Single Chemical stage UI
#' @param id Module id.
#' @keywords internal
single_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::uiOutput(ns("errors")),
      shiny::actionButton(ns("fit"), "Fit", class = "btn-primary")
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Dose-response curve"),
                  plotly::plotlyOutput(ns("dr"))),
      bslib::card(bslib::card_header("Observed vs predicted"),
                  plotly::plotlyOutput(ns("op")))
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel("About this model", model_help_single())
    ),
    bslib::card(bslib::card_header("Parameters"), DT::DTOutput(ns("params")))
  )
}

#' Single Chemical stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
single_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() paste0("single_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("single", input$response), file, row.names = FALSE)
    )

    parsed <- shiny::reactive({
      shiny::req(input$file)
      read_upload(input$file$datapath)
    })

    errs <- shiny::reactive({
      shiny::req(input$file)
      validate_upload(parsed(), "single", input$response)
    })

    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    fit_df <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      to_engine_df(parsed(), "single")
    })

    fit_r <- shiny::eventReactive(input$fit, {
      shiny::req(length(errs()) == 0)
      tryCatch(analyse_single(fit_df()), error = function(e) {
        shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
        NULL
      })
    })

    output$dr <- plotly::renderPlotly({
      shiny::req(fit_r())
      p <- plot_dose_response(fit_r(), fit_df())
      plotly::layout(p, xaxis = list(title = axis_label(meta)),
                     yaxis = list(title = if (!is.null(meta$endpoint) && nzchar(meta$endpoint))
                                            meta$endpoint else "Response"))
    })
    output$op <- plotly::renderPlotly({
      shiny::req(fit_r())
      plot_obs_pred(fit_r(), fit_df())
    })
    output$params <- DT::renderDT({
      shiny::req(fit_r())
      p <- fit_r()$par
      DT::datatable(
        data.frame(Parameter = c("max", "slope", "ec50", "SSR", "n"),
                   Value = c(round(unname(p[c("max", "slope", "ec50")]), 4),
                             round(fit_r()$ssr, 2), nrow(fit_df()))),
        rownames = FALSE, options = list(dom = "t"))
    })
  })
}
