# Single Chemical stage: upload one chemical's dose-response data, then either
# Autofit a three-parameter log-logistic curve via analyse_single() or Simulate
# the curve for caller-entered parameter values via eval_single(). Shows the
# curve, observed-vs-predicted, an editable parameter grid with bounds, and a
# live SSR readout.

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
      shiny::helpText(
        "Continuous = a measured amount (column ", shiny::tags$code("Res"),
        "). Quantal = counts (", shiny::tags$code("Affected"), " out of ",
        shiny::tags$code("Exposed"), "); the curve is fit to the proportion."),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::uiOutput(ns("errors"))
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Dose-response curve"),
                  plotly::plotlyOutput(ns("dr"))),
      bslib::card(bslib::card_header("Observed vs predicted"),
                  plotly::plotlyOutput(ns("op")))
    ),
    bslib::card(
      bslib::card_header("Parameters"),
      single_model_equation(),
      shiny::fluidRow(
        shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
        shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Value")))
      ),
      param_row(ns, "max", "max", "Response at C = 0 (control / upper plateau)."),
      param_row(ns, "slope", "slope", "Steepness of the decline (> 0 = decreasing).",
                hi_default = 50),
      param_row(ns, "ec50", "EC50", "Concentration that halves the response."),
      shiny::tags$small(shiny::HTML(
        "SSR = &Sigma; (y &minus; &#375;)&sup2; &nbsp;&nbsp;(&#375; = model prediction)")),
      shiny::uiOutput(ns("diagnostics")),
      shiny::div(
        shiny::actionButton(ns("autofit"), "Autofit parameters", class = "btn-primary"),
        shiny::actionButton(ns("simulate"), "Simulate")
      )
    )
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

    fit_df <- shiny::reactive({
      shiny::req(input$file, length(errs()) == 0)
      to_engine_df(parsed(), "single")
    })

    current_fit <- shiny::reactiveVal(NULL)

    # The Value column read as a named numeric (blank -> NA).
    current_values <- function() {
      raw <- list(max = input$val_max, slope = input$val_slope, ec50 = input$val_ec50)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

    shiny::observeEvent(input$autofit, {
      shiny::req(length(errs()) == 0)
      b <- collect_bounds(shiny::reactiveValuesToList(input), c("max", "slope", "ec50"))
      if (!is.null(b$lower) && !is.null(b$upper)) {
        common <- intersect(names(b$lower), names(b$upper))
        if (length(common) && any(b$lower[common] > b$upper[common])) {
          shiny::showNotification("Lower bound exceeds upper bound.", type = "error")
          return()
        }
      }
      vals <- current_values()
      start <- vals[!is.na(vals)]
      if (!length(start)) start <- NULL
      fit <- tryCatch(
        analyse_single(fit_df(), lower = b$lower, upper = b$upper, start = start),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      shiny::updateNumericInput(session, "val_max",   value = round(fit$par[["max"]], 4))
      shiny::updateNumericInput(session, "val_slope", value = round(fit$par[["slope"]], 4))
      shiny::updateNumericInput(session, "val_ec50",  value = round(fit$par[["ec50"]], 4))
      current_fit(fit)
    })

    shiny::observeEvent(input$simulate, {
      shiny::req(input$file, length(errs()) == 0)
      vals <- current_values()
      if (any(is.na(vals))) {
        shiny::showNotification("Enter max, slope and EC50 to simulate.", type = "warning")
        return()
      }
      resp <- obs_response(fit_df())
      current_fit(eval_single(fit_df()$C1, resp,
                              vals[["max"]], vals[["slope"]], vals[["ec50"]]))
    })

    output$dr <- plotly::renderPlotly({
      shiny::req(current_fit())
      p <- plot_dose_response(current_fit(), fit_df())
      plotly::layout(p, xaxis = list(title = axis_label(meta)),
                     yaxis = list(title = if (!is.null(meta$endpoint) && nzchar(meta$endpoint))
                                            meta$endpoint else "Response"))
    })
    output$op <- plotly::renderPlotly({
      shiny::req(current_fit())
      plot_obs_pred(current_fit(), fit_df())
    })
    output$diagnostics <- shiny::renderUI({
      shiny::req(current_fit())
      shiny::tags$p(
        shiny::tags$b("SSR: "), round(current_fit()$ssr, 2),
        "   |   ", shiny::tags$b("n: "), nrow(fit_df())
      )
    })
  })
}
