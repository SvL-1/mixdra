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
    curve_fit_ui(ns("curve"))
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

    curve_fit_server("curve", fit_df = fit_df, meta = meta)
  })
}
