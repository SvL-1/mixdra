# Single Stressor stage: upload one stressor's dose-response data, then either
# Autofit a three-parameter log-logistic curve via analyse_single() or Simulate
# the curve for caller-entered parameter values via eval_single(). Shows the
# curve, observed-vs-predicted, an editable parameter grid with bounds, and a
# live SSR readout.

#' Concentration-axis label from the shared meta store
#'
#' With `chem_field` (e.g. "chem1") uses that stressor's name and its matching
#' per-stressor unit (`unit1`/`unit2`/`unit3`). The Single Stressor tab passes
#' `"chem0"`/`"unit0"`, the name and unit entered on that tab. When the name is
#' blank (any field) it falls back to a generic "Concentration"; when only the
#' unit is set the unit is still appended. Used by the single stage and the
#' binary/ternary stages (per stressor).
#' @keywords internal
axis_label <- function(meta, chem_field = NULL) {
  if (is.null(chem_field)) return("Concentration")
  nm   <- meta[[chem_field]]
  unit <- meta[[sub("^chem", "unit", chem_field)]]   # "chem1" -> "unit1"
  base <- if (!is.null(nm) && nzchar(nm)) nm else "Concentration"
  if (!is.null(unit) && nzchar(unit)) paste0(base, " (", unit, ")") else base
}

#' Single Stressor stage UI
#' @param id Module id.
#' @keywords internal
single_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::textInput(ns("name"), "Stressor name"),
      shiny::textInput(ns("unit"), "Stressor unit (e.g. mg/L)"),
      shiny::helpText("Used to label the dose-response axis."),
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

#' Single Stressor stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
single_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    # The single tab's own stressor name + unit feed the dose-response axis via
    # the "chem0"/"unit0" meta fields (axis_label's per-stressor convention).
    shiny::observe({
      meta$chem0 <- input$name
      meta$unit0 <- input$unit
    })

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

    curve_fit_server("curve", fit_df = fit_df, meta = meta, chem_field = "chem0")
  })
}
