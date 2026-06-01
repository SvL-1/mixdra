# Introduction stage: collects experiment metadata into the shared `meta` store.
# `meta` is a reactiveValues created by app_server() and shared with every stage,
# which read it only for plot titles and axis labels.

#' Introduction stage UI
#' @param id Module id.
#' @keywords internal
intro_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::textInput(ns("expID"), "Experiment ID"),
      shiny::textInput(ns("species"), "Species"),
      shiny::textInput(ns("endpoint"), "Endpoint"),
      shiny::helpText("Enter chemical names in the order used throughout the app."),
      shiny::textInput(ns("chem1"), "Chemical 1 name"),
      shiny::textInput(ns("chem2"), "Chemical 2 name"),
      shiny::textInput(ns("unit"), "Concentration unit (e.g. mg/L)")
    ),
    bslib::card(
      bslib::card_header("Welcome"),
      shiny::p("Upload concentration-response data on the Single Chemical and ",
               "Binary Mixture tabs. The experiment details entered here are used ",
               "to label plots and axes.")
    )
  )
}

#' Introduction stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
intro_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      meta$expID    <- input$expID
      meta$species  <- input$species
      meta$endpoint <- input$endpoint
      meta$chem1    <- input$chem1
      meta$chem2    <- input$chem2
      meta$unit     <- input$unit
    })
  })
}
