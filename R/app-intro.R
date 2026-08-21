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
      shiny::helpText("Enter stressor names in the order used throughout the app, ",
                      "with the unit for each."),
      shiny::textInput(ns("chem1"), "Stressor 1 name"),
      shiny::textInput(ns("unit1"), "Stressor 1 unit (e.g. mg/L)"),
      shiny::textInput(ns("chem2"), "Stressor 2 name"),
      shiny::textInput(ns("unit2"), "Stressor 2 unit (e.g. mg/L)"),
      shiny::textInput(ns("chem3"), "Stressor 3 name (ternary)"),
      shiny::textInput(ns("unit3"), "Stressor 3 unit (e.g. mg/L)")
    ),
    bslib::card(
      bslib::card_header("Welcome"),
      shiny::p("Upload concentration-response data on the Single Stressor and ",
               "Campaign tabs. The experiment details entered here are used ",
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
      meta$chem3    <- input$chem3
      meta$unit1    <- input$unit1
      meta$unit2    <- input$unit2
      meta$unit3    <- input$unit3
    })
  })
}
