# App assembly + launcher. app_ui()/app_server are internal; run_app() is the
# single exported entry point and guards the Suggested UI dependencies.

#' Assemble the navbar UI from the three stage modules
#' @keywords internal
app_ui <- function() {
  bslib::page_navbar(
    title = "Mixture Toxicity (mixdra)",
    # Tabs use normal document flow and scroll rather than being squeezed into
    # one viewport (the staged Binary tab in particular is tall: two curve
    # panels + two gated stages).
    fillable = FALSE,
    bslib::nav_panel("Introduction", intro_ui("intro")),
    bslib::nav_panel("Single Chemical", single_ui("single")),
    bslib::nav_panel("Binary Mixture", binary_ui("binary")),
    bslib::nav_panel("Ternary Mixture", ternary_ui("ternary"))
  )
}

#' Wire the three stage module servers around a shared meta store
#' @keywords internal
app_server <- function(input, output, session) {
  meta <- shiny::reactiveValues()
  intro_server("intro", meta)
  single_server("single", meta)
  binary_server("binary", meta)
  ternary_server("ternary", meta)
}

#' Launch the mixdra Shiny app
#'
#' Starts the interactive app (Introduction, Single Chemical, Binary Mixture).
#' The UI stack (`shiny`, `bslib`, `plotly`, `DT`) is a set of Suggested
#' dependencies; this function stops with an install hint if any are missing.
#' @param ... Passed to [shiny::runApp()] (e.g. `launch.browser`, `port`).
#' @return Invisibly; called for the side effect of running the app.
#' @export
run_app <- function(...) {
  needed <- c("shiny", "bslib", "plotly", "DT")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("run_app() needs these packages: ", paste(missing, collapse = ", "),
         ". Install with install.packages(c(",
         paste0("\"", missing, "\"", collapse = ", "), ")).", call. = FALSE)
  shiny::runApp(shiny::shinyApp(ui = app_ui(), server = app_server), ...)
}
