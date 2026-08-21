# App assembly + launcher. app_ui()/app_server are internal; run_app() is the
# single exported entry point and guards the Suggested UI dependencies.

#' Assemble the navbar UI from the stage modules
#'
#' Introduction, Single Stressor and Campaign are the app's three tabs. The
#' legacy standalone Binary Mixture and Ternary Mixture tabs were retired in
#' favour of the Campaign's per-pair sub-tabs and ternary sub-tab, which reuse
#' the same [pair_workspace_ui()]/[pair_workspace_server()] and
#' [ternary_ui()]/[ternary_server()] modules.
#' @keywords internal
app_ui <- function() {
  bslib::page_navbar(
    title = "Mixture Toxicity (mixdra)",
    # Tabs use normal document flow and scroll rather than being squeezed into
    # one viewport (the Campaign tab in particular is tall: two curve panels +
    # two gated stages per pair).
    fillable = FALSE,
    bslib::nav_panel("Introduction", intro_ui("intro")),
    bslib::nav_panel("Single Stressor", single_ui("single")),
    bslib::nav_panel("Campaign", campaign_ui("campaign"))
  )
}

#' Wire the stage module servers around a shared meta store
#' @keywords internal
app_server <- function(input, output, session) {
  meta <- shiny::reactiveValues()
  intro_server("intro", meta)
  single_server("single", meta)
  campaign_server("campaign", meta)
}

#' Launch the mixdra Shiny app
#'
#' Starts the interactive app: Introduction, Single Stressor and Campaign are
#' the app's three tabs.
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
