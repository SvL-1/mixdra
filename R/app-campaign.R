# Campaign stage: ONE upload covering 2-3 stressors, one campaign-wide response
# and reference choice, and the store every downstream stage reads. The stages
# (Singles / Binary pairs / Ternary) are sub-tabs of this one navbar entry --
# Sam's requested layout (issue #10). Single-stressor curves are fitted once on
# the Singles page and held fixed everywhere below it.

#' Default for a NULL value
#'
#' Returns `y` when `x` is `NULL`, otherwise `x`. Used to give campaign store
#' fields a sensible value before the first upload has been read.
#' @param x Value to test.
#' @param y Fallback used when `x` is `NULL`.
#' @return `x`, or `y` when `x` is `NULL`.
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Increment the base-parameter version stamp
#'
#' Downstream stages record the `base_version` they were fitted at; when the
#' stamp no longer matches they show a stale banner instead of silently
#' presenting numbers computed from a superseded curve. Call this whenever a
#' single-stressor fit changes.
#' @param store The campaign store (a [shiny::reactiveValues()]).
#' @return Invisibly, the new version.
#' @keywords internal
campaign_bump_base <- function(store) {
  # isolate() the read: this is called both from inside reactive contexts
  # (campaign_server's observeEvent) and directly in unit tests, and a bare
  # reactiveValues read outside a reactive context always errors in shiny.
  new_version <- (shiny::isolate(store$base_version) %||% 0L) + 1L
  store$base_version <- new_version
  invisible(new_version)
}

#' Campaign stage UI
#' @param id Module id.
#' @keywords internal
campaign_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA)" = "IA")),
      shiny::helpText(shiny::tags$small(
        "One file for the whole campaign: single-stressor rows carry zeros in ",
        "the other concentrations, pair rows carry zero in the third.")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload campaign CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example campaign (FBSA + CPF + IMI, continuous) is loaded until ",
        "you upload your own.")),
      shiny::uiOutput(ns("errors")),
      shiny::uiOutput(ns("summary"))
    ),
    shiny::uiOutput(ns("stages"))
  )
}

#' Campaign stage server
#' @param id Module id.
#' @param meta Shared reactiveValues; the campaign store is built on it so the
#'   flat `chem1..3`/`unit1..3` fields written by [intro_server()] stay
#'   available to [axis_label()] and the curve-fit panels.
#' @return The campaign store, invisibly.
#' @keywords internal
campaign_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {
    store <- meta
    store$base_version <- 0L
    store$singles <- list()
    store$pairs   <- list()

    output$template <- shiny::downloadHandler(
      filename = function() paste0("campaign_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("campaign", input$response), file,
                         row.names = FALSE)
    )

    # The bundled example stands in until the user uploads their own file.
    parsed <- shiny::reactive({
      if (is.null(input$file))
        read_upload(system.file("extdata",
                                "ternary_ca_fbsa_cpf_imi_continuous.csv",
                                package = "mixdra"))
      else read_upload(input$file$datapath)
    })

    errs <- shiny::reactive(
      validate_upload(parsed(), "campaign", input$response %||% "continuous"))

    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    shiny::observe({
      shiny::req(length(errs()) == 0)
      df <- to_engine_df(parsed(), "campaign")
      store$raw       <- df
      store$chems     <- campaign_chems(df)
      store$n_chem    <- campaign_n_chem(df)
      store$response  <- input$response
      store$reference <- input$reference
    })

    # A new file or a changed response/reference invalidates every fit below.
    shiny::observeEvent(list(input$file, input$response, input$reference), {
      store$singles <- list()
      store$pairs   <- list()
      store$ternary <- NULL
      campaign_bump_base(store)
    }, ignoreInit = TRUE)

    output$summary <- shiny::renderUI({
      shiny::req(store$raw)
      cls <- classify_rows(store$raw)
      shiny::tags$small(shiny::HTML(paste0(
        "<b>", store$n_chem, "</b> stressors &middot; ",
        sum(cls == "control"), " control, ", sum(cls == "single"), " single, ",
        sum(cls == "binary"), " pair, ", sum(cls == "ternary"), " ternary rows.")))
    })

    output$stages <- shiny::renderUI(campaign_stage_nav(session$ns, store))

    invisible(store)
  })
}

#' Sub-navigation for the campaign stages
#'
#' Built server-side because which sub-tabs exist depends on the data: the
#' Ternary sub-tab is absent for a two-stressor campaign, and pairs with no
#' mixture rows are disabled.
#' @param ns The module's namespace function.
#' @param store The campaign store.
#' @keywords internal
campaign_stage_nav <- function(ns, store) {
  bslib::navset_card_tab(
    bslib::nav_panel("Singles", singles_ui(ns("singles")))
  )
}
