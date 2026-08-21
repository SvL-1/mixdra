# Singles stage: every stressor's dose-response curve, fitted ONCE, on one page.
# These curves are the campaign's single source of truth -- every pair workspace
# and the ternary stage hold them fixed. Each stressor gets its own curve_fit
# panel, reused verbatim from the standalone Single Stressor tab.

#' The campaign's frozen base-parameter vector
#'
#' `NULL` until every stressor in the campaign has a fitted curve. Returns the
#' binary parameter names (`ec501`/`ec502`) for a two-stressor campaign and the
#' ternary ones (`ec50_1`..`ec50_3`) for a three-stressor campaign, matching the
#' registries the downstream fits use.
#' @param store The campaign store.
#' @return A named numeric vector, or `NULL`.
#' @keywords internal
campaign_base <- function(store) {
  # isolate() the reads: this is called both from reactive contexts and
  # directly in unit tests, and a bare reactiveValues read outside a reactive
  # context always errors in shiny (same rationale as campaign_bump_base).
  chems <- shiny::isolate(store$chems)
  if (is.null(chems)) return(NULL)
  fits <- lapply(as.character(chems), function(k) shiny::isolate(store$singles[[k]]))
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)
  if (length(fits) == 3)
    assemble_curve_params3(fits[[1]], fits[[2]], fits[[3]])
  else
    assemble_curve_params(fits[[1]], fits[[2]])
}

#' Singles stage UI
#' @param id Module id.
#' @keywords internal
singles_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header("Single-stressor curves"),
    shiny::p("Fit each stressor's dose-response curve (Autofit or Simulate). ",
             "These curves are fitted once here and held fixed by every stage ",
             "below, so one stressor has exactly one EC50 for the whole campaign."),
    shiny::uiOutput(ns("panels"))
  )
}

#' Singles stage server
#' @param id Module id.
#' @param store The campaign store.
#' @keywords internal
singles_server <- function(id, store) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$panels <- shiny::renderUI({
      shiny::req(store$chems)
      shiny::tagList(lapply(store$chems, function(k) {
        shiny::div(
          class = if (k == store$chems[1]) "" else "mt-4",
          shiny::h5(axis_label(store, paste0("chem", k))),
          curve_fit_ui(ns(paste0("chem", k))))
      }))
    })

    # One curve_fit_server per stressor. Registered for all three slots up front
    # (a module server cannot be created inside renderUI); slots the campaign
    # does not use simply never receive data.
    for (k in 1:3) local({
      kk <- k
      fit <- curve_fit_server(
        paste0("chem", kk),
        fit_df = shiny::reactive({
          shiny::req(store$raw, kk %in% store$chems)
          single_df(store$raw, kk)
        }),
        meta = store, chem_field = paste0("chem", kk))

      shiny::observeEvent(fit(), {
        store$singles[[as.character(kk)]] <- fit()
        campaign_bump_base(store)
      }, ignoreNULL = TRUE)
    })
  })
}
