# Singles stage: every stressor's dose-response curve, fitted ONCE, on one page.
# These curves are the campaign's single source of truth -- every pair workspace
# and the ternary stage hold them fixed. Each stressor gets its own curve_fit
# panel, reused verbatim from the standalone Single Stressor tab.

#' Joint single-stressor fit for a campaign
#'
#' Fits every single-stressor arm together with ONE shared `max`, exactly as the
#' engine's own Stage 1 does (`fit_ternary_asa()`), rather than fitting each
#' stressor separately and averaging their asymptotes. A campaign has one control
#' group, so the upper asymptote is a single quantity estimated once from all the
#' single-stressor arms — and this is what the reference workbook does, which
#' keeps campaign numbers comparable with Excel.
#' @param df Campaign engine frame (`C1`, `C2`, optionally `C3`, response cols).
#' @param chems Integer stressor indices, from [campaign_chems()].
#' @param reference "CA" or "IA".
#' @param response "continuous" or "quantal" (mapped to the engine's "binary").
#' @return A named numeric vector of base curve parameters.
#' @keywords internal
campaign_fit_base <- function(df, chems, reference, response) {
  # The engine infers the stressor count from the concentration columns PRESENT,
  # so drop any column this campaign never doses -- otherwise a leftover all-zero
  # C3 makes a two-stressor campaign fit as a three-stressor one.
  keep <- c(paste0("C", chems),
            intersect(c("Res", "Affected", "Exposed"), names(df)))
  d <- df[keep]
  names(d)[seq_along(chems)] <- paste0("C", seq_along(chems))

  resp <- if (identical(response, "quantal")) "binary" else "continuous"
  cls  <- classify_rows(d)
  singles <- d[cls %in% c("control", "single"), , drop = FALSE]

  params <- model_spec(reference, "reference", length(chems))$params
  seed   <- seed_from_singles(d, resp)
  fit <- fit_model(singles, reference, "reference", resp, start = seed,
                   n_starts = 1, time_limit = 30)
  fit$par[params]
}

#' The campaign's frozen base-parameter vector
#'
#' The result of the joint single-stressor fit, stored by [singles_server()] when
#' the user runs it. `NULL` until then. Downstream stages read this and never
#' re-derive it.
#' @param store The campaign store.
#' @return A named numeric vector, or `NULL`.
#' @keywords internal
campaign_base <- function(store) store$base

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
    shiny::actionButton(ns("fit_singles"), "Fit single-stressor curves",
                        class = "btn-primary"),
    shiny::helpText(shiny::tags$small(
      "Fits all single-stressor arms together with one shared upper asymptote ",
      "(the campaign has one control group). These values are held fixed by ",
      "every stage below.")),
    shiny::uiOutput(ns("base_readout")),
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

      # DISPLAY STATE ONLY. The panels are exploratory: a scientist may Autofit
      # or Simulate one to interrogate a stressor, and that must never become
      # the campaign's number. The campaign base is the joint fit stored by the
      # "Fit single-stressor curves" button below, so a panel fit deliberately
      # does NOT touch store$base and does NOT bump base_version -- bumping it
      # would invalidate every pair and lock the user out of the ternary for an
      # action the design promises is inert (design doc, section 5a).
      shiny::observeEvent(fit(), {
        store$singles[[as.character(kk)]] <- fit()
      }, ignoreNULL = TRUE)
    })

    shiny::observeEvent(input$fit_singles, {
      shiny::req(store$raw, store$chems)
      b <- tryCatch(
        campaign_fit_base(store$raw, store$chems, store$reference %||% "CA",
                          store$response %||% "continuous"),
        error = function(e) {
          shiny::showNotification(paste0("Single-stressor fit failed: ",
                                         conditionMessage(e)), type = "error")
          NULL
        })
      if (!is.null(b)) {
        store$base <- b
        campaign_bump_base(store)
      }
    })

    output$base_readout <- shiny::renderUI({
      b <- campaign_base(store)
      if (is.null(b)) return(shiny::tags$small("Not fitted yet."))
      shiny::tags$small(shiny::HTML(paste0(
        "<b>Campaign base (shared max):</b> ",
        paste(sprintf("%s = %.4g", names(b), b), collapse = " &middot; "))))
    })
  })
}
