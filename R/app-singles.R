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
#' @param bounds Optional list with `lower`/`upper` named numeric vectors over
#'   the base parameters, from [campaign_base_bounds()]. A parameter whose lower
#'   and upper bound are equal is held FIXED at that value. This is how the
#'   Lower/Upper cells and "Fix" boxes on the Singles page reach the campaign
#'   fit; without it the panels would be purely exploratory (issue #12).
#' @return A named numeric vector of base curve parameters.
#' @keywords internal
campaign_fit_base <- function(df, chems, reference, response, bounds = NULL) {
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
  # Route pins (lower == upper) through fit_model()'s `fixed` rather than
  # handing L-BFGS-B a zero-width box, which it rejects.
  sp <- split_fixed_bounds(bounds$lower, bounds$upper, seed)
  fit <- fit_model(singles, reference, "reference", resp, start = sp$start,
                   fixed = sp$fixed, lower = sp$lower, upper = sp$upper,
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

#' Campaign-wide constraint row for the shared upper asymptote
#'
#' `max` is one quantity for the whole campaign (one control group, one
#' asymptote estimated jointly from all single-stressor arms), so it gets ONE
#' constraint row here rather than one per stressor panel: three panels each
#' naming their own `max` bound would be contradictory. The `max` cells inside
#' the panels still constrain those panels' own exploratory Autofit.
#'
#' The input ids match [panel_constraints()]'s `lo_`/`hi_`/`val_`/`pin_`
#' convention so the same reader works on this row and on a panel.
#' @param ns The module's namespace function.
#' @return A `shiny` tag.
#' @keywords internal
campaign_max_row <- function(ns) {
  shiny::div(
    class = "mt-3 mb-2",
    shiny::tags$b("Shared upper asymptote (max)"),
    shiny::helpText(shiny::tags$small(
      "Constrain or fix the experiment-wide max here. The max cells inside ",
      "each ",
      "stressor panel below constrain only that panel's own Autofit.")),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("lo_max"), "Lower constraint",
                                           value = NA)),
      shiny::column(3, shiny::numericInput(ns("hi_max"), "Upper constraint",
                                           value = NA)),
      shiny::column(3, shiny::numericInput(ns("val_max"), "Value", value = NA)),
      shiny::column(3, shiny::checkboxInput(ns("pin_max"), "Fix",
                                            value = FALSE))))
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
             "below, so one stressor has exactly one EC50 for the whole ",
             "experiment."),
    shiny::actionButton(ns("fit_singles"), "Fit single-stressor curves",
                        class = "btn-primary"),
    shiny::helpText(shiny::tags$small(
      "Fits all single-stressor arms together with one shared upper asymptote ",
      "(the experiment has one control group). These values are held fixed by ",
      "every stage below.")),
    campaign_max_row(ns),
    shiny::uiOutput(ns("base_readout")),
    shiny::uiOutput(ns("ec50_note")),
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
    # Each panel's constraint reactive is kept here so the campaign fit below can
    # read what the user typed into that panel -- the module boundary means
    # there is no other way to reach a child module's inputs.
    panel_bounds <- list()
    for (k in 1:3) local({
      kk <- k
      panel <- curve_fit_server(
        paste0("chem", kk),
        fit_df = shiny::reactive({
          shiny::req(store$raw, kk %in% store$chems)
          single_df(store$raw, kk)
        }),
        meta = store, chem_field = paste0("chem", kk))
      fit <- panel$fit
      panel_bounds[[as.character(kk)]] <<- panel$constraints

      # DISPLAY STATE ONLY -- for the fitted VALUES. The panels are exploratory:
      # a scientist may Autofit or Simulate one to interrogate a stressor, and
      # that must never become the campaign's number. The campaign base is the
      # joint fit stored by the "Fit single-stressor curves" button below, so a
      # panel fit deliberately does NOT touch store$base and does NOT bump
      # base_version -- bumping it would invalidate every pair and lock the user
      # out of the ternary for an action the design promises is inert.
      #
      # A panel's CONSTRAINTS are the opposite: a bound or a ticked Fix box is an
      # instruction about how to fit, not a result, so it IS carried into the
      # campaign fit (panel_bounds above; design doc, section 5a as amended).
      shiny::observeEvent(fit(), {
        store$singles[[as.character(kk)]] <- fit()
      }, ignoreNULL = TRUE)
    })

    shiny::observeEvent(input$fit_singles, {
      shiny::req(store$raw, store$chems)
      bounds <- campaign_base_bounds(
        lapply(panel_bounds, function(f) f()),
        panel_constraints(shiny::reactiveValuesToList(input), "max"),
        store$chems)
      b <- tryCatch(
        campaign_fit_base(store$raw, store$chems, store$reference %||% "CA",
                          store$response %||% "continuous", bounds = bounds),
        error = function(e) {
          shiny::showNotification(paste0("Single-stressor fit failed: ",
                                         conditionMessage(e)), type = "error")
          NULL
        })
      if (!is.null(b)) {
        store$base <- b
        # Published alongside the base so each pair's joint "Optimize all"
        # honours the same constraints; without it the joint stage would be the
        # one place that quietly frees a parameter the user pinned.
        store$base_bounds <- bounds
        campaign_bump_base(store)
      }
    })

    output$base_readout <- shiny::renderUI({
      b <- campaign_base(store)
      if (is.null(b)) return(shiny::tags$small("Not fitted yet."))
      pinned <- split_fixed_bounds(store$base_bounds$lower,
                                   store$base_bounds$upper, b)$fixed
      shiny::tagList(
        shiny::tags$small(shiny::HTML(paste0(
          "<b>Fitted base (shared max):</b> ",
          paste(sprintf("%s = %.4g", names(b), b), collapse = " &middot; ")))),
        if (length(pinned))
          shiny::tags$small(shiny::HTML(paste0(
            "<br><b>Held fixed:</b> ", paste(pinned, collapse = ", ")))))
    })

    # An EC50 the joint fit put outside the tested range is flagged HERE, on the
    # page that can repair it: the Fix boxes are a few lines down, and a user who
    # works through the pairs without opening the ternary hub would otherwise
    # never be told that the campaign's toxic units rest on an extrapolation
    # (issue #11).
    output$ec50_note <- shiny::renderUI({
      b <- campaign_base(store)
      shiny::req(b, store$raw, store$chems)
      ec50_range_note(
        store$raw, base_ec50s(b, length(store$chems)),
        vapply(store$chems, function(k) axis_label(store, paste0("chem", k)),
               character(1)),
        chems = store$chems, on_singles = TRUE)
    })
  })
}
