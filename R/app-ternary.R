# Ternary Mixture stage: a sub-tab of the Campaign, appearing once all three
# pair tabs have been fitted. It owns NO upload, curve fitting, or pairwise
# fitting of its own -- those are the Singles page and the three pair tabs.
# Instead it consumes the campaign's frozen single-stressor curves
# (campaign_base()) and each pair's frozen S/A interaction term
# (campaign_pairwise()), and runs the two remaining stages of the staged
# Advanced-S/A fit: an overall three-way A4 from all the data, and an
# individual A4 per ternary mixture ratio. This reuse was requested by the
# project's scientist and proven exact: each pair's fitted `a` equals the
# ternary model's corresponding A-term to 5-6 significant figures on real
# data (see the design doc, section 3). CA + continuous only --
# campaign_stage_nav() shows an explanatory panel instead of this UI for
# IA / quantal / a not-yet-fully-fitted campaign.

#' Ternary Mixture stage UI
#' @param id Module id.
#' @keywords internal
ternary_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("sa_note")),

    bslib::card(
      bslib::card_header("Advanced S/A (frozen base + pairwise)"),
      shiny::p("Uses the experiment's frozen single-stressor curves and each ",
               "pair's frozen S/A interaction term (A1 = pair 1-2, A2 = pair ",
               "1-3, A3 = pair 2-3) to fit the overall three-way A4 and an ",
               "individual A4 per ternary mixture ratio."),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced fitting options",
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1)
        )
      ),
      shiny::uiOutput(ns("base_pairwise")),
      shiny::helpText(
        "Each row is one ternary ratio; the ", shiny::tags$b("Overall"),
        " row is the single pooled A4. When per-ratio A4 values differ in sign ",
        "but the overall A4 is near zero, real interaction is being averaged ",
        "away \u2014 that contrast is the point of the ternary analysis. ",
        "Click a ratio row to inspect it below."),
      shiny::helpText(shiny::tags$small(
        shiny::tags$b("Ratio"), " is the concentration ratio as dosed; ",
        shiny::tags$b("TU ratio"), " is the same mixture in toxic units ",
        "(TU = C / EC50 from the singles), i.e. each stressor's share of the ",
        "mixture's toxicity rather than of its mass. The two differ whenever ",
        "the stressors differ in potency; the models are fitted on the TU ",
        "shares.")),
      DT::DTOutput(ns("hub")),
      shiny::uiOutput(ns("ec50_note"))
    ),

    bslib::card(
      bslib::card_header("Inspect"),
      bslib::layout_columns(
        bslib::card(bslib::card_header("EC50 isoplane (simplex)"),
                    plotly::plotlyOutput(ns("isoplane"), height = "520px")),
        bslib::card(bslib::card_header("\u03a3TU vs z"),
                    plotly::plotlyOutput(ns("sigma_tu"), height = "520px"))
      ),
      bslib::card(bslib::card_header("Selected ratio \u2014 effect size"),
                  shiny::uiOutput(ns("effect")))
    )
  )
}

#' Ternary Mixture stage server
#'
#' Consumes the campaign store: [campaign_base()] for the frozen
#' single-stressor curves and [campaign_pairwise()] for the frozen per-pair
#' S/A terms. Fits automatically once both are available -- there is no
#' local upload or curve/pairwise fitting to trigger first; Singles and the
#' pair tabs already did that. `store$raw`/`store$reference`/`store$response`
#' re-gate the fit to CA + continuous even though `campaign_stage_nav()`
#' already shows an explanatory panel instead of this UI otherwise, so a
#' stale mounted panel never silently fits.
#' @param id Module id.
#' @param store The campaign store.
#' @keywords internal
ternary_server <- function(id, store) {
  shiny::moduleServer(id, function(input, output, session) {

    # Notes which pair(s) did NOT select S/A as their winning model -- that
    # can be DR, DL, or "reference" (no interaction at all). The three-way
    # model has no DR/DL form and no no-interaction shortcut, so it always
    # uses that pair's S/A term regardless (agreed with the scientist). Lives
    # here, not on the pair tab, because it describes the ternary fit.
    output$sa_note <- shiny::renderUI({
      odd <- names(Filter(function(p) !identical(p$chosen, "SA"), store$pairs))
      if (length(odd))
        shiny::div(class = "alert alert-info",
                   "Pair(s) ", paste(odd, collapse = ", "),
                   " were not best described by the Synergism/Antagonism ",
                   "(S/A) model (a dose-ratio or dose-level model fit better, ",
                   "or no interaction at all). The ternary fit uses their S/A ",
                   "term anyway, which is the only interaction form the ",
                   "three-way model has.")
    })

    asa_res <- shiny::reactive({
      shiny::req(store$raw, store$reference == "CA",
                 store$response == "continuous")
      base <- campaign_base(store); shiny::req(base)
      pw   <- campaign_pairwise(store); shiny::req(pw)
      analyse_ternary(store$raw, reference = "CA", response = "continuous",
                      base = base, pairwise = pw,
                      n_starts = input$n_starts %||% 1)
    })

    sel_row <- shiny::reactiveVal(integer(0))
    shiny::observeEvent(asa_res(), sel_row(1L), ignoreNULL = TRUE)  # default: Overall
    shiny::observeEvent(input$hub_rows_selected, {
      sel_row(input$hub_rows_selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    effect_tbl <- shiny::reactive({
      res <- asa_res(); shiny::req(res)
      ternary_effect_table(res, store$raw)
    })

    # Hub: Overall row + one row per ternary ratio (joined with the effect table).
    hub_df <- shiny::reactive({
      res <- asa_res(); shiny::req(res)
      ind <- res$individual
      eff <- effect_tbl()
      m <- merge(ind, eff[, c("ratio", "pred_CA", "pred_SA", "pred_ASA", "a4_effect")],
                 by = "ratio", all.x = TRUE, sort = FALSE)
      ec <- c(res$base[["ec50_1"]], res$base[["ec50_2"]], res$base[["ec50_3"]])
      tu <- vapply(seq_len(nrow(m)),
                   function(i) ratio_label(
                     tu_shares(c(m$C1[i], m$C2[i], m$C3[i]), ec)),
                   character(1))
      per <- data.frame(
        Ratio    = sprintf("%.2f:%.2f:%.2f", m$C1, m$C2, m$C3),
        `TU ratio` = tu,
        n        = m$n,
        A4       = round(m$A4, 4),
        pred_CA  = round(m$pred_CA, 1),
        pred_SA  = round(m$pred_SA, 1),
        pred_ASA = round(m$pred_ASA, 1),
        a4_effect = round(m$a4_effect, 1),
        .key = m$ratio, stringsAsFactors = FALSE, check.names = FALSE)
      overall <- data.frame(
        Ratio = "Overall", `TU ratio` = "\u2014",
        n = sum(m$n), A4 = round(res$A4_overall, 4),
        pred_CA = NA_real_, pred_SA = NA_real_, pred_ASA = NA_real_,
        a4_effect = NA_real_, .key = NA_character_,
        stringsAsFactors = FALSE, check.names = FALSE)
      rbind(overall, per)
    })

    selected_ratio <- shiny::reactive({
      h <- hub_df(); r <- sel_row()
      if (!length(r) || r[[1]] < 1 || r[[1]] > nrow(h)) return(NULL)
      k <- h$.key[r[[1]]]
      if (is.na(k)) NULL else k
    })

    output$base_pairwise <- shiny::renderUI({
      res <- asa_res(); shiny::req(res)
      b <- res$base; p <- res$pairwise
      shiny::tags$p(
        shiny::tags$b("Base: "),
        sprintf("max %.1f | slope %.2f/%.2f/%.2f | EC50 %.3f/%.3f/%.3f",
                b[["max"]], b[["slope1"]], b[["slope2"]], b[["slope3"]],
                b[["ec50_1"]], b[["ec50_2"]], b[["ec50_3"]]),
        shiny::tags$br(),
        shiny::tags$b("Pairwise: "),
        sprintf("A1 %.3f | A2 %.3f | A3 %.3f", p[["A1"]], p[["A2"]], p[["A3"]]))
    })

    # A TU share is only as good as the EC50 in its denominator: when a fitted
    # EC50 sits above every dose actually tested for that stressor it is an
    # extrapolation, and so is its contribution to the TU ratio.
    output$ec50_note <- shiny::renderUI({
      res <- asa_res(); shiny::req(res, store$raw)
      ec  <- c(res$base[["ec50_1"]], res$base[["ec50_2"]], res$base[["ec50_3"]])
      bad <- which(ec50_out_of_range(store$raw, ec, 1:3))
      if (!length(bad)) return(NULL)
      shiny::div(
        class = "alert alert-warning mt-2",
        shiny::tags$small(
          "Extrapolated EC50 for ",
          paste(vapply(bad, function(k) axis_label(store, paste0("chem", k)),
                       character(1)), collapse = ", "),
          ": the fitted EC50 lies above every dose tested for that stressor, ",
          "so its share of the TU ratio is a model extrapolation rather than a ",
          "measured quantity."))
    })

    output$hub <- DT::renderDT({
      h <- hub_df()
      hidekey <- which(names(h) == ".key") - 1L     # 0-based
      dt <- DT::datatable(
        h, rownames = FALSE,
        selection = list(mode = "single", target = "row",
                         selected = shiny::isolate(sel_row())),
        options = list(dom = "t", ordering = FALSE, scrollX = TRUE,
                       columnDefs = list(list(visible = FALSE, targets = hidekey))))
      DT::formatStyle(dt, "Ratio", valueColumns = "Ratio", target = "row",
                      fontWeight = DT::styleEqual("Overall", "bold"),
                      backgroundColor = DT::styleEqual("Overall", "#d8f0d8"))
    })

    output$isoplane <- plotly::renderPlotly({
      res <- asa_res(); shiny::req(res)
      plot_isoplane(res, store$raw, labels = list(
        x = axis_label(store, "chem1"),
        y = axis_label(store, "chem2"),
        z = axis_label(store, "chem3")))
    })
    output$sigma_tu <- plotly::renderPlotly({
      res <- asa_res(); shiny::req(res); plot_sigma_tu(res)
    })

    output$effect <- shiny::renderUI({
      res <- asa_res(); shiny::req(res)
      k <- selected_ratio()
      if (is.null(k))
        return(shiny::tags$p(shiny::tags$small(
          "Select a ratio row above to see its near-EC50 effect size. The ",
          "Overall row pools every ternary mixture.")))
      eff <- effect_tbl(); row <- eff[eff$ratio == k, , drop = FALSE]
      shiny::req(nrow(row) == 1)
      shiny::tags$p(
        shiny::tags$b("At the near-EC50 point "),
        sprintf("(C1=%.3g, C2=%.3g, C3=%.3g): ", row$C1, row$C2, row$C3),
        sprintf("CA %.1f \u2192 CA+S/A %.1f \u2192 CA+S/A+S/A %.1f  |  A4 effect %.1f",
                row$pred_CA, row$pred_SA, row$pred_ASA, row$a4_effect))
    })

    list(asa_res = asa_res, hub_df = hub_df, effect_tbl = effect_tbl,
         selected_ratio = selected_ratio)  # for testability
  })
}
