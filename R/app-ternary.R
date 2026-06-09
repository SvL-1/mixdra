# Ternary Mixture stage. Reuses the binary tab's UX skeleton: a sidebar, three
# single-curve panels (Stage 1), a per-ratio A4 "hub" table (Stage 2), and an
# inspect stage (Stage 3 -- isoplane + Sigma-TU + the selected ratio's effect).
# The science differs from binary: there is no SA/DR/DL model selection. The
# staged Advanced-S/A fit (analyse_ternary) is the method; its punchline is the
# overall-A4 vs per-ratio-A4 contrast (averaging-out). CA + continuous only.

#' Client-side JS to disable one option of a Shiny radio group
#'
#' Shiny renders radio options as `<input name=... value=...>`; jQuery (bundled
#' with Shiny) flips the unsupported one to disabled on load. Used to show the
#' Quantal / IA options greyed-out without adding a JS dependency.
#' @param input_name The namespaced input id (e.g. `ns("response")`).
#' @param value The option value to disable (e.g. `"quantal"`).
#' @keywords internal
disable_radio_js <- function(input_name, value) {
  shiny::tags$script(shiny::HTML(sprintf(
    "$(function(){$('input[name=\"%s\"][value=\"%s\"]').prop('disabled', true);});",
    input_name, value)))
}

#' Ternary Mixture stage UI
#' @param id Module id.
#' @keywords internal
ternary_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380,
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous",
                            "Quantal (coming soon)" = "quantal"),
                          selected = "continuous"),
      disable_radio_js(ns("response"), "quantal"),
      shiny::radioButtons(ns("reference"), "Reference model",
                          c("Concentration addition (CA)" = "CA",
                            "Independent action (IA) (coming soon)" = "IA"),
                          selected = "CA"),
      disable_radio_js(ns("reference"), "IA"),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example dataset (FBSA + CPF + IMI, continuous) is loaded until you upload your own.")),
      shiny::uiOutput(ns("errors")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced fitting options",
          shiny::helpText("Apply to the staged Advanced-S/A fit."),
          shiny::numericInput(ns("n_starts"), "n_starts", value = 1, min = 1),
          shiny::numericInput(ns("time_limit"), "time_limit (s/stage)", value = 30, min = 1),
          shiny::checkboxInput(ns("thorough"), "Thorough fit (multi-start, slower)", FALSE),
          shiny::uiOutput(ns("thorough_note"))
        )
      )
    ),

    # Stage 1 -- three single-chemical curve panels.
    bslib::card(
      bslib::card_header("Stage 1 · Single curves"),
      shiny::p("Fit each stressor's dose-response curve (Autofit or Simulate). ",
               "These curves are held fixed as the base for the staged ",
               "Advanced-S/A fit below, which appears once all three are fitted."),
      shiny::div(shiny::h5(shiny::textOutput(ns("chem1_title"))),
                 curve_fit_ui(ns("chem1"))),
      shiny::div(class = "mt-4",
                 shiny::h5(shiny::textOutput(ns("chem2_title"))),
                 curve_fit_ui(ns("chem2"))),
      shiny::div(class = "mt-4",
                 shiny::h5(shiny::textOutput(ns("chem3_title"))),
                 curve_fit_ui(ns("chem3"))),
      shiny::uiOutput(ns("reveal_note"))
    ),

    # Stage 2 -- staged Advanced-S/A fit + the per-ratio A4 hub table.
    bslib::card(
      bslib::card_header("Stage 2 · Advanced S/A"),
      shiny::p("Fit all three single curves in Stage 1, then ",
               shiny::tags$b("Fit Advanced S/A"), " to run the staged fit ",
               "(base → pairwise A1/A2/A3 from the binaries → overall A4 ",
               "→ an individual A4 per ternary ratio)."),
      shiny::div(shiny::actionButton(ns("fit_asa"), "Fit Advanced S/A",
                                     class = "btn-primary")),
      shiny::uiOutput(ns("base_pairwise")),
      shiny::helpText(
        "Each row is one ternary ratio; the ", shiny::tags$b("Overall"),
        " row is the single pooled A4. When per-ratio A4 values differ in sign ",
        "but the overall A4 is near zero, real interaction is being averaged ",
        "away — that contrast is the point of the ternary analysis. ",
        "Click a ratio row to inspect it below."),
      DT::DTOutput(ns("hub"))
    ),

    # Stage 3 -- inspect.
    bslib::card(
      bslib::card_header("Stage 3 · Inspect"),
      bslib::layout_columns(
        bslib::card(bslib::card_header("EC50 isoplane (simplex)"),
                    plotly::plotlyOutput(ns("isoplane"), height = "520px")),
        bslib::card(bslib::card_header("ΣTU vs z"),
                    plotly::plotlyOutput(ns("sigma_tu"), height = "520px"))
      ),
      bslib::card(bslib::card_header("Selected ratio — effect size"),
                  shiny::uiOutput(ns("effect")))
    )
  )
}

#' Ternary Mixture stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
ternary_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() "ternary_continuous_template.csv",
      content  = function(file)
        utils::write.csv(template_df("ternary", "continuous"), file, row.names = FALSE)
    )

    output$thorough_note <- shiny::renderUI({
      if (isTRUE(input$thorough))
        shiny::div(class = "text-warning",
                   shiny::tags$small("Multi-start fitting may take several minutes."))
    })

    chem_title <- function(field, n) shiny::renderText({
      nm <- meta[[field]]; if (!is.null(nm) && nzchar(nm)) nm else paste("Stressor", n)
    })
    output$chem1_title <- chem_title("chem1", 1)
    output$chem2_title <- chem_title("chem2", 2)
    output$chem3_title <- chem_title("chem3", 3)

    # Data source: upload, else the bundled FBSA/CPF/IMI example.
    upload_path <- shiny::reactive({
      if (!is.null(input$file)) return(input$file$datapath)
      ex <- system.file("extdata", "ternary_ca_fbsa_cpf_imi_continuous.csv", package = "mixdra")
      if (nzchar(ex)) ex else NULL
    })
    parsed <- shiny::reactive({ shiny::req(upload_path()); read_upload(upload_path()) })
    errs <- shiny::reactive({
      shiny::req(upload_path()); validate_upload(parsed(), "ternary", "continuous")
    })
    output$errors <- shiny::renderUI({
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    engine_df <- shiny::reactive({
      shiny::req(upload_path(), length(errs()) == 0)
      to_engine_df(parsed(), "ternary")
    })
    m1 <- shiny::reactive(marginal_df3(engine_df(), 1))
    m2 <- shiny::reactive(marginal_df3(engine_df(), 2))
    m3 <- shiny::reactive(marginal_df3(engine_df(), 3))

    fit1 <- curve_fit_server("chem1", fit_df = m1, meta = meta, chem_field = "chem1")
    fit2 <- curve_fit_server("chem2", fit_df = m2, meta = meta, chem_field = "chem2")
    fit3 <- curve_fit_server("chem3", fit_df = m3, meta = meta, chem_field = "chem3")

    frozen <- shiny::reactive(!is.null(fit1()) && !is.null(fit2()) && !is.null(fit3()))
    curve_params <- shiny::reactive({
      shiny::req(fit1(), fit2(), fit3())
      assemble_curve_params3(fit1(), fit2(), fit3())
    })

    n_starts_eff <- shiny::reactive(
      if (isTRUE(input$thorough)) max(input$n_starts, 20) else input$n_starts)

    output$reveal_note <- shiny::renderUI({
      if (!isTRUE(frozen()))
        shiny::div(class = "text-muted",
                   shiny::tags$small("Fit all three single curves (Autofit or ",
                                     "Simulate) to enable the Advanced-S/A fit."))
    })

    # The staged result; cleared when the curves / reference / response change.
    asa_res <- shiny::reactiveVal(NULL)
    sel_row <- shiny::reactiveVal(integer(0))
    shiny::observeEvent(list(curve_params(), input$reference, input$response), {
      asa_res(NULL); sel_row(integer(0))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$fit_asa, {
      if (!isTRUE(frozen())) {
        shiny::showNotification(
          "Fit all three single curves first (Autofit / Simulate in Stage 1).",
          type = "warning")
        return()
      }
      res <- shiny::withProgress(message = "Fitting Advanced S/A", value = 0.5, {
        tryCatch(
          analyse_ternary(engine_df(), reference = "CA", response = "continuous",
                          base = curve_params(), n_starts = n_starts_eff(),
                          time_limit = input$time_limit),
          error = function(e) {
            shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
            NULL
          })
      })
      asa_res(res)
      sel_row(1L)                       # default-select the Overall row
    })

    shiny::observeEvent(input$hub_rows_selected, {
      sel_row(input$hub_rows_selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    effect_tbl <- shiny::reactive({
      res <- asa_res(); shiny::req(res)
      ternary_effect_table(res, engine_df())
    })

    # Hub: Overall row + one row per ternary ratio (joined with the effect table).
    hub_df <- shiny::reactive({
      res <- asa_res(); shiny::req(res)
      ind <- res$individual
      eff <- effect_tbl()
      m <- merge(ind, eff[, c("ratio", "pred_CA", "pred_SA", "pred_ASA", "a4_effect")],
                 by = "ratio", all.x = TRUE, sort = FALSE)
      per <- data.frame(
        Ratio    = sprintf("%.2f:%.2f:%.2f", m$C1, m$C2, m$C3),
        n        = m$n,
        A4       = round(m$A4, 4),
        pred_CA  = round(m$pred_CA, 1),
        pred_SA  = round(m$pred_SA, 1),
        pred_ASA = round(m$pred_ASA, 1),
        a4_effect = round(m$a4_effect, 1),
        .key = m$ratio, stringsAsFactors = FALSE, check.names = FALSE)
      overall <- data.frame(
        Ratio = "Overall", n = sum(m$n), A4 = round(res$A4_overall, 4),
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
      res <- asa_res(); shiny::req(res); plot_isoplane(res, engine_df())
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
        sprintf("CA %.1f → CA+S/A %.1f → CA+S/A+S/A %.1f  |  A4 effect %.1f",
                row$pred_CA, row$pred_SA, row$pred_ASA, row$a4_effect))
    })

    list(asa_res = asa_res, hub_df = hub_df, effect_tbl = effect_tbl,
         selected_ratio = selected_ratio, frozen = frozen)  # for testability
  })
}
