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
      shiny::p("Fit each chemical's dose-response curve (Autofit or Simulate). ",
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
