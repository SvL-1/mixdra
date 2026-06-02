# Reusable interactive curve-fit panel. Given an injected reactive `fit_df`
# (a single chemical's `C1` + response columns), it fits a 3-parameter
# log-logistic curve via analyse_single() (Autofit) or evaluates caller-entered
# values via eval_single() (Simulate), and shows the curve, observed-vs-predicted,
# an editable parameter grid with bounds, and a live SSR/n readout. It returns a
# reactive of the current fit so a parent (Single or Binary tab) can read it.
# The panel owns no data source and knows nothing about upload validation.

#' Visible model-equation header for the curve-fit panel
#' @return A `shiny` tag.
#' @keywords internal
single_model_equation <- function() {
  shiny::tags$p(
    shiny::tags$b("Model: "),
    shiny::tags$code("Y = max / (1 + (C / EC50)", shiny::tags$sup("slope"), ")")
  )
}

#' One parameter row: label, plain-English meaning, and lower/upper/value inputs
#' @param ns Module namespace function.
#' @param param Parameter key (`max`/`slope`/`ec50`); drives input ids.
#' @param label Display label.
#' @param meaning One-line explanation.
#' @param hi_default Default for the upper-bound input (NA = blank).
#' @keywords internal
param_row <- function(ns, param, label, meaning, hi_default = NA) {
  shiny::fluidRow(
    shiny::column(2, shiny::tags$b(label)),
    shiny::column(4, shiny::tags$small(meaning)),
    shiny::column(2, shiny::numericInput(ns(paste0("lo_", param)), NULL, value = NA)),
    shiny::column(2, shiny::numericInput(ns(paste0("hi_", param)), NULL, value = hi_default)),
    shiny::column(2, shiny::numericInput(ns(paste0("val_", param)), NULL, value = NA))
  )
}

#' One interaction-parameter row: label, meaning, inert bound cells, value input
#'
#' Same five-column layout as [param_row()] for visual consistency with the
#' curve grid, but `a`/`b` are unconstrained by the engine, so the Lower/Upper
#' columns are inert "—" placeholders and only the Value input is editable.
#' @param ns Module namespace function.
#' @param param Parameter key (`a`/`b`); drives the `val_<param>` input id.
#' @param label Display label.
#' @param meaning One-line explanation.
#' @keywords internal
interaction_param_row <- function(ns, param, label, meaning) {
  shiny::fluidRow(
    shiny::column(2, shiny::tags$b(label)),
    shiny::column(4, shiny::tags$small(meaning)),
    shiny::column(2, shiny::tags$small("—")),
    shiny::column(2, shiny::tags$small("—")),
    shiny::column(2, shiny::numericInput(ns(paste0("val_", param)), NULL, value = NA))
  )
}

#' Curve-fit panel UI (plots + parameter grid + Autofit/Simulate)
#' @param id Module id.
#' @keywords internal
curve_fit_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    bslib::layout_columns(
      bslib::card(bslib::card_header("Dose-response curve"),
                  plotly::plotlyOutput(ns("dr"))),
      bslib::card(bslib::card_header("Observed vs predicted"),
                  plotly::plotlyOutput(ns("op")))
    ),
    bslib::card(
      bslib::card_header("Parameters"),
      single_model_equation(),
      shiny::fluidRow(
        shiny::column(2, shiny::tags$small(shiny::tags$b("Parameter"))),
        shiny::column(4, shiny::tags$small(shiny::tags$b("Meaning"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Lower"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Upper"))),
        shiny::column(2, shiny::tags$small(shiny::tags$b("Value")))
      ),
      param_row(ns, "max", "max", "Response at C = 0 (control / upper plateau)."),
      param_row(ns, "slope", "slope", "Steepness of the decline (> 0 = decreasing).",
                hi_default = 50),
      param_row(ns, "ec50", "EC50", "Concentration that halves the response."),
      shiny::tags$small(shiny::HTML(
        "SSR = &Sigma; (y &minus; &#375;)&sup2; &nbsp;&nbsp;(&#375; = model prediction)")),
      shiny::uiOutput(ns("diagnostics")),
      shiny::div(
        shiny::actionButton(ns("autofit"), "Autofit parameters", class = "btn-primary"),
        shiny::actionButton(ns("simulate"), "Simulate")
      )
    )
  )
}

#' Curve-fit panel server
#' @param id Module id.
#' @param fit_df A reactive returning the fit data frame (`C1` + response columns).
#'   The reactive must self-gate (suspend via `req()` / return no rows) when the
#'   data is invalid; the panel fits whatever non-empty frame it is given.
#' @param meta Shared reactiveValues for experiment metadata (axis labels).
#' @param chem_field Optional meta field for the x-axis label (e.g. "chem1"); NULL
#'   uses a generic "Concentration" label.
#' @return A reactive returning the current fit (a [analyse_single()]/[eval_single()]
#'   result), or NULL before any fit.
#' @keywords internal
curve_fit_server <- function(id, fit_df, meta, chem_field = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    current_fit <- shiny::reactiveVal(NULL)

    # The Value column read as a named numeric (blank -> NA).
    current_values <- function() {
      raw <- list(max = input$val_max, slope = input$val_slope, ec50 = input$val_ec50)
      vapply(raw, function(x)
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x), numeric(1))
    }

    shiny::observeEvent(input$autofit, {
      df <- fit_df()
      shiny::req(!is.null(df), nrow(df) > 0)
      b <- collect_bounds(shiny::reactiveValuesToList(input), c("max", "slope", "ec50"))
      if (!is.null(b$lower) && !is.null(b$upper)) {
        common <- intersect(names(b$lower), names(b$upper))
        if (length(common) && any(b$lower[common] > b$upper[common])) {
          shiny::showNotification("Lower bound exceeds upper bound.", type = "error")
          return()
        }
      }
      vals <- current_values()
      start <- vals[!is.na(vals)]
      if (!length(start)) start <- NULL
      fit <- tryCatch(
        analyse_single(df, lower = b$lower, upper = b$upper, start = start),
        error = function(e) {
          shiny::showNotification(paste("Fit failed:", conditionMessage(e)), type = "error")
          NULL
        })
      if (is.null(fit)) return()
      shiny::updateNumericInput(session, "val_max",   value = round(fit$par[["max"]], 4))
      shiny::updateNumericInput(session, "val_slope", value = round(fit$par[["slope"]], 4))
      shiny::updateNumericInput(session, "val_ec50",  value = round(fit$par[["ec50"]], 4))
      current_fit(fit)
    })

    shiny::observeEvent(input$simulate, {
      df <- fit_df()
      shiny::req(!is.null(df), nrow(df) > 0)
      vals <- current_values()
      if (any(is.na(vals))) {
        shiny::showNotification("Enter max, slope and EC50 to simulate.", type = "warning")
        return()
      }
      resp <- obs_response(df)
      current_fit(eval_single(df$C1, resp,
                              vals[["max"]], vals[["slope"]], vals[["ec50"]]))
    })

    output$dr <- plotly::renderPlotly({
      shiny::req(current_fit())
      p <- plot_dose_response(current_fit(), fit_df())
      plotly::layout(p, xaxis = list(title = axis_label(meta, chem_field)),
                     yaxis = list(title = if (!is.null(meta$endpoint) && nzchar(meta$endpoint))
                                            meta$endpoint else "Response"))
    })
    output$op <- plotly::renderPlotly({
      shiny::req(current_fit())
      plot_obs_pred(current_fit(), fit_df())
    })
    output$diagnostics <- shiny::renderUI({
      shiny::req(current_fit())
      shiny::tags$p(
        shiny::tags$b("SSR: "), round(current_fit()$ssr, 2),
        "   |   ", shiny::tags$b("n: "), nrow(fit_df())
      )
    })

    current_fit
  })
}
