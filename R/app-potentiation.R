# (De)potentiation tab: does a second stressor, applied at one or a few fixed
# levels, shift the EC50 of stressor 1? Upload one CSV (Series, C1, Res or
# Exposed/Affected); analyse_potentiation() re-runs whenever the data, the
# reference series or the response type changes (no Fit button).

#' (De)potentiation stage UI
#' @param id Module id.
#' @keywords internal
potentiation_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::textInput(ns("name1"), "Stressor 1 name"),
      shiny::textInput(ns("unit1"), "Stressor 1 unit (e.g. mg/L)"),
      shiny::textInput(ns("name2"), "Stressor 2 name"),
      shiny::helpText("Used to label the axis and the conclusion."),
      shiny::radioButtons(ns("response"), "Response type",
                          c("Continuous" = "continuous", "Quantal" = "quantal")),
      shiny::helpText(
        "Continuous = a measured amount (column ", shiny::tags$code("Res"),
        "). Quantal = counts (", shiny::tags$code("Affected"), " out of ",
        shiny::tags$code("Exposed"), ")."),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::uiOutput(ns("errors")),
      shiny::selectInput(ns("reference"), "Reference series", choices = character(0))
    ),
    shiny::uiOutput(ns("source_note")),
    bslib::card(bslib::card_header("Dose-response curves"),
                plotly::plotlyOutput(ns("plot"))),
    bslib::card(bslib::card_header("Curve parameters per series"),
                DT::DTOutput(ns("curves"))),
    bslib::card(bslib::card_header("Tests against the reference"),
                DT::DTOutput(ns("tests"))),
    bslib::card(bslib::card_header("Controls"), shiny::uiOutput(ns("controls"))),
    bslib::card(bslib::card_header("Conclusion"), shiny::uiOutput(ns("conclusion")))
  )
}

#' Controls note of a (de)potentiation result
#'
#' One sentence per comparison that borrowed the reference's controls; a single
#' "own controls" sentence when none did.
#' @param res Result of [analyse_potentiation()].
#' @return Character vector.
#' @keywords internal
#' @noRd
potentiation_controls_note <- function(res) {
  cmp <- res$comparisons
  b <- which(cmp$controls_borrowed)
  if (!length(b))
    return("Controls: each series uses its own control observations.")
  k <- sum(res$curves[[res$reference]]$data$C1 == 0)
  vapply(b, function(i) paste0(
    "Controls: the ", k, " control observations of ", res$reference,
    " are used for both curves of the comparison with ", cmp$series[[i]],
    ", as in the Excel workbook, so they are counted twice (n = ", cmp$n[[i]],
    "). This is the approach used here; it makes the test slightly more ",
    "liberal than counting them once."), character(1))
}

#' (De)potentiation stage server
#' @param id Module id.
#' @param meta Shared reactiveValues for experiment metadata.
#' @keywords internal
potentiation_server <- function(id, meta) {
  shiny::moduleServer(id, function(input, output, session) {

    output$template <- shiny::downloadHandler(
      filename = function() paste0("potentiation_", input$response, "_template.csv"),
      content  = function(file)
        utils::write.csv(template_df("potentiation", input$response), file,
                         row.names = FALSE)
    )

    use_example <- shiny::reactive(is.null(input$file))

    # The bundled (continuous) example stands in until the user uploads a file.
    parsed <- shiny::reactive({
      if (use_example())
        read_upload(system.file("extdata", "potentiation_example_continuous.csv",
                                package = "mixdra"))
      else read_upload(input$file$datapath)
    })

    # With the example loaded the data are continuous whatever the radio says.
    response <- shiny::reactive(
      if (use_example()) "continuous" else input$response %||% "continuous")

    errs <- shiny::reactive(validate_upload(
      parsed(), "potentiation", input$response %||% "continuous"))

    output$errors <- shiny::renderUI({
      if (use_example()) return(NULL)
      e <- errs()
      if (length(e)) shiny::div(class = "text-danger",
                                lapply(e, function(x) shiny::tags$p(x)))
    })

    output$source_note <- shiny::renderUI({
      if (use_example())
        shiny::helpText("Showing the bundled example (continuous data) until ",
                        "you upload a file.")
    })

    valid <- shiny::reactive(use_example() || length(errs()) == 0)

    series_labels <- shiny::reactive({
      shiny::req(valid())
      unique(trimws(as.character(parsed()$Series)))
    })

    shiny::observe({
      s <- series_labels()
      keep <- shiny::isolate(input$reference)
      shiny::updateSelectInput(session, "reference", choices = s,
        selected = if (!is.null(keep) && keep %in% s) keep else s[[1]])
    })

    name1 <- shiny::reactive(
      if (nzchar(input$name1 %||% "")) input$name1 else "Stressor 1")

    result <- shiny::reactive({
      shiny::req(valid())
      s <- series_labels()
      ref <- input$reference
      if (is.null(ref) || !ref %in% s) ref <- s[[1]]
      analyse_potentiation(parsed(), reference = ref, response = response(),
                           stressor1 = name1())
    })

    controls_note <- shiny::reactive(potentiation_controls_note(result()))

    axis <- shiny::reactive({
      nm <- name1(); u <- input$unit1 %||% ""
      if (nzchar(u)) paste0(nm, " (", u, ")") else nm
    })

    output$plot <- plotly::renderPlotly(
      plot_potentiation(result(), xlab = axis(),
                        ylab = if (response() == "quantal") "proportion affected"
                               else "response"))

    dt <- function(d) DT::datatable(d, rownames = FALSE,
                                    options = list(dom = "t", paging = FALSE))

    output$curves <- DT::renderDT({
      res <- result()
      ref_ec50 <- res$curves[[res$reference]]$par[["ec50"]]
      d <- do.call(rbind, lapply(res$series, function(s) {
        cv <- res$curves[[s]]
        data.frame(Series = s, max = signif(cv$par[["max"]], 4),
                   slope = signif(cv$par[["slope"]], 4),
                   EC50 = signif(cv$par[["ec50"]], 4),
                   `EC50 ratio vs reference` =
                     if (s == res$reference) "reference"
                     else format_num(cv$par[["ec50"]] / ref_ec50),
                   Note = if (isTRUE(cv$ec50_extrapolated)) "extrapolated" else "",
                   check.names = FALSE, stringsAsFactors = FALSE)
      }))
      dt(d)
    })

    output$tests <- DT::renderDT({
      res <- result(); cmp <- res$comparisons
      what <- if (response() == "quantal") "deviance" else "SS"
      d <- data.frame(Series = cmp$series, n = cmp$n,
                      signif(cmp$obj1, 4), signif(cmp$obj2, 4),
                      signif(cmp$chi, 4), p = vapply(cmp$p, format_p, ""),
                      Verdict = cmp$direction,
                      check.names = FALSE, stringsAsFactors = FALSE)
      names(d)[3:5] <- c(paste(what, "run 1"), paste(what, "run 2"), "X²")
      dt(d)
    })

    output$controls <- shiny::renderUI(
      lapply(controls_note(), function(x) shiny::tags$p(x)))

    output$conclusion <- shiny::renderUI(shiny::tags$p(result()$summary))
  })
}
