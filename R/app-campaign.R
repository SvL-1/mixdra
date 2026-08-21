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

    singles_server("singles", store)

    # All three slots are registered up front -- a module server cannot be
    # created inside renderUI. Unused slots never receive data. The return
    # value (current_fit/last_compare/... reactives) is kept in `pair_workspaces`
    # purely so tests can reach across the module boundary via testServer's
    # environment access -- production code never reads it back. `on_fit`
    # publishes the pair's S/A value into store$pairs[[k]], which
    # campaign_pairwise() reads to freeze the ternary stage's A1/A2/A3.
    pair_workspaces <- list()
    for (p in list(c(1, 2), c(1, 3), c(2, 3))) local({
      pp <- p
      k  <- pair_key(pp[1], pp[2])
      pair_workspaces[[k]] <<- pair_workspace_server(
        paste0("pair", k),
        fit_df = shiny::reactive({
          shiny::req(store$raw, pair_has_rows(store$raw, pp[1], pp[2]))
          pair_df(store$raw, pp[1], pp[2])
        }),
        base = shiny::reactive({
          b <- campaign_base(store)
          if (is.null(b)) NULL else pair_base(b, pp[1], pp[2])
        }),
        reference    = shiny::reactive(store$reference),
        response     = shiny::reactive(store$response),
        base_version = shiny::reactive(store$base_version),
        on_fit = function(res) store$pairs[[k]] <- res)
    })

    ternary_server("ternary", store)

    invisible(store)
  })
}

#' The pairs of a campaign, in A-B, A-C, B-C order
#' @param chems Integer vector of stressor indices.
#' @return A list of length-2 integer vectors.
#' @keywords internal
campaign_pairs <- function(chems) {
  if (length(chems) < 2) return(list())
  utils::combn(sort(chems), 2, simplify = FALSE)
}

#' Store key for a pair
#' @param i,j Stressor indices.
#' @return A character key, e.g. `"23"`.
#' @keywords internal
pair_key <- function(i, j) paste0(i, j)

#' The campaign's frozen pairwise interaction terms
#'
#' `A1` is the 1-2 pair, `A2` the 1-3 pair, `A3` the 2-3 pair -- the term order
#' `ca_asa_tri()` uses. Each value is that pair's fitted S/A `a`, which is the
#' same quantity (see the design doc, section 3). `NULL` until all three pairs
#' have been fitted, AND `NULL` again if any pair's `a` was fitted against a
#' since-superseded base (its recorded `base_version` no longer matches
#' `store$base_version`) -- otherwise a refit single-stressor curve would
#' silently combine a fresh base with stale pairwise terms, voiding the
#' exactness the campaign's staged reuse depends on (design doc, section 6).
#' @param store The campaign store.
#' @return A named numeric `c(A1, A2, A3)`, or `NULL`.
#' @keywords internal
campaign_pairwise <- function(store) {
  if (length(store$chems) != 3) return(NULL)
  keys <- c("12", "13", "23")
  vals <- lapply(keys, function(k) store$pairs[[k]]$a)
  if (any(vapply(vals, is.null, logical(1)))) return(NULL)
  if (any(vapply(keys, function(k)
        !identical(store$pairs[[k]]$base_version, store$base_version),
        logical(1)))) return(NULL)
  stats::setNames(as.numeric(unlist(vals)), c("A1", "A2", "A3"))
}

#' Are all three pairs fitted, even if stale against the current base?
#'
#' Used only to word the ternary tab's gating message: "fit the pairs" (none
#' fitted yet) vs. "re-fit the pairs" (fitted, but a single-stressor curve
#' changed since). [campaign_pairwise()] alone can't distinguish those --
#' it returns `NULL` for both.
#' @param store The campaign store.
#' @return `TRUE` if all three pairs have a recorded `a`, regardless of version.
#' @keywords internal
campaign_pairs_fitted <- function(store) {
  if (length(store$chems) != 3) return(FALSE)
  all(vapply(c("12", "13", "23"),
             function(k) !is.null(store$pairs[[k]]$a), logical(1)))
}

#' Does this pair have any mixture rows?
#'
#' A campaign may cover a pair's singles without ever dosing them together; that
#' pair's sub-tab is disabled rather than treated as an upload error.
#' @param df Campaign engine frame.
#' @param i,j Stressor indices.
#' @return `TRUE` if at least one row has both concentrations positive.
#' @keywords internal
pair_has_rows <- function(df, i, j) {
  any(df[[paste0("C", i)]] > 0 & df[[paste0("C", j)]] > 0, na.rm = TRUE)
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
  panels <- list(bslib::nav_panel("Singles", singles_ui(ns("singles"))))

  ready <- !is.null(campaign_base(store))
  for (p in campaign_pairs(store$chems)) {
    k     <- pair_key(p[1], p[2])
    title <- paste(axis_label(store, paste0("chem", p[1])), "×",
                   axis_label(store, paste0("chem", p[2])))
    body <- if (!pair_has_rows(store$raw, p[1], p[2])) {
      shiny::div(class = "p-3 text-muted",
                 "This campaign has no mixture rows for this pair.")
    } else if (!ready) {
      shiny::div(class = "p-3 text-muted",
                 sprintf("Fit all %d single-stressor curves first.",
                         store$n_chem))
    } else {
      pair_workspace_ui(ns(paste0("pair", k)))
    }
    panels <- c(panels, list(bslib::nav_panel(title, body)))
  }

  if (length(store$chems) == 3) {
    body <- if (!identical(store$reference, "CA")) {
      shiny::div(class = "p-3 text-muted",
                 "The ternary stage is not yet supported for Independent ",
                 "Action. Switch the campaign reference model to Concentration ",
                 "Addition to use it.")
    } else if (!identical(store$response, "continuous")) {
      shiny::div(class = "p-3 text-muted",
                 "The ternary stage is not yet supported for quantal data.")
    } else if (is.null(campaign_pairwise(store))) {
      shiny::div(class = "p-3 text-muted",
                 if (campaign_pairs_fitted(store))
                   paste("A single-stressor curve changed since these pairs ",
                         "were fitted. Re-run \"Fit interaction models\" on ",
                         "each pair tab to bring them up to date.")
                 else
                   "Fit the interaction models on all three pair tabs first.")
    } else {
      ternary_ui(ns("ternary"))
    }
    panels <- c(panels, list(bslib::nav_panel("Ternary", body)))
  }

  do.call(bslib::navset_card_tab, panels)
}
