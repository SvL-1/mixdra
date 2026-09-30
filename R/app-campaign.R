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
#' @noRd
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
        "One file for the whole experiment: single-stressor rows carry zeros in ",
        "the other concentrations, pair rows carry zero in the third.")),
      shiny::downloadButton(ns("template"), "Download template"),
      shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
      shiny::helpText(shiny::tags$small(
        "An example experiment (FBSA + CPF + IMI, continuous) is loaded until ",
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
      filename = function() paste0("multiple_stressors_", input$response,
                                   "_template.csv"),
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

    # A new file or a changed response/reference invalidates EVERY fit below,
    # `store$base` included: the campaign base is derived from the upload, so a
    # base that outlives the file it was fitted on would silently freeze the old
    # file's EC50s into the new file's pair and ternary fits.
    #
    # `priority` puts this ahead of the observer above that writes the new
    # `store$raw`/`chems`/`response`/`reference`, so no downstream stage can
    # observe the new data next to the old base: the whole derived state is
    # gone before a single new value lands.
    shiny::observeEvent(list(input$file, input$response, input$reference), {
      store$base        <- NULL
      store$base_bounds <- NULL
      store$singles <- list()
      store$pairs   <- list()
      store$ternary <- NULL
      campaign_bump_base(store)
    }, ignoreInit = TRUE, priority = 1000)

    output$summary <- shiny::renderUI({
      shiny::req(store$raw)
      cls <- classify_rows(store$raw)
      shiny::tags$small(shiny::HTML(paste0(
        "<b>", store$n_chem, "</b> stressors &middot; ",
        sum(cls == "control"), " control, ", sum(cls == "single"), " single, ",
        sum(cls == "binary"), " pair, ", sum(cls == "ternary"), " ternary rows.")))
    })

    # The sub-navigation's structure depends only on which stressors the file
    # doses, so it is rebuilt only for a new file (see campaign_stage_nav()).
    # Re-read the currently open sub-tab with isolate() and hand it back as
    # `selected` so a rebuild does not snap the user back to Singles: isolate,
    # because a plain read would make every tab click rebuild the whole navset.
    output$stages <- shiny::renderUI(
      campaign_stage_nav(session$ns, store, selected = shiny::isolate(input$stage)))

    # Pair titles follow stressor renames without rebuilding the navset, and
    # each tab's content is re-rendered only when it actually changes
    # (a reactiveVal ignores an identical value): see campaign_stage_nav().
    for (p in list(c(1, 2), c(1, 3), c(2, 3))) local({
      pp <- p
      output[[paste0("title_", pair_key(pp[1], pp[2]))]] <- shiny::renderText(
        paste(axis_label(store, paste0("chem", pp[1])), "×",
              axis_label(store, paste0("chem", pp[2]))))
    })
    for (slot in c("12", "13", "23", "ternary")) local({
      sl   <- slot
      state <- shiny::reactiveVal(NULL)
      shiny::observe(state(campaign_stage_state(store, sl)))
      output[[paste0("body_", sl)]] <- shiny::renderUI({
        st <- state()
        shiny::req(st)
        shiny::isolate(campaign_stage_body(session$ns, store, sl, st))
      })
    })

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
        # The campaign's constraints, renamed onto this pair's binary parameter
        # names, so "Optimize all params (joint)" cannot free a curve parameter
        # the user pinned on the Singles page (issue #12).
        bounds = shiny::reactive(pair_bounds(store$base_bounds, pp[1], pp[2])),
        reference    = shiny::reactive(store$reference),
        # The store holds the USER-facing response key ("continuous"/"quantal")
        # because the gating and the template download read it; the engine's key
        # is "continuous"/"binary". Map here, at the wiring, exactly as
        # campaign_fit_base() does for the singles fit -- fit_model() match.arg()s
        # its `response` and errors on "quantal".
        response     = shiny::reactive(
          if (identical(store$response, "quantal")) "binary" else "continuous"),
        base_version = shiny::reactive(store$base_version),
        on_fit = function(res) store$pairs[[k]] <- res,
        # Label the pair's plots with the stressor names/units and the endpoint
        # from the Introduction tab, as that tab promises.
        labels = shiny::reactive(list(
          x = axis_label(store, paste0("chem", pp[1])),
          y = axis_label(store, paste0("chem", pp[2])),
          z = if (!is.null(store$endpoint) && nzchar(store$endpoint))
                store$endpoint else "response")))
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

#' Does this campaign have any three-stressor mixture rows?
#'
#' The ternary stage's per-ratio `A4` step needs at least one row with all three
#' concentrations positive. Without one no row activates the `A4` term, so
#' `optim` returns its start point and the app would present `A4 = 0` as a
#' fitted result. A campaign covering only the pairwise designs is a legitimate
#' dataset, so this disables the Ternary sub-tab rather than rejecting the
#' upload (the standalone ternary stage's hard error, kept in
#' [validate_upload()], applies to a ternary-only upload).
#' @param df Campaign engine frame.
#' @return `TRUE` if at least one row has `C1`, `C2` and `C3` all positive.
#' @keywords internal
campaign_has_ternary_rows <- function(df) {
  cols <- c("C1", "C2", "C3")
  if (!all(cols %in% names(df))) return(FALSE)
  any(rowSums(as.matrix(df[cols]) > 0) == 3, na.rm = TRUE)
}

#' Sub-navigation for the campaign stages
#'
#' Built server-side because which sub-tabs exist depends on the data: pair
#' tabs for every pair of dosed stressors, and a Ternary tab for a
#' three-stressor campaign.
#'
#' Only the STRUCTURE lives here, and it reads nothing but `store$chems`, so it
#' is rebuilt only when a new file changes the set of stressors. The Singles
#' page is placed statically; the pair titles are `textOutput`s and every other
#' tab's content is a `uiOutput` filled from [campaign_stage_body()]. A
#' rebuild recreates every input inside the navset, so an earlier version that
#' re-rendered on any store write (a fit, a rename) silently cleared the Fix
#' boxes and constraints the user had just set on the Singles page.
#' @param ns The module's namespace function.
#' @param store The campaign store.
#' @param selected Value of the sub-tab to open (`"Singles"`, `"pair12"`, ...,
#'   `"Ternary"`), or `NULL` for the first one. A value that does not exist in
#'   this campaign is ignored.
#' @keywords internal
campaign_stage_nav <- function(ns, store, selected = NULL) {
  panels <- list(bslib::nav_panel("Singles", singles_ui(ns("singles")),
                                  value = "Singles"))
  for (p in campaign_pairs(store$chems)) {
    k <- pair_key(p[1], p[2])
    panels <- c(panels, list(bslib::nav_panel(
      shiny::textOutput(ns(paste0("title_", k)), inline = TRUE),
      shiny::uiOutput(ns(paste0("body_", k))),
      value = paste0("pair", k))))
  }
  if (length(store$chems) == 3)
    panels <- c(panels, list(bslib::nav_panel(
      "Ternary", shiny::uiOutput(ns("body_ternary")), value = "Ternary")))

  values <- c("Singles", paste0("pair", vapply(
    campaign_pairs(store$chems), function(p) pair_key(p[1], p[2]),
    character(1))), if (length(store$chems) == 3) "Ternary")
  sel <- if (!is.null(selected) && selected %in% values) selected else NULL
  do.call(bslib::navset_card_tab,
          c(panels, list(id = ns("stage"), selected = sel)))
}

#' Which state one campaign sub-tab is in
#'
#' A short key for what the tab should show: its workspace (`"ready"`) or one
#' of the reasons it is not available yet. [campaign_server()] keeps this key
#' in a `reactiveVal`, which ignores an identical value, so a tab's content is
#' re-rendered (and its inputs reset) only when the key changes -- not on every
#' store write. The key, not the rendered UI, is compared because the UI
#' carries random tabset ids and is never identical twice.
#' @param store The campaign store.
#' @param slot A pair key (`"12"`, `"13"`, `"23"`) or `"ternary"`.
#' @return A character key; `"absent"` for a slot this campaign does not have.
#' @keywords internal
campaign_stage_state <- function(store, slot) {
  if (identical(slot, "ternary")) {
    if (length(store$chems) != 3) return("absent")
    if (!campaign_has_ternary_rows(store$raw)) return("no_rows")
    if (!identical(store$reference, "CA")) return("not_ca")
    if (!identical(store$response, "continuous")) return("quantal")
    if (is.null(campaign_pairwise(store)))
      return(if (campaign_pairs_fitted(store)) "stale" else "need_pairs")
    return("ready")
  }
  ij <- as.integer(strsplit(slot, "")[[1]])
  if (!all(ij %in% store$chems)) return("absent")
  if (!pair_has_rows(store$raw, ij[1], ij[2])) return("no_rows")
  if (is.null(campaign_base(store))) return(paste0("need_base:", store$n_chem))
  "ready"
}

#' Content of one campaign sub-tab
#'
#' Either the stage's workspace or a note saying why it is not available yet.
#' @param ns The module's namespace function.
#' @param store The campaign store.
#' @param slot A pair key (`"12"`, `"13"`, `"23"`) or `"ternary"`.
#' @param state The tab's [campaign_stage_state()].
#' @return A `shiny` tag, or `NULL` for a slot this campaign does not have.
#' @keywords internal
campaign_stage_body <- function(ns, store, slot,
                                state = campaign_stage_state(store, slot)) {
  note <- function(...) shiny::div(class = "p-3 text-muted", ...)
  if (identical(state, "absent")) return(NULL)
  if (identical(slot, "ternary"))
    return(switch(state,
      no_rows = note("This experiment has no ternary rows (all three stressors ",
                     "dosed together), so the per-ratio A4 step has nothing to fit."),
      not_ca  = note("The ternary stage is not yet supported for Independent ",
                     "Action. Switch the reference model to Concentration ",
                     "Addition to use it."),
      quantal = note("The ternary stage is not yet supported for quantal data."),
      stale   = note("A single-stressor curve changed since these pairs ",
                     "were fitted. Re-run \"Fit interaction models\" on ",
                     "each pair tab to bring them up to date."),
      need_pairs = note("Fit the interaction models on all three pair tabs first."),
      ternary_ui(ns("ternary"))))
  if (identical(state, "no_rows"))
    return(note("This experiment has no mixture rows for this pair."))
  if (startsWith(state, "need_base"))
    return(note(sprintf("Fit all %d single-stressor curves first.",
                        store$n_chem)))
  pair_workspace_ui(ns(paste0("pair", slot)))
}
