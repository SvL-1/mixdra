# ============================================================
# Single Chemical Module
# ============================================================
# This module creates and controls the Single Chemical page.
#
# The page allows the user to:
#   1. Paste concentration-response data for each chemical
#   2. Convert the pasted data into an editable table
#   3. Fit a dose-response model to the data
#   4. View the fitted dose-response curve
#   5. View a summary of the fitted model
#   6. View a plot comparing observed and predicted responses
#
# The module is repeated for three chemicals using a loop.
# Each chemical gets its own input area, plot, summary table, and residual plot.
# ============================================================


# ============================================================
# Paste parser
# ============================================================
# Turns the raw text a user pastes into a clean two-column data frame with
# columns "Conc." and "Res.". It is deliberately forgiving, because users paste
# in several ways:
#   - two columns straight from Excel (tab-separated)
#   - a single comma-separated column copied from a CSV
#   - with or without a header row
#   - with header names that are not exactly "Conc."/"Res."
#
# Strategy: auto-detect the separator from the first line, drop a non-numeric
# first row as a header, then take the first two columns by POSITION and force
# the names to "Conc."/"Res." so the downstream model always finds them.
# On any unrecoverable problem it stop()s with a plain-English message that the
# caller shows to the user.
parse_pasted_dr <- function(text) {
  text <- trimws(text %||% "")
  if (!nzchar(text)) stop("Nothing was pasted. Paste two columns: concentration and response.")

  lines <- strsplit(text, "\r?\n")[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("Nothing was pasted. Paste two columns: concentration and response.")

  # Detect the separator from the first line: tab first, then comma, then any
  # run of whitespace (read.table's sep = "" default).
  first_line <- lines[1]
  sep <- if (grepl("\t", first_line)) "\t" else if (grepl(",", first_line)) "," else ""

  df <- tryCatch(
    read.table(text = paste(lines, collapse = "\n"), header = FALSE,
               sep = sep, stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) stop("Could not read the pasted data. Make sure it is two columns."))

  if (ncol(df) < 2)
    stop("Could not split the data into two columns. Paste a concentration and a response column (tab- or comma-separated).")

  # A first row that is not fully numeric is treated as a header and dropped.
  first_row_numeric <- suppressWarnings(
    all(!is.na(as.numeric(as.character(unlist(df[1, 1:2, drop = TRUE]))))))
  if (!first_row_numeric) df <- df[-1, , drop = FALSE]

  df <- df[, 1:2, drop = FALSE]
  names(df) <- c("Conc.", "Res.")
  df$Conc. <- suppressWarnings(as.numeric(as.character(df$Conc.)))
  df$Res.  <- suppressWarnings(as.numeric(as.character(df$Res.)))
  if (anyNA(df$Conc.) || anyNA(df$Res.))
    stop("The concentration and response columns must contain only numbers (and one optional header row).")

  rownames(df) <- NULL
  df
}

# Local null-coalescing helper (shiny exports `%||%`, but keep the module
# self-contained so the parser works when sourced on its own).
`%||%` <- function(a, b) if (is.null(a)) b else a


# ============================================================
# Single Chemical UI
# ============================================================
# This function defines what the Single Chemical page looks like.
# The same layout is created three times, once for each chemical.

singleChemicalUI <- function(id = "single") {
  
  # Create a namespace for this module.
  # This prevents input and output IDs from conflicting with IDs in other modules.
  ns <- NS(id)
  
  # ---- Single Chemical tab layout ----
  # This tab contains three repeated sections, one for each chemical.
  tabItem(tabName = id,
          
          # Enables shinys functions such as a hide/show inside this module
          useShinyjs(),
          
          # ---- Create one row per chemical ----
          # lapply(1:3, ...) repeats the same UI layout for Chemical 1, 2, and 3.
          # The index idx is used to create unique input and output IDs.
          lapply(1:3, function(idx) {
            fluidRow(
              
              # ---- Left column: data input ----
              # The input box title depends on the chemical name entered in the Introduction tab.
              # Because that name is reactive, this part is rendered with uiOutput().
              column(4,
                     uiOutput(ns(paste0("sgl_input_box_title_", idx)))
              ),
              
              # ---- Middle column: dose-response plot ----
              # This plot shows the observed data points and the fitted model curve.
              column(4,
                     box(title = "Dose–Response Curve", status = "success", width = 12, solidHeader = TRUE, height = "500px",
                         plotOutput(ns(paste0("sgl_plot", idx)))
                     )
              ),
              
              # ---- Right column: model summary and residual plot ----
              # The summary table shows the fitted parameter values.
              # The residual plot compares observed and predicted responses.
              column(4,
                     box(title = "Model Summary", status = "info", width = 12, solidHeader = TRUE, height = "250px",
                         tableOutput(ns(paste0("sgl_summary", idx)))
                     ),
                     box(title = "Residual Plot", status = "warning", width = 12, solidHeader = TRUE, height = "500px",
                         plotOutput(ns(paste0("sgl_resid", idx)))
                     )
              )
            )
          })
  )
}


# ============================================================
# Single Chemical Server
# ============================================================
# This function controls all behaviour on the Single Chemical page.
#
# It handles:
#   1. Dynamic input boxes for each chemical
#   2. Raw data paste-and-convert functionality
#   3. Editable data tables
#   4. Dose-response model fitting
#   5. Plot and summary generation
#   6. Storage of results in sgl_values
#
# The sgl_values object is shared with other modules.
# For example, the binary mixture module can later reuse fitted single-chemical
# parameters as starting values.

singleChemicalServer <- function(id, chem_names, sgl_values) {
  moduleServer(id, function(input, output, session) {
    
    # Use the module namespace for dynamically created IDs.
    ns <- session$ns
    
    
    # ========================================================
    # Dynamic input boxes
    # ========================================================
    # This observe block creates the data input box for each chemical.
    # The box title uses the chemical name entered in the Introduction tab.
    #
    # local({ ... }) is important inside loops.
    # It makes sure each output keeps the correct value of idx.
    # Without local(), all outputs may accidentally use the final loop value.
    observe({
      for (i in 1:3) {
        local({
          idx <- i
          output[[paste0("sgl_input_box_title_", idx)]] <- renderUI({
            
            # Get the chemical name for the current chemical.
            # For example:
            #   idx = 1 uses chem_names$chem1
            #   idx = 2 uses chem_names$chem2
            #   idx = 3 uses chem_names$chem3
            chem_label <- chem_names[[paste0("chem", idx)]]
            
            # ---- Data input box ----
            # This box allows the user to paste raw data, convert it into a table,
            # and edit the table if needed, and fit a dose-response model.
            box(title = paste(chem_label, "Data Input"), status = "primary", width = 12, solidHeader = TRUE,
                div(class = "scroll-box",
                    helpText("Paste data below with two columns: 'Conc.' and 'Res.'"),
                    
                    # Button to show or hide the raw data paste area
                    div(style = "margin-bottom: 10px;",
                        actionButton(ns(paste0("toggle", idx)), "Show/Hide Raw Data")
                    ),
                    
                    # Raw data paste area
                    # The user pastes tab-separated data here.
                    # The expected columns are:
                    #   Conc. = concentration
                    #   Res.  = response
                    div(id = ns(paste0("paste_block", idx)),
                        textAreaInput(ns(paste0("paste", idx)), label = NULL, rows = 60),
                        actionButton(ns(paste0("convert", idx)), "Convert to Table")
                    ),
                    br(), br(),
                    
                    # Editable table created after the pasted data are converted
                    rHandsontableOutput(ns(paste0("sgl_table", idx))),
                    actionButton(ns(paste0("fit", idx)), "Fit Model")
                )
            )
          })
        })
      }
    })
    
    # ========================================================
    # Main server logic for each chemical
    # ========================================================
    # This loop creates the actual actions for each chemical:
    #   - show/hide raw data area
    #   - convert pasted data to a table
    #   - fit the dose-response model
    #   - create plots and summary tables
    #
    # Again, local({ ... }) is used so each chemical keeps its own IDs.
    
    for (i in 1:3) {
      local({
        idx <- i
        
        # ---- Create ID names for the current chemical ----
        # These IDs are used to connect inputs and outputs for each chemical.
        # For example, when idx = 1:
        #   paste_id becomes "paste1"
        #   table becomes "sgl_table1"
        #   plot becomes "sgl_plot1"
        paste_id <- paste0("paste", idx)
        convert <- paste0("convert", idx)
        table <- paste0("sgl_table", idx)
        fit <- paste0("fit", idx)
        plot <- paste0("sgl_plot", idx)
        summary <- paste0("sgl_summary", idx)
        resid <- paste0("sgl_resid", idx)
        toggle <- paste0("toggle", idx)
        paste_block <- paste0("paste_block", idx)
        
        # ====================================================
        # Show or hide raw data paste area
        # ====================================================
        # When the user clicks the toggle button, the area to paste raw data is shown or hidden.
        observeEvent(input[[toggle]], {
          shinyjs::toggle(paste_block)
        })
        
        # ====================================================
        # Convert pasted data into an editable table
        # ====================================================
        # When the user clicks "Convert to Table":
        #   1. The pasted text is read as a data frame
        #   2. The data are stored in sgl_values
        #   3. The raw paste area is hidden
        #   4. The data are displayed as an editable table
        observeEvent(input[[convert]], {

          # Make sure the user pasted something before continuing.
          req(input[[paste_id]])

          # Parse the pasted text into clean Conc./Res. columns. The parser
          # auto-detects tab vs comma and an optional header; on failure it
          # stop()s with a message we show to the user instead of failing
          # silently.
          df <- tryCatch(
            parse_pasted_dr(input[[paste_id]]),
            error = function(e) {
              showNotification(conditionMessage(e), type = "error", duration = 8)
              NULL
            })
          if (is.null(df)) return()

          # Store the data for the current chemical.
          # The index is converted to character because reactiveValues stores named elements.
          sgl_values[[as.character(idx)]]$df <- df

          # Hide the raw paste area after conversion to keep the UI clean.
          shinyjs::hide(paste_block)

          # Show the converted data as an editable table.
          output[[table]] <- renderRHandsontable({
            rhandsontable(df)
          })
        })
        
        # ====================================================
        # Fit dose-response model
        # ====================================================
        # When the user clicks "Fit Model":
        #   1. The current table is read
        #   2. A three-parameter log-logistic model is fitted (by default as we did in the excel model)
        #   3. Predictions are generated
        #   4. Plots and summary tables are created
        #   5. Results are stored in sgl_values
        observeEvent(input[[fit]], {

          # The user must convert the pasted data into a table first.
          if (is.null(input[[table]])) {
            showNotification("Convert the pasted data into a table before fitting.",
                             type = "warning", duration = 8)
            return()
          }

          # Read the current editable table.
          # This means any manual edits made by the user are included.
          df <- tryCatch(hot_to_r(input[[table]]), error = function(e) NULL)

          # Validate the table before fitting so problems are reported, not swallowed.
          if (is.null(df) || !all(c("Conc.", "Res.") %in% names(df))) {
            showNotification("The table needs columns 'Conc.' and 'Res.'. Re-paste and convert the data.",
                             type = "error", duration = 8)
            return()
          }
          df$Conc. <- suppressWarnings(as.numeric(df$Conc.))
          df$Res.  <- suppressWarnings(as.numeric(df$Res.))
          if (anyNA(df$Conc.) || anyNA(df$Res.)) {
            showNotification("Concentration and response values must be numeric.",
                             type = "error", duration = 8)
            return()
          }

          # ---- Fit the dose-response model ----
          # The model uses:
          #   Res.  as the response variable
          #   Conc. as the concentration variable
          #
          # tryCatch() prevents the whole app from crashing if the model fails;
          # we keep the error message so we can show it to the user below.
          fit_err <- NULL
          fit_result <- tryCatch(
            drm(Res. ~ Conc., data = df, fct = LL.3()),
            error = function(e) { fit_err <<- conditionMessage(e); NULL })

          # Extract the model coefficients.
          # If this step fails, coefs becomes NULL.
          coefs <- tryCatch(
            coef(summary(fit_result)),
            error = function(e) { if (is.null(fit_err)) fit_err <<- conditionMessage(e); NULL })

          # Continue only if model coefficients were successfully extracted.
          # The LL.3 model should return three parameters.
          if (!is.null(coefs) && nrow(coefs) >= 3) {
            
            # ---- Get labels for the plots ----
            # Chemical name and unit are taken from the Introduction tab.
            chem_label <- chem_names[[paste0("chem", idx)]]
            unit_label <- chem_names$unit
            
            # Create x-axis label.
            # If both chemical name and unit are available, use both.
            # Otherwise, use a generic label.
            x_axis_label <- if (!is.null(chem_label) && nzchar(chem_label) && !is.null(unit_label) && nzchar(unit_label)) {
              paste0(chem_label, " (", unit_label, ")")
            } else {
              "Concentration"
            }
            
            # ---- Generate predicted dose-response curves ----
            # pred_df contains 100 concentration values between the minimum and maximum tested concentration.
            # These predictions are used to draw a smooth fitted curve.
            pred_df <- expand.grid(Conc. = seq(min(df$Conc.), max(df$Conc.), length.out = 100))
            pred_df$Predicted <- predict(fit_result, newdata = pred_df)
            
            # Add predicted responses for the original observations.
            # These values are used for the observed-versus-predicted plot.
            df$Predicted <- predict(fit_result)
            
            # Store the updated data frame and refresh the editable table.
            sgl_values[[as.character(idx)]]$df <- df
            
            output[[table]] <- renderRHandsontable({
              rhandsontable(df)
            })
            
            
            # ==================================================
            # Dose-response plot
            # ==================================================
            # This plot shows:
            #   - black points: observed responses
            #   - blue dashed line: fitted model prediction
            #
            # The x-axis is shown on a log scale.
            dose_plot <- ggplot(df, aes(x = Conc.)) +
              geom_point(aes(y = Res.), color = "black") +
              geom_line(data = pred_df, aes(x = Conc., y = Predicted), color = "blue", linetype = "dashed") +
              scale_x_log10() +
              scale_y_continuous(limits = c(0, max(df$Res.))) +
              labs(x = x_axis_label, y = if (!is.null(chem_names$endpoint) && nzchar(chem_names$endpoint)) chem_names$endpoint else "Response", title = paste("Dose–Response for", chem_label)) +
              theme_minimal() +
              theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 14)
              )
            
            # ==================================================
            # Observed-versus-predicted plot
            # ==================================================
            # This plot compares the observed response values with the
            # model-predicted response values.
            #
            # The black diagonal line represents 1:1 line
            #
            # Note:
            # Although the box title says "Residual Plot", 
            # this is technically an observed-versus-predicted plot.
            resid_plot <- ggplot(df, aes(x = Res., y = Predicted)) +
              geom_point(color = "darkred") +
              geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
              labs(x = "Observed", y = "Predicted", title = paste("Residual Plot for", chem_label)) +
              xlim(0, 1.2 * max(c(df$Res., df$Predicted), na.rm = TRUE)) +
              ylim(0, 1.2 * max(c(df$Res., df$Predicted), na.rm = TRUE)) +
              theme_minimal() +
              theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 14)
              )
            
            # ==================================================
            # Model summary table
            # ==================================================
            # This table gives the main fitted model information in a
            # user-friendly format.
            #
            # For the LL.3 model used here:
            #   coefs[1, 1] = slope
            #   coefs[2, 1] = upper response limit / maximum response
            #   coefs[3, 1] = EC50
            summary_table <- data.frame(
              Parameter = c("Chemical Name", "Max", "EC50", "Slope", "N (Observations)"),
              Value = c(
                chem_label,
                round(coefs[2, 1], 3),
                round(coefs[3, 1], 3),
                round(coefs[1, 1], 3),
                nrow(df)
              ),
              stringsAsFactors = FALSE
            )
            
            # ==================================================
            # Store all useful results
            # ==================================================
            # These saved values can be displayed in this module and reused globally by other modules.
            sgl_values[[as.character(idx)]]$sgl_dose_plot <- dose_plot
            sgl_values[[as.character(idx)]]$sgl_resid_plot <- resid_plot
            sgl_values[[as.character(idx)]]$sgl_summary_table <- summary_table
            sgl_values[[as.character(idx)]]$df <- df
            
            # Store the fitted parameters in a named vector.
            # These names make it easier to reuse the values later.
            sgl_values[[as.character(idx)]]$fit <- c(Max = coefs[2, 1], EC50 = coefs[3, 1], Slope = coefs[1, 1])
            
            # Store additional model outputs.
            sgl_values[[as.character(idx)]]$n <- nrow(df)
            sgl_values[[as.character(idx)]]$pred_df <- pred_df
            sgl_values[[as.character(idx)]]$residuals <- residuals(fit_result)
          } else {
            
            # ---- If model fitting fails ----
            # Store simple failure values instead of crashing the app.
            # This makes it easier for other modules to check whether a valid model fit is available.
            sgl_values[[as.character(idx)]]$fit <- "FAIL"
            sgl_values[[as.character(idx)]]$n <- NULL
            sgl_values[[as.character(idx)]]$pred_df <- NULL
            sgl_values[[as.character(idx)]]$residuals <- NULL

            # Tell the user why nothing appeared instead of failing silently.
            showNotification(
              paste0("Model fitting failed",
                     if (!is.null(fit_err)) paste0(": ", fit_err) else
                       " (the log-logistic curve could not be fitted to these data)."),
              type = "error", duration = 10)
          }
        })
        
        # ====================================================
        # Render dose-response plot
        # ====================================================
        # This output displays the saved dose-response plot.
        # req() makes sure the plot exists before trying to show it.
        output[[plot]] <- renderPlot({
          req(sgl_values[[as.character(idx)]]$sgl_dose_plot)
          sgl_values[[as.character(idx)]]$sgl_dose_plot
        })
        
        # ====================================================
        # Render model summary table
        # ====================================================
        # This output displays the saved model summary table.
        output[[summary]] <- renderTable({
          req(sgl_values[[as.character(idx)]]$sgl_summary_table)
          sgl_values[[as.character(idx)]]$sgl_summary_table
        })
        
        # ====================================================
        # Render observed-versus-predicted plot
        # ====================================================
        # This output displays the saved model diagnostic plot.
        output[[resid]] <- renderPlot({
          req(sgl_values[[as.character(idx)]]$sgl_resid_plot)
          sgl_values[[as.character(idx)]]$sgl_resid_plot
        })
      })
    }
  })
}
