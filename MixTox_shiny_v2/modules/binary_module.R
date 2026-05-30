# ============================================================
# Binary Mixture Module
# ============================================================
# This file contains the UI and server logic for binary mixture analysis.
#
# The binary mixture module handles three pairwise chemical combinations:
#   1. Chemical 1 + Chemical 2  -> "12"
#   2. Chemical 1 + Chemical 3  -> "13"
#   3. Chemical 2 + Chemical 3  -> "23"
#
# For each binary mixture, the user can:
#   - Paste binary mixture data
#   - Convert the pasted data into an editable table
#   - Enter or import model parameter values
#   - Fit mixture models
#
# NOTE:
# This module is still under development.
# Some parts are already functional, while others are placeholders or are not fully connected yet.
# ============================================================


# ============================================================
# Get starting values from fitted single-chemical models
# ============================================================
# This helper function retrieves the fitted parameters from the
# Single Chemical module and uses them as starting values for binary mixture model fitting, as we did in the excel model
#
# The binary mixture model needs starting values for:
#   - max    = maximum response
#   - slope1 = slope for the first chemical in the mixture
#   - ec501  = EC50 for the first chemical in the mixture
#   - slope2 = slope for the second chemical in the mixture
#   - ec502  = EC50 for the second chemical in the mixture
#
# The function uses bin_idx to determine which two single-chemical models should be used.
#
# For example:
#   bin_idx = "12" uses fitted values from Chemical 1 and Chemical 2
#   bin_idx = "13" uses fitted values from Chemical 1 and Chemical 3
#   bin_idx = "23" uses fitted values from Chemical 2 and Chemical 3

get_initials_from_singles <- function(bin_idx, sgl_values) {
  
  # Map each binary mixture index to the corresponding two single chemicals.
  idx_map <- list("12" = c(1, 2), "13" = c(1, 3), "23" = c(2, 3))
  
  # Extract the two chemical indices for the current binary mixture.
  s1 <- idx_map[[bin_idx]][1]
  s2 <- idx_map[[bin_idx]][2]
  
  # ---- Check whether fitted single-chemical values exist ----
  # If the single-chemical model has not been fitted yet,
  # there are no starting values available.
  if (is.null(sgl_values[[as.character(s1)]]$fit)) {
    cat("❌ sgl_values[[", s1, "]]$fit is NULL\n")
    return(NULL)
  }
  if (is.null(sgl_values[[as.character(s2)]]$fit)) {
    cat("❌ sgl_values[[", s2, "]]$fit is NULL\n")
    return(NULL)
  }
  
  # Retrieve fitted values for the two single chemicals.
  fit1 <- sgl_values[[as.character(s1)]]$fit
  fit2 <- sgl_values[[as.character(s2)]]$fit
  
  # Print diagnostic information to the R console.
  # This is useful while debugging, because it confirms which values
  # are being pulled from the single-chemical module.
  cat(paste0("[get_initials_from_singles] bin_idx: ", bin_idx, "\n"))
  cat("  fit1:", toString(fit1), "\n")
  cat("  fit2:", toString(fit2), "\n")
  
  # ---- Build initial values for the binary model ----
  # The max value is taken as the average of the two single-chemical max values.
  # Slopes and EC50 values are taken directly from the corresponding single-chemical model fits.
  initials <- list(
    max    = mean(c(as.numeric(fit1["Max"]), as.numeric(fit2["Max"])), na.rm = TRUE),
    slope1 = as.numeric(fit1["Slope"]),
    ec501  = as.numeric(fit1["EC50"]),
    slope2 = as.numeric(fit2["Slope"]),
    ec502  = as.numeric(fit2["EC50"])
  )
  
  # Print the final initial values to the console.
  cat("  initials (named):\n")
  for (param in names(initials)) {
    cat(sprintf("    %s: %g\n", param, initials[[param]]))
  }
  
  return(initials)
}

# ============================================================
# Vectorized binary CA model
# ============================================================
# CA_bi is the binary concentration addition model defined in model_functions.R.
#
# Vectorize() allows CA_bi to work on vectors of C1 and C2 values,
# rather than only one pair of concentrations at a time.
#
# This is needed because the model has to predict responses for all rows in the dataset at once.
CA_bi_vec <- Vectorize(CA_bi, vectorize.args = c("C1", "C2"))


# ============================================================
# Print SSR comparison
# ============================================================
# SSR means sum of squared residuals.
#
# In simple terms:
#   SSR measures how far the model predictions are from the observed data.
#   Smaller SSR means the model fits the data better.
#
# This helper function compares:
#   - SSR using the starting parameter values
#   - SSR using the optimized parameter values
#
# It is mainly a diagnostic function for checking whether model fitting actually improved the fit.
print_ssr_comparison <- function(start_vals, optimized_vals, C1, C2, Res, model_fun) {
  
  # Make sure the parameter vectors have consistent names.
  names(start_vals)     <- names(optimized_vals) <- c("max","slope1","slope2","ec501","ec502")
  
  # Generate predictions using the starting values.
  pred_start <- do.call(model_fun, c(list(C1=C1, C2=C2), as.list(start_vals)))
  
  # Generate predictions using the optimized values.
  pred_optim <- do.call(model_fun, c(list(C1=C1, C2=C2), as.list(optimized_vals)))
  
  # Print both SSR values to the console.
  cat(sprintf("Initial SSR:   %.4f\n", sum((Res - pred_start)^2)))
  cat(sprintf("Optimized SSR: %.4f\n", sum((Res - pred_optim)^2)))
}

# ============================================================
# Fit the Single Chemical model inside the Binary Mixture module
# ============================================================
# Despite the function name, this does not refit the original single-chemical dose-response curves.
#
# Instead, it fits a binary mixture model using the single-chemical parameters as a starting point.
#
# The fitting is currently done in two steps:
#
# Step 1:
#   Fit slope1, slope2, ec501, and ec502 using only treatment data.
#   The max parameter is kept fixed.
#
# Step 2:
#   Fit max using all data.
#   The slope and EC50 parameters from Step 1 are kept fixed.
#
# A toggle switch is define for the user to choose if parameters need to be fixed
#
# This part contains a lot of diagnostic steps for coding

fit_single_model <- function(df, init_vals, fix_flags, model_fun) {
  
  # ---- Diagnostic printout ----
  # Print the starting values received by the function.
  # This is useful for checking whether values were correctly imported
  # from the single-chemical module or manually overridden by the user.
  cat(">>> fit_single_model() received init_vals:\n")
  print(init_vals)
  print(str(init_vals))
  
  # Convert the initial values from a list into a named numeric vector.
  start_vals <- unlist(init_vals)
  
  # Identify which parameters are free to optimize.
  # NOTE:
  # This is currently calculated, but not actually used later.
  # So right now, the fix_flags are not fully respected.
  free_idx   <- which(!unlist(fix_flags))
  par_names  <- names(start_vals)
  
  # ---- Find the position of each model parameter ----
  # These indices make it easier to update specific parameters
  # during optimization.
  idx_max    <- which(par_names == "max")
  idx_slope1 <- which(par_names == "slope1")
  idx_slope2 <- which(par_names == "slope2")
  idx_ec501  <- which(par_names == "ec501")
  idx_ec502  <- which(par_names == "ec502")
  
  # ---- Split the dataset into controls and treatments ----
  # Controls are rows where both chemical concentrations are zero, or non-zero if user uses measured concentration.
  # Should be updated to accommdate this possibility.
  # Treatments are all other rows.
  is_control    <- (df$C1 == 0 & df$C2 == 0)
  is_treatment  <- !is_control
  
  # ========================================================
  # Step 1: Fit slopes and EC50 values using treatment data
  # ========================================================
  # In this step:
  #   - max is fixed
  #   - slope1, slope2, ec501, and ec502 are optimized
  #
  # Only treatment rows are used, meaning control rows are excluded.
  par1         <- start_vals
  
  # Parameters to optimize in Step 1.
  opt_idx1     <- c(idx_slope1, idx_slope2, idx_ec501, idx_ec502)
  
  # Starting values for the Step 1 optimizer.
  theta0_1     <- start_vals[opt_idx1]
  
  # Set lower and upper bounds for optimization.
  # The current bounds allow each parameter to vary from 10% to 1000% of its starting value.
  lb1          <- pmax(1e-8, theta0_1 * 0.1)
  ub1          <- pmax(lb1 * 1.01, theta0_1 * 10)
  
  # ---- Print Step 1 starting values ----
  cat("\n[Step 1: Slopes/EC50s Only] Initial values:\n")
  cat(sprintf("  max    = %.4f\n", par1["max"]))
  cat(sprintf("  slope1 = %.4f\n  slope2 = %.4f\n  ec501  = %.4f\n  ec502  = %.4f\n", 
              par1["slope1"], par1["slope2"], par1["ec501"], par1["ec502"]))
  
  # Calculate SSR before optimization for Step 1.
  preds_init1 <- model_fun(
    C1 = df$C1[is_treatment], C2 = df$C2[is_treatment],
    max = par1["max"], slope1 = par1["slope1"], slope2 = par1["slope2"],
    ec501 = par1["ec501"], ec502 = par1["ec502"]
  )
  SSR_init1 <- sum((df$Res.[is_treatment] - preds_init1)^2)
  cat(sprintf("  Initial SSR (treatments): %.3f\n", SSR_init1))
  
  # ---- Objective function for Step 1 ----
  # The optimizer tries to find parameter values that minimize this function.
  # Here, it minimizes the SSR for treatment data only.
  obj_step1 <- function(theta) {
    
    # Replace the current slope and EC50 values with the optimizer's trial values.
    par1[opt_idx1] <- theta
    
    # Predict responses using the current parameter values.
    preds <- model_fun(
      C1 = df$C1[is_treatment], C2 = df$C2[is_treatment],
      max = par1["max"], slope1 = par1["slope1"], slope2 = par1["slope2"],
      ec501 = par1["ec501"], ec502 = par1["ec502"]
    )
    
    # Return SSR.
    # This is the value the optimizer tries to make as small as possible.
    sum((df$Res.[is_treatment] - preds)^2)
  }
  
  # ---- Run Step 1 optimization ----
  # NLOPT_LN_SBPLX is a derivative-free optimizer.
  # This can be useful when the model surface is irregular or difficult.
  # This is different from the optimizer function/calculation we used in Excel.
  # I was stuck here trying to find a way to find a comparable function or recreate the "GRG Non-linear" solving function used in Excel
  res1 <- nloptr::nloptr(
    x0 = theta0_1,
    eval_f = obj_step1,
    lb = lb1,
    ub = ub1,
    opts = list(
      algorithm = "NLOPT_LN_SBPLX",
      xtol_rel = 1e-6,
      maxeval = 500
    )
  )
  
  # Update the parameter vector with the optimized Step 1 values.
  par1[opt_idx1] <- res1$solution
  
  # ---- Print Step 1 optimized values ----
  cat("[Step 1] Optimized values (treatments):\n")
  cat(sprintf("  max    = %.4f (fixed)\n", par1["max"]))
  cat(sprintf("  slope1 = %.4f\n  slope2 = %.4f\n  ec501  = %.4f\n  ec502  = %.4f\n", 
              par1["slope1"], par1["slope2"], par1["ec501"], par1["ec502"]))
  
  # Calculate SSR after Step 1 optimization.
  preds_opt1 <- model_fun(
    C1 = df$C1[is_treatment], C2 = df$C2[is_treatment],
    max = par1["max"], slope1 = par1["slope1"], slope2 = par1["slope2"],
    ec501 = par1["ec501"], ec502 = par1["ec502"]
  )
  
  SSR_opt1 <- sum((df$Res.[is_treatment] - preds_opt1)^2)
  cat(sprintf("  Optimized SSR (treatments): %.3f\n\n", SSR_opt1))
  
  # ========================================================
  # Step 2: Fit max using all data
  # ========================================================
  # In this step:
  #   - max is optimized
  #   - slope1, slope2, ec501, and ec502 are fixed from Step 1
  #
  # All rows are used, including controls.
  
  par2         <- par1
  
  # Only max is optimized in Step 2.
  opt_idx2     <- idx_max
  
  # Starting value for max.
  theta0_2     <- par1[opt_idx2]
  
  # Bounds for max.
  lb2          <- pmax(1e-8, theta0_2 * 0.1)
  ub2          <- pmax(lb2 * 1.01, theta0_2 * 10)
  
  # ---- Print Step 2 starting values ----
  cat("[Step 2: Max Only] Initial values (all data):\n")
  cat(sprintf("  max    = %.4f\n", par2["max"]))
  cat(sprintf("  slope1 = %.4f (fixed)\n  slope2 = %.4f (fixed)\n  ec501  = %.4f (fixed)\n  ec502  = %.4f (fixed)\n", 
              par2["slope1"], par2["slope2"], par2["ec501"], par2["ec502"]))
  
  # Calculate SSR before optimizing max.
  preds_init2 <- model_fun(
    C1 = df$C1, C2 = df$C2,
    max = par2["max"], slope1 = par2["slope1"], slope2 = par2["slope2"],
    ec501 = par2["ec501"], ec502 = par2["ec502"]
  )
  SSR_init2 <- sum((df$Res. - preds_init2)^2)
  cat(sprintf("  Initial SSR (all data): %.3f\n", SSR_init2))
  
  # ---- Objective function for Step 2 ----
  # This time, only max is changed.
  obj_step2 <- function(theta) {
    par2[opt_idx2] <- theta
    preds <- model_fun(
      C1 = df$C1, C2 = df$C2,
      max = par2["max"], slope1 = par2["slope1"], slope2 = par2["slope2"],
      ec501 = par2["ec501"], ec502 = par2["ec502"]
    )
    sum((df$Res. - preds)^2)
  }
  
  # ---- Run Step 2 optimization ----
  res2 <- nloptr::nloptr(
    x0 = theta0_2,
    eval_f = obj_step2,
    lb = lb2,
    ub = ub2,
    opts = list(
      algorithm = "NLOPT_LN_SBPLX",
      xtol_rel = 1e-6,
      maxeval = 200
    )
  )
  
  # Update max with the optimized value.
  par2[opt_idx2] <- res2$solution
  
  # ---- Print final optimized values ----
  cat("[Step 2] Optimized values (all data):\n")
  cat(sprintf("  max    = %.4f\n", par2["max"]))
  cat(sprintf("  slope1 = %.4f\n  slope2 = %.4f\n  ec501  = %.4f\n  ec502  = %.4f\n", 
              par2["slope1"], par2["slope2"], par2["ec501"], par2["ec502"]))
  
  # Calculate final SSR using all data.
  preds_opt2 <- model_fun(
    C1 = df$C1, C2 = df$C2,
    max = par2["max"], slope1 = par2["slope1"], slope2 = par2["slope2"],
    ec501 = par2["ec501"], ec502 = par2["ec502"]
  )
  SSR_opt2 <- sum((df$Res. - preds_opt2)^2)
  cat(sprintf("  Optimized SSR (all data): %.3f\n\n", SSR_opt2))
  
  # ---- Return final fitted parameter values ----
  # The result is returned as a named vector.
  #
  # NOTE:
  # Later in the code, observeFitSingleModel() expects a list containing
  # params, preds, residuals, and ssr. This function currently only returns
  # the parameter vector. That mismatch needs to be fixed.
  return(par2)
}


# ============================================================
# Get manually entered parameter overrides
# ============================================================
# This helper function checks whether the user manually entered parameter values in the UI.
#
# If a value exists, it is returned as an override.
# If the value is missing, NULL is returned.
#
# Example:
#   prefix = "single"
#   param = "ec501"
#   bin_idx = "12"
#
# The function looks for:
#   input$single_ec501_12

get_param_overrides <- function(prefix, params, bin_idx, input) {
  sapply(params, function(p) {
    val <- input[[paste0(prefix, "_", p, "_", bin_idx)]]
    if (!is.null(val) && !is.na(val)) val else NULL
  }, simplify = FALSE)
}

# ============================================================
# Observe and fit the "Single Chemical" binary model
# ============================================================
# This function sets up an observeEvent for the "Fit Singles Model" button.
#
# When that button is clicked:
#   1. The binary mixture data are retrieved
#   2. Starting values are imported from the fitted single-chemical models
#   3. User-entered overrides are applied
#   4. Fixed/free parameter switches are read
#   5. The model is fitted
#   6. Results are stored in bin_values
#   7. The UI parameter boxes are updated with fitted values
observeFitSingleModel <- function(input, output, session, bin_idx, bin_values, sgl_values, model_fun) {
  observeEvent(input[[paste0("fit_single_model_", bin_idx)]], {
    
    # Make sure binary mixture data exist before fitting.
    req(bin_values[[bin_idx]]$data)
    
    df <- bin_values[[bin_idx]]$data
    
    # Print dataset information to the console for debugging.
    cat(paste0("[observeFitSingleModel] bin_idx=", bin_idx, " | df has ", nrow(df), " rows and cols: ", paste(colnames(df), collapse=", "), "\n"))
    
    # Get starting values from the corresponding fitted single-chemical models.
    initials <- get_initials_from_singles(bin_idx, sgl_values)
    
    # Read any manual parameter values entered by the user.
    overrides <- get_param_overrides("single", names(initials), bin_idx, input)
    
    # Replace imported starting values with user-entered values where available.
    for (param in names(overrides)) {
      if (!is.null(overrides[[param]])) {
        initials[[param]] <- overrides[[param]]
      }
    }
    
    # Read whether each parameter should be fixed during fitting.
    # TRUE means the parameter should be fixed.
    # FALSE means the parameter is allowed to change.
    fix_flags <- sapply(names(initials), function(param) {
      isTRUE(input[[paste0("fix_single_", param, "_", bin_idx)]])
    })
    
    # Fit the model.
    # NOTE:
    # The function receives model_fun as an argument, but currently calls CA_bi_vec directly instead of using model_fun. 
    fit_result <- fit_single_model(df, initials, fix_flags, CA_bi_vec)
    if (is.null(fit_result)) return()
    
    print("✅ Fit complete for bin_idx: ", quote = FALSE)

    # NOTE:
    # This currently assumes fit_result has a $params element.
    # But fit_single_model() currently returns only a named vector.
    # This will break unless fit_single_model() is changed.
    print(fit_result$params)
    showNotification(paste0("Single model fit completed for combo ", bin_idx), type = "message")
    
    
    # ---- Store fitted model results ----
    # NOTE:
    # These lines assume fit_result is a list with:
    #   fit_result$params
    #   fit_result$preds
    #   fit_result$residuals
    #   fit_result$ssr
    #
    # But the current fit_single_model() does not return that structure yet.
    bin_values[[bin_idx]]$single_model <- fit_result$params
    bin_values[[bin_idx]]$single_preds <- fit_result$preds
    bin_values[[bin_idx]]$single_residuals <- fit_result$residuals
    bin_values[[bin_idx]]$single_ssr <- fit_result$ssr
    
    # Update the parameter input boxes with the fitted parameter values.
    for (param in names(fit_result$params)) {
      updateNumericInput(
        session,
        inputId = paste0("single_", param, "_", bin_idx),
        value = round(fit_result$params[[param]], 4)
      )
    }
  })
}



# ============================================================
# Binary Mixture UI wrapper
# ============================================================
# This function creates the outer UI for the binary mixture module.
#
# The actual UI is generated inside the server using renderUI(),
# because it depends on the chemical names entered in the Introduction tab.
binaryMixtureModuleUI <- function(id = "binary") {
  ns <- NS(id)
  uiOutput(ns("binary_ui"))
}

# ============================================================
# Binary Mixture Server
# ============================================================
# This function controls the binary mixture module.
#
# It creates:
#   - Sub-tabs for the three binary combinations
#   - Parameter input matrices for different model types
#   - Raw data input and conversion
#   - Buttons for importing parameter values and fitting models
binaryMixtureModuleServer <- function(id, chem_names, sgl_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Binary mixture combinations.
    # Each value indicates which two chemicals are combined.
    bin_combos <- c("12", "13", "23")
    
    # Store binary mixture data and model results.
    #
    # NOTE:
    # This object is created inside the binary module.
    # If the main app also creates bin_values, this local version is separate
    # and will not be accessible outside this module.
    bin_values <- reactiveValues()
    
    # ========================================================
    # Dynamically create binary mixture sub-tabs
    # ========================================================
    # The UI is generated only after the three chemical names are available.
    # This allows the tab titles and labels to use the actual chemical names.
    output$binary_ui <- renderUI({
      req(chem_names$chem1, chem_names$chem2, chem_names$chem3)
      
      # Map binary mixture indices to chemical names.
      chem_map <- list(
        "12" = list(chem1 = chem_names$chem1, chem2 = chem_names$chem2),
        "13" = list(chem1 = chem_names$chem1, chem2 = chem_names$chem3),
        "23" = list(chem1 = chem_names$chem2, chem2 = chem_names$chem3)
      )
      
      # Create a tabsetPanel containing one tab per binary mixture.
      do.call(tabsetPanel, c(list(id = ns("binary_subtabs")),
                             lapply(bin_combos, function(bin_idx) {
                               chem1_name <- chem_map[[bin_idx]]$chem1
                               chem2_name <- chem_map[[bin_idx]]$chem2
                               
                               # ---- One binary mixture sub-tab ----
                               tabPanel(
                                 title = paste(chem1_name, "&", chem2_name),
                                 tagList(
                                   shinyjs::useShinyjs(),
                                   h3(paste("Binary Mixture of", chem1_name, "&", chem2_name)),
                                   
                                   # ---- CA / IA mode switch ----
                                   # This switch changes the label between CA and IA.
                                   #
                                   # NOTE:
                                   # At the moment, this only changes the labels.
                                   # The actual model function has not been linked to switch yet.
                                   fluidRow(
                                     column(12,
                                            switchInput(ns(paste0("mode_switch_", bin_idx)),
                                                        label = NULL,
                                                        onLabel = "CA",
                                                        offLabel = "IA",
                                                        value = TRUE))
                                   ),
                                   
                                   br(),
                                   
                                   # ==================================================
                                   # Parameter input matrix
                                   # ==================================================
                                   # This block creates the parameter table where the user can
                                   # enter, import, fix, or fit model parameters.
                                   #
                                   # Columns represent model types:
                                   #   - single = single-chemical-based binary model
                                   #   - ca     = CA or IA model, depending on the switch
                                   #   - sa     = S/A model
                                   #   - dr     = DR model
                                   #   - dl     = DL model
                                   #
                                   # Rows represent model parameters:
                                   #   - max
                                   #   - slope1
                                   #   - slope2
                                   #   - ec501
                                   #   - ec502
                                   #   - a
                                   #   - b
                                   div(
                                     style = "margin-bottom: 30px;",
                                     
                                     # ---- Header row of parameter matrix ----
                                     div(style = "display: flex; font-weight: bold; align-items: center; margin-bottom: 10px;",
                                         div(style = "width: 180px; font-size: 16px;", "Model Parameters"),
                                         
                                         # Create one header for each model type.
                                         lapply(c("single","ca", "sa", "dr", "dl"), function(prefix) {
                                           
                                           # The CA column label is dynamic.
                                           # It shows either CA or IA depending on the mode switch.
                                           if (prefix == "ca") {
                                             return(div(
                                               style = "width: 200px; display: flex; flex-direction: column; align-items: center; justify-content: center;",
                                               uiOutput(ns(paste0("ca_label_", bin_idx)))
                                             ))
                                           }
                                           
                                           # Static labels for the other model types.
                                           label <- switch(prefix,
                                                           single = "Single Chemical",
                                                           sa = "S/A",
                                                           dr = "DR",
                                                           dl = "DL")
                                           div(style = "width: 200px; display: flex; flex-direction: column; align-items: center; justify-content: center;", 
                                               div(style = "font-size: 16px; margin-left: -30px", label))
                                         })
                                     ),
                                     
                                     # ---- Parameter rows ----
                                     # For each parameter, create input boxes across the relevant model columns.
                                     lapply(
                                       list(
                                         list(label = "Max", param = "max"),
                                         list(label = paste0("Slope ", chem1_name), param = "slope1"),
                                         list(label = paste0("Slope ", chem2_name), param = "slope2"),
                                         list(label = paste0("EC50 ", chem1_name), param = "ec501"),
                                         list(label = paste0("EC50 ", chem2_name), param = "ec502"),
                                         list(label = "a", param = "a"),
                                         list(label = "b", param = "b")
                                       ),
                                       
                                       function(p) {
                                         div(style = "display: flex; align-items: center; margin-bottom: 10px;",
                                             
                                             # Parameter name on the left.
                                             div(style = "width: 180px; display: flex; align-items: center;", strong(p$label)),
                                             
                                             # Create input boxes for each model column.
                                             lapply(c("single", "ca", "sa", "dr", "dl"), function(prefix) {
                                               
                                               # Decide whether this parameter should appear
                                               # for the current model type.
                                               show_param <- switch(prefix,
                                                                    
                                                                    # Single and CA/IA models use only the basic max, slope, and EC50 parameters.
                                                                    single = p$param %in% c("max", "slope1", "slope2", "ec501", "ec502"),
                                                                    ca = p$param %in% c("max", "slope1", "slope2", "ec501", "ec502"),
                                                                    
                                                                    # S/A uses the basic parameters plus parameter a.
                                                                    sa = p$param %in% c("max", "slope1", "slope2", "ec501", "ec502", "a"),
                                                                    
                                                                    # DR and DL currently show all parameters.
                                                                    dr = TRUE,
                                                                    dl = TRUE
                                               )
                                               if (show_param) {
                                                 
                                                 # Numeric parameter input plus a fixed/free switch.
                                                 div(style = "width: 200px; display: flex; align-items: center; gap: 5px; justify-content: center;",
                                                     numericInput(ns(paste0(prefix, "_", p$param, "_", bin_idx)), label = NULL, value = NULL, width = "200px"),
                                                     switchInput(ns(paste0("fix_", prefix, "_", p$param, "_", bin_idx)), 
                                                                 label = NULL, 
                                                                 onLabel = "Fixed",
                                                                 offLabel = "",
                                                                 value = FALSE, size = "mini")
                                                 )
                                               } else {
                                                 
                                                 # Empty space used when a parameter does not apply
                                                 # to a certain model type.
                                                 div(style = "width: 200px;")
                                               }
                                             })
                                         )
                                       }
                                     )
                                   ),
                                   
                                   # ---- Model action buttons ----
                                   # Buttons are rendered separately so their labels can change
                                   # depending on the selected CA/IA mode.
                                   uiOutput(ns(paste0("model_buttons_", bin_idx))),
                                   
                                   # ==================================================
                                   # Binary mixture data input
                                   # ==================================================
                                   # The user pastes data with three columns:
                                   #   C1   = concentration of the first chemical
                                   #   C2   = concentration of the second chemical
                                   #   Res. = observed response
                                   fluidRow(
                                     column(12,
                                            br(),
                                            box(title = paste0(chem1_name, " & ", chem2_name, " Data Input"), status = "primary", solidHeader = TRUE, width = 12,
                                                div(class = "scroll-box",
                                                    helpText("Paste data below with three columns: 'C1', 'C2', and 'Res.'"),
                                                    div(style = "margin-bottom: 10px;",
                                                        actionButton(ns(paste0("toggle_data_", bin_idx)), "Show/Hide Raw Data")
                                                    ),
                                                    div(id = ns(paste0("paste_block_", bin_idx)),
                                                        textAreaInput(ns(paste0("paste_data_", bin_idx)), label = NULL, rows = 60),
                                                        actionButton(ns(paste0("convert_data_", bin_idx)), "Convert to Table")
                                                    ),
                                                    br(), br(),
                                                    rHandsontableOutput(ns(paste0("data_table_", bin_idx)))
                                                )
                                            )
                                     )
                                   )
                                 )
                               )
                             })
      ))
    })
    
    # ========================================================
    # Server logic for each binary combination
    # ========================================================
    # This loop creates server-side behaviour for each binary mixture tab.
    #
    # local({ ... }) is used so each binary mixture keeps its own bin_idx.
    lapply(bin_combos, function(idx) {
      local({
        bin_idx <- idx
        
        # ====================================================
        # Dynamic CA / IA label
        # ====================================================
        # The parameter matrix header shows either CA or IA depending on the mode switch.
        output[[paste0("ca_label_", bin_idx)]] <- renderUI({
          label <- if (isTruthy(input[[paste0("mode_switch_", bin_idx)]])) "CA" else "IA"
          div(style = "font-size: 16px; margin-left: -30px", label)
        })
        
        # ====================================================
        # Model button row
        # ====================================================
        # This creates the "Get Values" and "Fit Model" buttons underneath the parameter matrix.
        #
        # The labels change depending on the model type.
        output[[paste0("model_buttons_", bin_idx)]] <- renderUI({
          div(style = "display: flex; justify-content: flex-start; align-items: flex-start",
              
              # Empty left spacer to align buttons with the parameter columns.
              div(style = "width: 160px;", ""),
              
              lapply(c("single", "ca", "sa", "dr", "dl"), function(prefix) {
                
                # ---- Get-values button label ----
                # Not every model type has a get-values button yet.
                get_label <- if (prefix == "ca") {
                  "Get Singles Values"
                } else {
                  switch(prefix,
                         
                         # For S/A, the imported values come from either CA or IA, depending on the mode switch.
                         sa = {
                           mode_label <- if (isTruthy(input[[paste0("mode_switch_", bin_idx)]])) "CA" else "IA"
                           paste("Get", mode_label, "Values")
                         },
                         
                         # DR and DL both import S/A values.
                         dr = "Get S/A Values",
                         dl = "Get S/A Values",
                         NULL)
                }
                
                # ---- Fit-model button label ----
                fit_label <- if (prefix == "ca") {
                  mode_label <- if (isTruthy(input[[paste0("mode_switch_", bin_idx)]])) "CA" else "IA"
                  paste("Fit", mode_label, "Model")
                } else {
                  switch(prefix,
                         single = "Fit Singles Model",
                         sa = "Fit S/A Model",
                         dr = "Fit DR Model",
                         dl = "Fit DL Model")
                }
                
                # Create the actual buttons.
                btns <- list()
                if (!is.null(get_label)) {
                  btns <- append(btns, list(
                    actionButton(ns(paste0("get_", prefix, "_values_", bin_idx)), get_label)
                  ))
                } else {
                  btns <- append(btns, list(div(style = "min-height: 33px; display: inline-block;")))
                }
                
                # Placeholder space so button alignment stays consistent.
                btns <- append(btns, list(
                  actionButton(ns(paste0("fit_", prefix, "_model_", bin_idx)), fit_label)
                ))
                
                # One vertical button stack per model type.
                div(style = "width: 200px; display: flex; flex-direction: column; align-items: center; row-gap: 10px;", btns)
              })
          )
        })
        
        # ====================================================
        # Create ID names for the current binary mixture
        # ====================================================
        # These IDs make it easier to refer to the correct inputs
        # and outputs for each binary combination.
        bin_idx <- idx
        
        get_single <- paste0("get_single_values_", bin_idx)
        get_ca     <- paste0("get_ca_values_", bin_idx)
        get_sa     <- paste0("get_sa_values_", bin_idx)
        get_dr     <- paste0("get_dr_values_", bin_idx)
        get_dl     <- paste0("get_dl_values_", bin_idx)
        
        fit_single <- paste0("fit_single_model_", bin_idx)
        fit_ca     <- paste0("fit_ca_model_", bin_idx)
        fit_sa     <- paste0("fit_sa_model_", bin_idx)
        fit_dr     <- paste0("fit_dr_model_", bin_idx)
        fit_dl     <- paste0("fit_dl_model_", bin_idx)
        
        toggle <- paste0("toggle_data_", bin_idx)
        convert <- paste0("convert_data_", bin_idx)
        paste_id <- paste0("paste_data_", bin_idx)
        table <- paste0("data_table_", bin_idx)
        paste_block <- paste0("paste_block_", bin_idx)
        output_toggle <- paste0("show_table_", bin_idx)
        
        # This reactive output checks whether data exist for the current binary mixture.
        #
        # NOTE:
        # At the moment, output_toggle does not appear to be used in the UI.
        output[[output_toggle]] <- reactive({
          !is.null(bin_values[[bin_idx]]$data)
        })
        outputOptions(output, output_toggle, suspendWhenHidden = FALSE)
        
        # ====================================================
        # Show or hide raw data paste area
        # ====================================================
        # When the user clicks the toggle button, the raw data input area is shown or hidden.
        observeEvent(input[[toggle]], {
          shinyjs::toggle(paste_block)
        })
        
        # ====================================================
        # Convert pasted binary mixture data into a table
        # ====================================================
        # When the user clicks "Convert to Table":
        #   1. The pasted text is read as a data frame
        #   2. The data are stored in bin_values
        #   3. The raw paste area is hidden
        #   4. The data are displayed as an editable table
        observeEvent(input[[convert]], {
          req(input[[paste_id]])
          
          # Read tab-separated data.
          # Expected columns:
          #   C1
          #   C2
          #   Res.
          df <- tryCatch({
            read.table(text = input[[paste_id]],
                       header = TRUE,
                       sep = "\t",
                       strip.white = TRUE,
                       stringsAsFactors = FALSE,
                       fill = TRUE,
                       comment.char = "")
          })
          
          # Store the binary mixture data for this combination.
          bin_values[[bin_idx]]$data <- df
          
          # Hide the raw paste area after conversion.
          shinyjs::hide(paste_block)
          
          # Display the converted data as an editable table.
          output[[table]] <- renderRHandsontable({
            rhandsontable(df)
          })
          
          # Set up the observer for fitting the Single Chemical model.
          #
          # NOTE:
          # This observer is currently created only after data conversion.
          # If the user converts data multiple times, this can create multiple observers for the same fit button.
          # It would be cleaner to move observeFitSingleModel() outside this observeEvent.
          observeFitSingleModel(
            input = input,
            output = output,
            session = session,
            bin_idx = bin_idx,
            bin_values = bin_values,
            sgl_values = sgl_values,
            model_fun = CA_bi_vec
          )
        })
      })
    })
  })
}