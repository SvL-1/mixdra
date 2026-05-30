# ============================================================
# Mixture Toxicity Dashboard - Main App File
# ============================================================
# This script is the starting point of the Shiny app.
# It loads all required packages, imports helper functions and modules,
# defines the dashboard layout, connects the server logic,
# and finally launches the app.
# ============================================================


library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(rhandsontable)
library(drc)
library(shinyjs)
library(shinyWidgets)
library(numDeriv)
library(Rsolnp)
library(nloptr)

# ---- Load external scripts ----
# These source files contain the model functions and the separate app modules.

source("functions/model_functions.R", local = TRUE)        # Contains the binary and ternary model functions 
source("modules/intro_module.R", local = TRUE)             # Introduction page module
source("modules/single_module.R")                          # Single chemical dose-response curve fitting module
source("modules/binary_module.R", local = TRUE)            # Binary mixture analysis module


# ============================================================
# User Interface
# ============================================================
# The UI controls what the user sees in the app.
# Here, the app is arranged as a dashboard with a header, sidebar menu,
# and main body area containing different tabs.


ui <- dashboardPage(
  
  # ---- Dashboard title ---- #
  dashboardHeader(title = "Mixture Toxicity Dashboard"),                               
  
  # ---- Sidebar menu ---- #
  # The sidebar contains the main navigation tabs of the app.
  # Each menu item links to a different part of the dashboard.
  # The users should start with the introduction, then to the single chemical, binary and then ternary mixture tabs,
  # As the previous model fit parameters will be used as the start point of the following model
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("Introduction", tabName = "intro", icon = icon("info-circle")),
                menuItem("Single Chemical", tabName = "single", icon = icon("vial")),
                menuItem("Binary Mixture", tabName = "binary", icon = icon("vials")),
                menuItem("Ternary Mixture", tabName = "ternary", icon = icon("flask"))
    )
  ),
  
  # ---- Dashboard body ---- #
  # This is the main content area where the selected tab is displayed.
  dashboardBody(
    
    # Enables shinyjs functions used elsewhere in the app
    useShinyjs(),
    
    # ---- Custom CSS styling ----
    # This adds a reusable scrollable box style.
    # Any UI element with class "scroll-box" will have a maximum height
    # and a vertical scrollbar when the content becomes too long.
    tags$head(
      tags$style(HTML(
        ".scroll-box { overflow-y: auto; max-height: 700px; }"
      ))
    ),
    
    # ---- Main tab contents ----
    # Each tabItem defines what should appear when a sidebar tab is selected.
    # Some tabs use module UI functions directly, while others are rendered dynamically.
    tabItems(
      tabItem(tabName = "intro", introUI("intro")),
      tabItem(tabName = "single", uiOutput("single_ui")),
      tabItem(tabName = "binary", uiOutput("binary_ui"))
      
      
      # Note:
      # The sidebar includes a "Ternary Mixture" tab, but no ternary tabItem is defined yet.
      # So now clicking on the Ternary tab will currently show an empty page.
      # My idea was to add a Summary tab, where the plots, model parameters can be retrieved,
      # and stored/saved as an output file. But this hasn't been done yet
    )
  )
)

# ============================================================
# Server
# ============================================================
# The server controls how the app behaves.
# It stores shared values, connects module logic, reacts to user input,
# and generates outputs such as plots, tables, and model results.

server <- function(input, output, session) {
  
  # ---- Create shared storage objects ----
  # These reactiveValues objects store information that needs to be shared across different modules.
  # And the reactiveValues should allow model parameter retrieve at the end
  
  # The following sgl_values stores results from the single chemical module, such as fitted models,
  # predictions, residuals, plots, and summary tables.
  sgl_values <- reactiveValues()
  
  # The following bin_values intended to store results from the binary mixture module.
  # At the moment, this object is created here but is not passed into
  # binaryMixtureModuleServer() below.
  bin_values <- reactiveValues()
  
  # chem_names stores general experiment information entered in the Introduction tab.
  # These values can then be reused by other modules.
  chem_names <- reactiveValues(chem1 = NULL, chem2 = NULL, chem3 = NULL, unit = NULL, expID = NULL, species = NULL, endpoint = NULL)
  
  # ---- Introduction module ----
  # This module collects experiment information and chemical names.
  # The information is stored in chem_names.
  introServer("intro", chem_names)
  
  # ---- Single chemical module ----
  # This module handles dose-response analysis for individual chemicals.
  # It uses the chemical names from chem_names and stores results in sgl_values.
  singleChemicalServer("single", chem_names, sgl_values)
  
  # ---- Render single chemical UI ----
  # The UI for the single chemical module is generated dynamically.
  output$single_ui <- renderUI({
    singleChemicalUI("single")
  })
  
  # ---- Binary mixture module ----
  # This module handles binary mixture analysis.
  # It uses chemical names and can also reuse fitted single-chemical results (as starting point),
  # from sgl_values as starting values for mixture model fitting.
  binaryMixtureModuleServer("binary", chem_names, sgl_values)
  
  # ---- Render binary mixture UI ----
  # The UI for the binary mixture module is generated dynamically.
  output$binary_ui <- renderUI({
  binaryMixtureModuleUI("binary")
  })
}

# ---- Launch the app ----
# This line starts the Shiny application using the UI and server defined above.
shinyApp(ui, server)

