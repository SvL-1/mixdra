# ============================================================
# Introduction Module
# ============================================================
# This module creates the Introduction page of the dashboard.
# It allows the user to enter basic experiment information,
# including the experiment ID, species, endpoint, chemical names,
# and concentration unit.
#
# The entered information is stored in chem_names, so it can be
# reused globally throughout the app.
# ============================================================


# ============================================================
# Introduction UI
# ============================================================
# This function defines what the Introduction page looks like.
# The page contains:
#   1. A welcome/instruction box on the left, for which text will be added
#   2. A chemical and experiment information input box on the right

introUI <- function(id = "intro") {
  
  # Create a namespace for this module.
  # This prevents input IDs in this module from accidentally conflicting
  # with input IDs in other modules.
  ns <- NS(id)
  
  # ---- Introduction tab layout ----
  # This tab is displayed when the user selects the Introduction page.
  tabItem(tabName = id,
          fluidRow(
            
            # ---- Left column: welcome/instruction text ----
            # This area is reserved for general instructions or background text.
            column(8,
                   box(title = "Welcome", status = "info", width = 12, solidHeader = TRUE,
                       # The scroll-box class makes this area scrollable if the text
                       # becomes longer than the available space.
                       div(class = "scroll-box", style = "height: 500px;",
                           # Placeholder text for now.
                           # The actual welcome or instruction text will be added later.
                           p("[Text to be added later]")
                       )
                   )
            ),
            
            # ---- Right column: experiment and chemical input ----
            # The user enters general information that will be used across the app.
            column(4,
                   box(title = "Chemical Input", status = "primary", width = 12, solidHeader = TRUE,
                       textInput(ns("exp_name"), "Experiment ID:"),
                       textInput(ns("species"), "Species:"),
                       textInput(ns("endpoint"), "Endpoint:"),
                       helpText("Please enter the chemical names in the correct order. This order will be used throughout the app."),
                       
                       # Chemical names entered by the user
                       textInput(ns("chem_label1"), "Name of Chemical 1:"),
                       textInput(ns("chem_label2"), "Name of Chemical 2:"),
                       textInput(ns("chem_label3"), "Name of Chemical 3:"),
                       
                       # Concentration unit used in the experiment
                       textInput(ns("chem_unit"), "Unit (e.g., mg/kg or mg/L)"),
                       br(),
                       
                       # ---- Navigation button ----
                       # This button sends the user to the Single Chemical tab.
                       div(style = "text-align: right;",
                           actionButton(ns("go_to_single"), "Next", icon = icon("arrow-right"))
                       )
                   )
            )
          )
  )
}


# ============================================================
# Introduction Server
# ============================================================
# This function controls the behaviour of the Introduction page.
# It does two main things:
#   1. Moves the user to the Single Chemical tab when they click "Next"
#   2. Saves the entered experiment and chemical information in chem_names

introServer <- function(id, chem_names) {
  moduleServer(id, function(input, output, session) {
    
    # ---- Move to Single Chemical tab ----
    # When the user clicks the "Next" button, the app switches from the
    # Introduction tab to the Single Chemical tab.
    
    observeEvent(input$go_to_single, {
      updateTabItems(session = session$rootScope(), inputId = "tabs", selected = "single")
    })
    
    # ---- Store experiment and chemical information ----
    # This observe block keeps chem_names updated with the values entered by the user.
    #
    # Because chem_names is a reactiveValues object created in the main app,
    # the stored values can be accessed by other modules.
    observe({
      
      # General experiment information
      chem_names$expID <- input$exp_name
      chem_names$species <- input$species
      chem_names$endpoint <- input$endpoint
      
      # Chemical names
      chem_names$chem1 <- input$chem_label1
      chem_names$chem2 <- input$chem_label2
      chem_names$chem3 <- input$chem_label3
      
      # Concentration unit
      chem_names$unit  <- input$chem_unit
    })
  })
}
