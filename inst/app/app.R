# Launcher for the mixdra Shiny app. Assumes the mixdra package is installed.
library(mixdra)
shiny::shinyApp(ui = mixdra:::app_ui(), server = mixdra:::app_server)
