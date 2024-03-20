library(shiny)

# Load your shiny app modules
source("modules/shiny_module_runFit.R")
source("modules/shiny_module_Survival.R")

# Define UI for shiny app
ui <- fluidPage(
  titlePanel("Choose Analysis"),
  
  # Action buttons to choose analysis
  actionButton("survival_button", "Analysis Survival"),
  actionButton("fit_button", "Analysis Fit"),
  
  # Placeholder for rendering the selected shiny app module
  uiOutput("selected_module")
)

# Define server logic
server <- function(input, output, session) {
  
  # Render the selected shiny app module based on button click
  observeEvent(input$survival_button, {
    output$selected_module <- renderUI({
      # Call the survival analysis shiny app module
      survivalModuleUI("survival")
    })
    callModule(survivalModule, "survival")
  })
  
  observeEvent(input$fit_button, {
    output$selected_module <- renderUI({
      # Call the INCOMMON analysis shiny app module
      incommonModuleUI("incommon")
    })
    callModule(incommonModule, "incommon")
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
