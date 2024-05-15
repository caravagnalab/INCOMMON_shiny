#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize = 30*1024^2)
library(shiny)
library(ggplot2)
library(waiter)
# library(INCOMMON)
dir_scripts <- "~/Documents/GitHub/INCOMMON/R/"
files.sources = list.files(path = dir_scripts,full.names = T)
sapply(files.sources, source)

# Define server logic required to draw a histogram
server = function(input, output, session) {

  # What happens when clicking on Classification button
  observeEvent(input$classification_button, {
    output$selected_module <- renderUI({
      # Call the classification shiny app module
      classification_ui("classification")
    })
    callModule(classification_module, "classification")
  })

  # What happens when clicking on Survival Analysis button

  observeEvent(input$survival_button, {
    output$selected_module <- renderUI({
      # Call the survival analysis shiny app module
      survival_analysis_ui("survival")
    })
    callModule(survival_analysis_module, "survival")
  })


  # What happens when clicking on Metastatic propensity button
  observeEvent(input$met_propensity_button, {
    output$selected_module <- renderUI({
      # Call the survival analysis shiny app module
      met_propensity_ui("met_propensity")
    })
    callModule(met_propensity_module, "met_propensity")
  })

  # What happens when clicking on Metastatic propensity button
  observeEvent(input$met_tropism_button, {
    output$selected_module <- renderUI({
      # Call the survival analysis shiny app module
      met_tropism_ui("met_tropism")
    })
    callModule(met_tropism_module, "met_tropism")
  })
}
