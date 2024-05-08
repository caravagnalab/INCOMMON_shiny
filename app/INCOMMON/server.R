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
library(dplyr)
library(ggplot2)

# library(INCOMMON)
dir_scripts <- "~/Documents/GitHub/INCOMMON/R/"
files.sources = list.files(path = dir_scripts,full.names = T)
sapply(files.sources, source)

# Define server logic required to draw a histogram
server = function(input, output, session) {

  observeEvent(input$survival_button, {
    output$selected_module <- renderUI({
      survival_panel
    })
  })

}
