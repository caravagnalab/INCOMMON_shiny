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
function(input, output, session) {
  
  data <- reactive({
    req(input$dataFile)
    readRDS(input$dataFile$datapath)
  })

  observe({
    req(data())

    updateSelectizeInput(session, "tumorType", choices = unique(data()$clinical_data$tumor_type), server = TRUE)
  })
  
  observe({
    req(data())

    updateSelectizeInput(session, "gene", choices = unique(data()$genomic_data$gene), server = TRUE)
  })
  

  # Define reactive expression for survival analysis
  
  
  do_figure = function(x, tumor_type, gene, covariates){
    
    x = kaplan_meier_fit(
      x = x,
      tumor_type = tumor_type,
      gene = gene,
      survival_time = 'OS_MONTHS',
      survival_status = 'OS_STATUS'
    )
    
    # Fit Cox
    x = cox_fit(
      x = x,
      tumor_type = tumor_type,
      gene = gene,
      survival_time = 'OS_MONTHS',
      survival_status = 'OS_STATUS',
      covariates = covariates,
      tmb_method = ">10"
    )
    
    # Plot surv analysis
    plot = plot_survival_analysis(
      x = x,
      tumor_type = tumor_type,
      gene = gene,
      cox_covariates = covariates
    )
    
    return(plot)
  }

  survival_plot <- reactive({
    if (!is.null(data())) {
      do_figure(data(), input$tumorType, input$gene, input$covariates)
    }
  })
  
  # Render plot
  output$survivalPlot <- renderPlot({
    plots <- survival_plot()
    plots  # Return the plot object directly
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("KM_Plot-", input$tumorType, '-', input$gene, '.pdf', sep = "")
    },
    content = function(file) {
      
      ggsave(file, plot = last_plot(), device = "pdf", width = 10, height = 8, units = "in", dpi = 300)
    }
  )
}
