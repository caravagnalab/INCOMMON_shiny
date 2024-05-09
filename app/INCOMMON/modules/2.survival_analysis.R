# UI for the survival analysis module
survival_analysis_ui = function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Survival Analysis"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("dataFile"), "Choose a Data File (RDS):", accept = c(".rds")),
        selectizeInput(ns("tumorType"), "Tumor Type", choices = NULL),
        selectizeInput(ns("gene"), "Gene", choices = NULL),
        actionButton(ns("plotButton"), "Plot"),
        downloadButton(ns("downloadPlot"), "Download Plot"),
        width = 3
      ),
      mainPanel(
        plotOutput(ns("kmPlot"), width = "70%", height = "800px")
      )
    )
  )
}

# Server logic for the survival analysis module
survival_analysis_module = function(input, output, session) {
  # Read data from file
  data <- reactive({
    req(input$dataFile)
    readRDS(input$dataFile$datapath)
  })

  # Choose tumour type from drop-down menu
  observe({
    req(data())
    updateSelectizeInput(session, "tumorType",
                         choices = unique(data()$clinical_data$tumor_type),
                         server = TRUE)
  })

  # Choose gene from drop-down menu
  observe({
    req(data())
    updateSelectizeInput(session, "gene",
                         choices = unique(data()$genomic_data$gene),
                         server = TRUE)
  })


  observeEvent(input$plotButton,{
    # Input validation
    req(input$tumorType, input$gene, data())

    # Use INCOMMON functions to perform survival analysis and make summary figure
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
  })


}


# Run the application
shinyApp(
  ui = survival_analysis_ui("incommon"),
  server = function(input, output, session) {
    callModule(survival_analysis_module, "incommon")
  }
)
