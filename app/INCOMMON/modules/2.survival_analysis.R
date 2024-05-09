library(INCOMMON)
choices_covariates <- c("PRIMARY_SITE",
                        "SEX",
                        "AGE_AT_DEATH",
                        "AGE_AT_SEQUENCING",
                        "TMB_NONSYNONYMOUS",
                        "FGA")
# UI for the survival analysis module
survival_analysis_ui = function(id) {
  ns <- NS(id)
  fluidPage(

    titlePanel("Survival Analysis"),

    sidebarLayout(
      sidebarPanel(
        fileInput("dataFile", label = HTML(paste0("Upload your own INCOMMON object with classified data\n",
                                                  "or use results of our analysis available at <a href='https://zenodo.org/records/10927218'>zenodo.org/records/10927218 </a>",
                                                  "\n('msk_classified_with_priors.rds')")), accept = c(".rds")),
        checkboxGroupInput("covariates", "Covariates:", choices = choices_covariates),
        selectizeInput("tumorType", "Tumor Type:", choices = NULL),
        selectizeInput("gene", "Gene:", choices = NULL),
        actionButton("plotButton", "Plot"),
        downloadButton("downloadPlot", "Download Plot"),
        width = 4
      ),

      mainPanel(
        plotOutput("survivalPlot",width = "100%", height = "600px")
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

    survival_plot = reactive({
      if (!is.null(data())) {
        do_figure(data(), input$tumorType, input$gene, input$covariates)
      }
    })

    # Render plot
    output$survival_plot <- renderPlot({
      plots <- survival_plot()
      plots  # Return the plot object directly
    })
  })

  # Download plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("KM_Plot-", input$tumorType, '-', input$gene, '.pdf', sep = "")
    },
    content = function(file) {
      ggsave(
        file, plot = last_plot(), device = "pdf",
        width = 8, height = 12, units = "in",
        dpi = 300
      )
    }
  )
}

# Run the application
shinyApp(
  ui = survival_analysis_ui("incommon"),
  server = function(input, output, session) {
    callModule(survival_analysis_module, "incommon")
  }
)
