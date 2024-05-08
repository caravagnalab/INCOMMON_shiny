#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
choices_covariates <- c("PRIMARY_SITE","SEX","AGE_AT_DEATH","AGE_AT_SEQUENCING","TMB_NONSYNONYMOUS","FGA")

ui <- fluidPage(
  # Title
  div(style = "background-color: #A9A9A7; padding: 10px;",
      titlePanel("Welcome to INCOMMON ShinyApp", windowTitle = "INCOMMON")),
  # Subtitle asking the user what they want to do
  div(style = "background-color: #B0BEC5; padding: 10px;",
      h3("Use INCOMMON for:")),

  # # Action buttons to choose what to do
  # actionButton("classification_button", "Classification", ),
  # actionButton("survival_button", "Survival Analysis"),
  # actionButton("met_propensity_button", "Metastatic Porpensity"),
  # actionButton("met_tropism_button", "Metastatic Tropism"),

  # CSS to style the action buttons
  tags$head(
    tags$style(HTML("
      .action-btn {
        width: 200px;
        height: 100px;
        margin: 20px;
        font-size: 18px;
        border-radius: 10px;
        position: relative;
      }
      .action-btn .tooltiptext {
        visibility: hidden;
        width: 220px;
        background-color: #555;
        color: #fff;
        text-align: center;
        border-radius: 6px;
        padding: 10px;
        position: absolute;
        z-index: 1;
        bottom: 125%;
        left: 50%;
        margin-left: -110px;
        opacity: 0;
        transition: opacity 0.3s;
      }
      .action-btn:hover .tooltiptext {
        visibility: visible;
        opacity: 1;
      }

      .classification-btn { background-color: #4CAF50; }
      .survival-btn { background-color: #2196F3; }
      .met-propensity-btn { background-color: #f44336; }
      .met-tropism-btn { background-color: #ff9800; }
    "))
  ),

  # Action buttons in two columns
  fluidRow(
    column(6, align = "center",
           actionButton("classification_button", "Classification",
                        class = "action-btn classification-btn",
                        title = "Input: Data for classification.\nOutput: Predicted classes for the data.\n")),
    column(6, align = "center",
           actionButton("survival_button", "Survival Analysis", class = "action-btn survival-btn",
                        title = "Input: Survival data with covariates.\nOutput: Survival curves and hazard estimates.")))
  ,
  fluidRow(
    column(6, align = "center",
           actionButton("met_propensity_button", "Metastatic Propensity", class = "action-btn met-propensity-btn",
                        title = "Input: Clinical and biological data.\nOutput: Predicted propensity of metastasis.")),
    column(6, align = "center",
           actionButton("met_tropism_button", "Metastatic Tropism", class = "action-btn met-tropism-btn",
                        title = "Input: Tumor imaging data.\nOutput: Prediction of metastatic tropism.")))
  ,

  # Placeholder for rendering the selected shiny app module
  uiOutput("selected_module")
)

survival_panel <- fluidPage(

  titlePanel("Survival Analysis"),

  data <- reactive({
    req(input$dataFile)
    readRDS(input$dataFile$datapath)
  }),

  # Choose tumour type
  observe({
    req(data())

    updateSelectizeInput(session, "tumorType",
                         choices = unique(data()$clinical_data$tumor_type),
                         server = TRUE)
  }),

  # Choose gene
  observe({
    req(data())

    updateSelectizeInput(session, "gene",
                         choices = unique(data()$genomic_data$gene),
                         server = TRUE)
  }),


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
  }),

  # Render plot
  output$survivalPlot <- renderPlot({
    plots <- survival_plot()
    plots  # Return the plot object directly
  }),

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("KM_Plot-", input$tumorType, '-', input$gene, '.pdf', sep = "")
    },
    content = function(file) {

      ggsave(file, plot = last_plot(), device = "pdf", width = 10, height = 8, units = "in", dpi = 300)
    }
  ),

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
