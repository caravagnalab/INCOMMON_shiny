library(dplyr)
library(DT)

met_propensity_ui = function(id) {
  ns <- NS(id)
  fluidPage(
    tabsetPanel(
      id = ns('tabs'),
      tabPanel(
        title = 'Metastatic propensity',
        fluidRow(
          sidebarLayout(
            sidebarPanel(
              fileInput(ns("dataFile"), label = HTML(paste0("Upload your own INCOMMON object with classified data\n",
                                                            "or use results of our analysis available at <a href='https://zenodo.org/records/10927218'>zenodo.org/records/10927218 </a>",
                                                            "\n('msk_classified_with_priors.rds')")), accept = c(".rds")),
              selectizeInput(ns("tumorType"), "Tumor Type:", choices = NULL),
              actionButton(ns("plotButton"), "Plot"),
              downloadButton(ns("downloadPlot"), "Download Plot"),
              width = 4
            ),
            
            mainPanel(
              plotOutput(ns("met_prop_plot"),width = "80%", height = "400px")
            )
          )
        )
      )
    )
    # titlePanel("Survival Analysis"),
    #
    # sidebarLayout(
    #   sidebarPanel(
    # fileInput(ns("dataFile"), label = HTML(paste0("Upload your own INCOMMON object with classified data\n",
    #                                           "or use results of our analysis available at <a href='https://zenodo.org/records/10927218'>zenodo.org/records/10927218 </a>",
    #                                           "\n('msk_classified_with_priors.rds')")), accept = c(".rds")),
    #     checkboxGroupInput(ns("covariates"), "Covariates:", choices = choices_covariates),
    #     selectizeInput(ns("tumorType"), "Tumor Type:", choices = NULL),
    #     selectizeInput(ns("gene"), "Gene:", choices = NULL),
    #     actionButton(ns("plotButton"), "Plot"),
    #     downloadButton(ns("downloadPlot"), "Download Plot"),
    #     width = 4
    #   ),
    #
    # mainPanel(
    #   plotOutput(ns("survival_plot"),width = "100%", height = "600px")
    # )
    # )
  )
  
#   ns <- NS(id)
#   fluidPage(
#     tabsetPanel(
#       id = ns('tabs'),
#       tabPanel(
#         title = 'Metastatic propensity',
#         fluidRow(
#           column(
#             3,
#             h1("Metastatic propensity score"),
#             p("Infer the metastatic propensity of primary tumor genomes containing a specific gene"),
#             fileInput(
#               ns("dataFile"),
#               label = HTML(paste0("Upload your own initialised INCOMMON object\n",
#                                   "or use results of our analysis available at <a href='https://zenodo.org/records/10927218'>zenodo.org/records/10927218 </a>",
#                                   "\n('msk_classified_with_priors.rds')")
#               ),
#               accept = c(".rds")
#             ),
#             # p("Need purity column"),
#             selectizeInput(ns("tumor_type"), "Tumor Type:", choices = NULL),
#             actionButton(ns("plotButton"), "Plot")
#             # selectInput(ns("gene"), "Gene", choices = NULL),
#         ),
#         mainPanel(
#           plotOutput(ns("met_prop_plot"),width = "100%", height = "600px")
#         ),
#         downloadButton(ns("downloadPlot"), "Download .pdf")
#       )
#     )
#   )
# )
}

met_propensity_module = function(input, output, session) {

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
  # observe({
  #   req(data())
  #   updateSelectizeInput(session, "gene",
  #                        choices = unique(data()$genomic_data$gene),
  #                        server = TRUE)
  # })
  
  observeEvent(input$plotButton,{
    # Input validation
    req(input$tumorType, data())
    
    # Use INCOMMON functions to perform survival analysis and make summary figure
    do_figure = function(x, tumor_type){
      top_genes = classification(x) %>% 
        dplyr::filter(state != 'Tier-2') %>% 
        dplyr::group_by(gene) %>% 
        dplyr::reframe(N = length(unique(sample))) %>% 
        dplyr::arrange(dplyr::desc(N)) %>% 
        dplyr::slice_head(n = 50) %>% 
        pull(gene)
      # run metastatic propensity analysis
      for(g in top_genes){
        x = met_propensity(x, tumor_type = tumor_type, gene = g)
      }
      
      # Plot surv analysis
      plot <- plot_met_volcano(x = x, tumor_type = tumor_type)
      
      return(plot)
    }
    
    met_prop_plot = reactive({
      if (!is.null(data())) {
        do_figure(data(), input$tumorType)
      }
    })
    
    # Render plot
    output$met_prop_plot <- renderPlot({
      plots <- met_prop_plot()
      plots  # Return the plot object directly
    })
  })
  
  # Download plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("Metastatic_Prop_Plot-", input$tumorType, '.pdf', sep = "")
    },
    content = function(file) {
      ggsave(
        file, plot = last_plot(), device = "pdf",
        width = 12, height = 8, units = "in",
        dpi = 300
      )
    }
  )
  
}