library(dplyr)
library(DT)

met_propensity_ui = function(id) {
  ns <- NS(id)
  fluidPage(
    # autoWaiter(),
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
              downloadButton(ns("downloadTable"), "Download Table"),
              width = 4
            ),
            
            mainPanel(
              plotOutput(ns("met_prop_plot"),width = "80%", height = "400px"),
              DTOutput(ns("met_prop_table"))
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
      results.names <- c("table", "plot") 
      results <-vector("list", length(results.names)) |> setNames(results.names)
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
        x$metastatic_propensity[tumor_type][[g]]
      }
      results$table <- do.call(rbind, x$metastatic_propensity[[tumor_type]])
      # Plot surv analysis
      results$plot <- plot_met_volcano(x = x, tumor_type = tumor_type)
      
      return(results)
    }
    
    met_prop_plot = reactive({
      if (!is.null(data())) {
        do_figure(data(), input$tumorType)
      }
    })
    
    # Render plot
    output$met_prop_plot <- renderPlot({
      plots <- met_prop_plot()$plot
      plots  # Return the plot object directly
    })
    
    # Render the output table of the fit
    output$met_prop_table <- renderDT({
      datatable((met_prop_plot()$table %>% 
                  dplyr::mutate(low = round(low, 4),
                                up = round(up, 4),
                                p.value = round(p.value, 4),
                                OR = round(OR,4))),
        options = list(scrollX = TRUE, scrollY = TRUE))
    })
  })
  
  # Download button for the output table
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("Metastatic_Prop_Results", input$tumorType,".csv", sep = "")
    },
    content = function(file) {
      write.csv(met_prop_plot()$table, file, row.names = FALSE,
                append = F, quote = F)
    }
  )
  
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