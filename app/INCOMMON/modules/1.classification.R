library(dplyr)
library(DT)

classification_ui = function(id) {
  ns <- NS(id)
  fluidPage(
    # autoWaiter(),
    tabsetPanel(
      id = ns('tabs'),
      tabPanel(
        title = 'Classification',
        fluidRow(
              column(
                3,
                h1("INCOMMON classification"),
                     p("Infer mutation copy number and multiplicity"),
                     fileInput(
                       ns("dataFile"),
                       label = HTML(paste0("Upload your own initialised INCOMMON object\n",
                                           "or use results of our analysis available at <a href='https://zenodo.org/records/10927218'>zenodo.org/records/10927218 </a>",
                                           "\n('msk_classified_with_priors.rds')")
                                    ),
                               accept = c(".rds")
                               ),
                     # p("Need purity column"),
                     selectInput(ns("sample"), "Sample", choices = NULL),
                     p("Select one sample"),
                     numericInput(ns("entropy"),
                                  label = "Entropy Cut-off",
                                  value = 1000),
                     p("Choose an entropy cut-off"),
                     numericInput(ns("rho"), label = "Rho", value = 0.01),
                     p("Over-dispersion parameter"),
                     textInput(ns("tumor_type"), label = "Tumor type", value = NA),
                     p("Tumor type"),
                     actionButton(ns("submit"), "Submit")),
              column(6, tableOutput(ns("data"))),
              column(6, plotOutput(ns("classification_plot")))
              ),
        downloadButton(ns("downloadPlot"), "Download .pdf"),
        downloadButton(ns("downloadRDS"), "Download .rds")
            # tabPanel(
            #   title = 'Output table',
            #   fluidRow(
            #     column(6, DTOutput(ns("samplesTable"))),
            #     column(6, plotOutput(ns("plot")))
            #   ),
            #   downloadButton(ns("downloadTable"), "Download .tsv"),
            #   downloadButton(ns("downloadPlot"), "Download .pdf")
            # )
        )
      )
    )
}

# Server logic for the survival analysis module
classification_module = function(input, output, session) {
  # Read data from file
  data <- reactive({
    req(input$dataFile)
    readRDS(input$dataFile$datapath)
  })

  # Select a single sample
  observe({
    updateSelectInput(session, "sample", choices = unique(data()$input$sample))
  })

  # Subset the original dataframe according to the select sample
  filtered_data <- reactive({
    req(input$sample)
    data() %>%
      INCOMMON:::subset_sample(x = ., sample = input$sample) %>%
      INCOMMON:::input() %>%
      dplyr::select(tumor_type, purity, gene, gene_role, chr, from, to, ref, alt, HGVSp_Short, NV, DP, VAF)
    # data()$input %>%
    #   filter(sample == input$sample) %>%
    #   select(chr, from, to, ref, alt, gene, NV, DP, VAF, purity)
  })
  
  filtered_data_incommon <- reactive({
    req(input$submit, input$sample)
    data() %>%
      INCOMMON:::subset_sample(x = ., sample = input$sample)
    # data()$input %>%
    #   filter(sample == input$sample) %>%
    #   select(chr, from, to, ref, alt, gene, NV, DP, VAF, purity)
  })
  
  output$data <- renderTable({
    req(input$sample)
    filtered_data()
  })
  
  # Step2 : run the fit
  # 
  do_figure = function(x, entropy, rho, sample){
    results.names <- c("rds", "plot") 
    results <-vector("list", length(results.names)) |> setNames(results.names)
    # run fit
    x =INCOMMON::classify(x = x,
                       priors = INCOMMON::pcawg_priors,
                       entropy_cutoff = entropy,
                       rho = rho)
    

    results$rds <- x 
    # Plot fit
    results$plot = plot_classification(x, 
                               sample = sample, assembly = T)
    
    return(results)
  }

  classification_plot = reactive({
    if (!is.null(filtered_data_incommon())) {
      do_figure(filtered_data_incommon(), input$entropy, input$rho, input$sample)
    }
  })
  # Render plot
  output$classification_plot <- renderPlot({
    plots <- classification_plot()$plot
    plots  # Return the plot object directly
  })
  
  # Download plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("Classification_Plot-", input$sample,'.pdf', sep = "")
    },
    content = function(file) {
      ggsave(
        file, plot = last_plot(), device = "pdf",
        width = 12, height = 8, units = "in",
        dpi = 300
      )
    }
  )
  # Download incommon fit
  output$downloadRDS <- downloadHandler(
    filename = function() {
      paste("Classification_Object-", input$sample,'.rds', sep = "")
    },
    content = function(file) {
      saveRDS(object = classification_plot()$rds,file = file)
    }
  )
  
}
