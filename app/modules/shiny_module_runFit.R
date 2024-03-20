library(shiny)
library(INCOMMON)
library(dplyr)
library(DT)
library(ggplot2)
incommonModuleUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    tabsetPanel(
      id = ns("tabs"),
      tabPanel("Import data",
               fluidRow(
                 column(3,
                        h1("INCOMMON classification"),
                        p("Classify mutations using a Beta-Binomial mixture"),
                        fileInput(ns("upload"), "Enter a csv file:"),
                        p("Need purity column"),
                        selectInput(ns("sample"), "Sample", choices = NULL),
                        p("Select sample of interest"),                      
                        sliderInput(ns("cutoff"),
                                    label = "Entropy cut-off",
                                    value = 0.75,
                                    min = 0,
                                    max = 1,
                                    step = 0.1),
                        p("Entropy cut-off for Tier-1 vs Tier-2 assignment"),
                        numericInput(ns("rho"), label = "RHO", value = 0.01),
                        p("Over-dispersion parameter"),
                        textInput(ns("tumor_type"), label = "Tumor type", value = NA),
                        p("Tumor type"),
                        actionButton(ns("submit"), "Submit")),
                 column(6, tableOutput(ns("data")))
               )
      ),
      tabPanel("Output prediction",
               fluidRow(
                 column(6, DTOutput(ns("samplesTable"))),
                 column(6, plotOutput(ns("plot")))
               ),
               downloadButton(ns("downloadTable"), "Download .tsv"),
               downloadButton(ns("downloadPlot"), "Download .pdf")
      )
    )
  )
}

incommonModule <- function(input, output, session) {
  selected_plot <- reactiveVal(NULL)
  
  # Reactive expression to read uploaded CSV file
  my_data <- reactive({
    req(input$upload)
    read.csv(input$upload$datapath)
  })
  
  # Select a single sample
  observe({
    updateSelectInput(session, "sample", choices = unique(my_data()$sample))
  })
  
  # Subset the original dataframe according to the select sample
  filtered_data <- reactive({
    req(input$submit, input$sample)
    my_data() %>% 
      filter(sample == input$sample) %>% 
      select(chr, from, to, ref, alt, gene, NV, DP, VAF, purity)
  })
  
  output$data <- renderTable({
    req(input$sample)
    filtered_data()
  })
  
  # INCOMMON part
  # Step1 : initialize the input data
  observeEvent(input$submit, {
    # Reactive expression for INCOMMON data
    my_data <- reactive({
      req(filtered_data())  # Ensure data is available
      
      input <- INCOMMON::init(
        mutations = filtered_data(),
        sample = input$sample,
        tumor_type = input$tumor_type,
        purity = filtered_data()$purity %>%  unique(),
        gene_roles = INCOMMON::cancer_gene_census
      )
      input$data <- na.omit(object = input$data)
      input
    })
    
    # Step2 : run the fit
    out <- reactive({
      req(my_data())  # Ensure INCOMMON data is available
      INCOMMON::classify(x = my_data(), 
                         priors = INCOMMON::pcawg_priors, 
                         entropy_cutoff = 1,
                         rho = input$rho)
    })
    
    # Render the output table of the fit
    output$samplesTable <- renderDT({
      datatable(
        out() %>% classification() %>% 
          dplyr::mutate(posterior = round(posterior, 2), 
                        entropy = round(entropy, 2),
                        VAF = round(VAF, 2)),
        selection = "single", filter = "top", options = list(scrollX = TRUE, scrollY = TRUE)) %>%
        formatStyle(
          'VAF',
          background = styleColorBar(range(0, 1), 'lightblue'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })
    
    updateTabsetPanel(session, "tabs", "Output prediction")
    
    # Plot only the gene which is selected from the output of the fit
    plots <- reactive({
      req(out())  # Ensure out is available
      plot_classification(out())
    })
    
    observeEvent(input$samplesTable_rows_selected, {
      info <- input$samplesTable_rows_selected
      if (length(info) > 0) {
        row <- info[[1]]
        if (row <= length(plots())) {
          selected_plot(plots()[[row]])
          output$plot <- renderPlot({
            selected_plot()
          })
        }
      }
    })
  })
  
  # Download button for the output table
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("incomm_prediction_", input$sample,".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE,
                append = F, quote = F)
    }
  )
  
  # Download button for the selected plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      plot <- selected_plot()
      gene_str =plot$labels$title %>% strsplit(" ") %>% unlist() 
      geneID <- gene_str[[1]]
      paste0("plot_fit_",input$sample,"_",geneID,".pdf")
    },
    content = function(file) {
      plot <- selected_plot()
      
      if (!is.null(plot)) {
        ggsave(file, plot = plot, device = "pdf", width = 8, height = 6, units = "in", dpi = 300)
      }
    }
  )
}

# Run the application
shinyApp(
  ui = incommonModuleUI("incommon"),
  server = function(input, output, session) {
    callModule(incommonModule, "incommon")
  }
)