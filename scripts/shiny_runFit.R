library(shiny)
library(INCOMMON)
library(dplyr)
library(DT)

# Define UI
ui <- fluidPage(
  tabsetPanel(
    id = "tabs",
    tabPanel("Import data",
             fluidRow(
               column(3,
                      h1("INCOMMON classification"),
                      p("Classify mutations using a Beta-Binomial mixture"),
                      fileInput("upload", "Enter a csv file:"),
                      p("Need purity column"),
                      #################
                      selectInput("sample", "Sample",choices = NULL),
                      p("Select sample of interest"),
                      #################                      
                      sliderInput(inputId = "cutoff",
                                  label = "Entropy cut-off",
                                  value = 0.75,
                                  min = 0,
                                  max = 1,
                                  step = 0.1),
                      p("Entropy cut-off for Tier-1 vs Tier-2 assignment"),
                      #################
                      numericInput(inputId = "rho", label = "RHO", value = 0.01),
                      p("Over-dispersion parameter"),
                      #################
                      textInput(inputId = "tumor_type",label = "Tumor type",value = NA),
                      p("Tumor type"),
                      #################
                      actionButton("submit", "Submit")),
               column(6,tableOutput("data"))
             )
    ),
    tabPanel("Output prediction",
             fluidRow(
               column(6, DTOutput("samplesTable")),
               column(6, plotOutput("plot"))
             ),
             downloadButton("downloadTable", "Download .tsv"),
             # downloadButton("downloadPlot", "Download .pdf")
    )
  )
)


server <- function(input, output, session) {
  
  # Reactive expression to read uploaded CSV file
  my_data <- reactive({
    req(input$upload)
    read.csv(input$upload$datapath)
  })
  # Select a single sample
  observe({
    updateSelectInput(session, "sample", choices = unique(my_data()$sample))
  })
  # Subset the original dataframe according to the
  # select sample
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
      
      input <- init(
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
      classify(x = my_data(), 
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
        # formatStyle(
        #           'VAF',
        #           backgroundColor = styleInterval(c(0.2,0.4,0.6,0.8),
        #                                           c("#FFBFBF","#FF9999","#FF7F7F","#FF3F3F","#FF0000")),
        #           fontWeight = 'bold')
    })
    # Switch to the "Output prediction" tab after submitting
    updateTabsetPanel(session, "tabs", "Output prediction")
    
    # Plot only the gene which
    # is selected from the output of the fit
    plots <- reactive({
      req(out())  # Ensure out is available
      plot_classification(out())
    })
    observeEvent(input$samplesTable_rows_selected, {
      info <- input$samplesTable_rows_selected
      if (length(info) > 0) {
        row <- info[[1]]
        if (row <= length(plots())) {
          output$plot <- renderPlot({
            plot <- plots()[[row]]
            # geneID <- str_extract(plot$labels$title, "\\b\\w+(?= \\(TSG\\))")
            if (!is.null(plot)) {
              plot
            }
          })
        }
      }
    })
  })
  
  # Download button for the output table
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("incomm_prediction", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE,
                append = F,quote = F)
    }
  )

}


# Run the application
shinyApp(ui, server)
