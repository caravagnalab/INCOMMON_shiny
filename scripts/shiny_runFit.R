library(shiny)
library(TAPACLOTH)
library(dplyr)
library(DT)
# Load cancer gene census data
## updated version of the cancer_gene_census table
## this one contains only a unique assigment gene --> gene_role 
cancer_gene_census_new <- readRDS("~/Desktop/dottorato/tapacloth/cancer_gene_census_new.rds")
default_prior <- readRDS(("~/Desktop/dottorato/tapacloth/priors_my_drivers.rds"))
# Define UI
ui <- fluidPage(
  tabsetPanel(
    id = "tabs",
    tabPanel("Import data",
             fileInput("upload", "Enter a csv file:"),
             selectInput("sample", "Sample",choices = NULL),
             actionButton("submit", "Submit"),
             tableOutput("data")
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
  data <- reactive({
    req(input$upload)
    read.csv(input$upload$datapath)
  })
  # Select a single sample
  observe({
    updateSelectInput(session, "sample", choices = unique(data()$sample))
  })
  # Subset the original dataframe according to the
  # select sample
  filtered_data <- reactive({
    req(input$submit, input$sample)
    data() %>% 
      filter(sample == input$sample)
  })
  
  output$data <- renderTable({
    req(input$sample)
    filtered_data()
  })
  # TAPACLOTH part
  # Step1 : initialize the input data
  observeEvent(input$submit, {
    # Reactive expression for tapacloth data
    tapacloth_data <- reactive({
      req(filtered_data())  # Ensure data is available
      
      input <- init(
        mutations = filtered_data(),
        sample = input$sample,
        purity = filtered_data()$purity %>%  unique(),
        gene_roles = cancer_gene_census_new
      )
      input$data <- na.omit(object = input$data)
      input
    })
    
    # Step2 : run the fit
    out <- reactive({
      req(tapacloth_data())  # Ensure tapacloth_data is available
      run_classifier(tapacloth_data(), cutoff = 0.75)
    })
    

    
    
    # Render the output table of the fit
    output$samplesTable <- renderDT({

      datatable(out()$classifier$data %>% dplyr::select(-c(density, p_assign_all,sample,id)) %>%
                  dplyr::mutate(p_assign=round(p_assign,2)) %>%
                  dplyr::mutate(VAF=round(VAF,2)) %>%
                  dplyr::mutate(entropy=round(entropy,2)) %>%
                  dplyr::mutate(mean_entropy=round(mean_entropy,2)),
                selection = "single",filter = "top",options = list(scrollX = TRUE, scrollY = TRUE))%>% 
        formatStyle('VAF',
                    background = styleColorBar(range(0,1), 'lightblue'),
                    backgroundSize = '98% 88%',
                    backgroundRepeat = 'no-repeat',
                    backgroundPosition = 'center')
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
      plot_test(out())
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
