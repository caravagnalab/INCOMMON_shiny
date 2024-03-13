# --- LIBRARIES ---

library(shiny)
library(survival)
library(survminer)
library(dplyr)
library(CNAqc)
library(patchwork)

# --- FILES --- PC IFO: C:\Users\3704\Dropbox\2023.TAPACLOTH\results\survival_data.rds

# --- SHINY APP ---

ui <- fluidPage(
  titlePanel("Kaplan-Meier Survival Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("dataPath", "Data File Path (RDS):", value = ""),
      textInput("tumorType", "Tumor Type:", value = ""),
      textInput("gene", "Gene:", value = ""),
      actionButton("plotButton", "Plot")
    ),
    
    mainPanel(
      plotOutput("kmPlot")
    )
  )
)

server <- function(input, output) {
  
  data <- reactive({
    req(input$dataPath)
    readRDS(input$dataPath) 
  })
  
  observeEvent(input$plotButton, { 
    
    # Input validation 
    req(input$tumorType, input$gene, data()) 
    
    output$kmPlot <- renderPlot({ 
      
      #0 Helper Functions
      
      mutant_samples = function(x, tumor_type, gene) {
        x %>%
          dplyr::filter(tumor_type == !!tumor_type) %>%
          dplyr::filter(gene == !!gene) %>%
          dplyr::mutate(genotype = dplyr::case_when(
            class == 'Other' ~ paste0('Tier 2'),
            grepl('\\+CNA', class) & gene_role == 'TSG' ~ paste0('Mutant ', gene, ' with LOH'),
            grepl('\\+CNA', class) & gene_role != 'TSG' ~ paste0('Mutant ', gene, ' with AMP'),
            grepl('\\-CNA', class) & gene_role == 'TSG' ~ paste0('Mutant ', gene),
            grepl('\\-CNA', class) & gene_role != 'TSG' ~ paste0('Mutant ', gene)
          ))
      }
      
      wt_samples = function(x, tumor_type, gene) {
        x %>%
          dplyr::filter(tumor_type == !!tumor_type) %>%
          dplyr::filter(!grepl(!!gene, genotype)) %>%
          dplyr::mutate(genotype = 'WT')
      }
      
      prepare_fit_input = function(x, tumor_type, gene) {
        rbind(
          mutant_samples(x = x, tumor_type = tumor_type, gene = gene) %>% 
            dplyr::select(sample, tumor_type, genotype, dplyr::starts_with('OS'), SEX, AGE_AT_DEATH, MSI_SCORE, TMB_NONSYNONYMOUS) %>% 
            dplyr::rename(Age = AGE_AT_DEATH, Sex = SEX, MSI = MSI_SCORE, TMB = TMB_NONSYNONYMOUS) %>% 
            dplyr::rename(group = genotype),
          wt_samples(x = x, tumor_type = tumor_type, gene = gene) %>% 
            dplyr::select(sample, tumor_type, genotype, dplyr::starts_with('OS'), SEX, AGE_AT_DEATH, MSI_SCORE, TMB_NONSYNONYMOUS) %>% 
            dplyr::rename(Age = AGE_AT_DEATH, Sex = SEX, MSI = MSI_SCORE, TMB = TMB_NONSYNONYMOUS) %>% 
            dplyr::rename(group = genotype)
        )
      }
      
      
      #1 prepare input
      
      x <- prepare_fit_input(data(), input$tumorType, input$gene) %>% 
        dplyr::filter(group != 'Tier 2')
      
      #2 Define palette
      
      gene_role = data() %>% 
        dplyr::filter(tumor_type == input$tumorType, gene == input$gene) %>%
        pull(gene_role) %>% 
        unique()
      
      colors = function(gene_role) {
        if(gene_role == 'TSG') 
          out = c('forestgreen','goldenrod2','gainsboro')
        else 
          out = c('forestgreen','purple3','gainsboro')
        return(out)
      }
      
      #3 Apply the fit function
      
      km_fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ group, data = x)
      
      #4 Modify aspects of the fit object
      
      names(km_fit$strata) <- gsub('group=', '', names(km_fit$strata))
      
      #5 Plot the fit
      
      km_plot <- survminer::ggsurvplot(
        km_fit,
        censor = F,
        conf.int = F,
        data = x,
        ylab = "Overall Survival",
        xlab = "Time (months)",
        fontsize = 4,
        risk.table = TRUE,
        risk.table.col = "strata",
        table.fontsize = 0.1,
        ggtheme = CNAqc:::my_ggplot_theme(cex = .8),
        tables.theme = CNAqc:::my_ggplot_theme(cex = .8),
        palette = colors(gene_role)
      ) 
      
      #6 Modify the plot object
      
      km_plot$plot$data$tumor_type = unique(x$tumor_type)
      km_plot$data.survplot$tumor_type = unique(x$tumor_type)
      
      km_plot$plot = km_plot$plot + xlab('') + guides(color = 'none') + facet_wrap(~tumor_type)
      km_plot$table = km_plot$table + ylab('') + theme(legend.position = 'none')
      
      km_plot$plot/
        km_plot$table+
        plot_layout(heights = c(3.5,1))
      
    })
  })
}

shinyApp(ui = ui, server = server)
