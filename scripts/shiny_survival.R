library(shiny)
library(survival)
library(survminer)
library(dplyr)
library(CNAqc)
library(patchwork)
library(tidyr)
library(tidyverse)
library(patchwork)
library(INCOMMON)

format_p = function(p){
  case_when(
    p < 0.0001 ~ "***",
    p < 0.001 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

ui <- fluidPage(
  titlePanel("Survival Analysis"),

  
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Choose a Data File (RDS):", accept = c(".rds")),
      selectizeInput("tumorType", "Tumor Type:", choices = NULL),
      selectizeInput("gene", "Gene:", choices = NULL),
      actionButton("plotButton", "Plot"),
      downloadButton("downloadPlot", "Download Plot"),
      width = 3
    ),

    mainPanel(
      plotOutput("kmPlot",width = "70%", height = "800px")
      # width = "80%"
    )
  )
)

server <- function(input, output, session) {
  
  data <- reactive({
    req(input$dataFile)
    readRDS(input$dataFile$datapath) 
  })
  
  
  observe({
    req(data())
    
    updateSelectizeInput(session, "tumorType", choices = unique(data()$tumor_type), server = TRUE)
  })
  
  genes <- sort(unique(INCOMMON::cancer_gene_census$gene))
  
  updateSelectizeInput(session, "gene", choices = genes, server = TRUE)
  
  
  observeEvent(input$plotButton, { 
    
    # Input validation 
    req(input$tumorType, input$gene, data()) 
    
    output$kmPlot <- renderPlot({ 
      
      #0 Helper Functions
      
      mutant_samples = function(x, tumor_type, gene) {
        x %>%
          dplyr::filter(tumor_type == !!tumor_type) %>%
          tidyr::separate_rows(genotype, sep = ",\\s*") %>%
          dplyr::filter(grepl(gene, genotype, ignore.case = TRUE)) %>%
          dplyr::group_by(sample) %>%
          dplyr::slice_head(n = 1) %>% 
          dplyr::summarise(group = toString(genotype), across(everything())) %>%
          dplyr::ungroup() %>% 
          dplyr::filter(!grepl('Tier-2', group))
      }
      
      
      
      wt_samples = function(x, tumor_type, gene) {
        x %>%
          dplyr::filter(tumor_type == !!tumor_type) %>%
          dplyr::filter(!grepl(!!gene, genotype)) %>%
          dplyr::mutate(group = paste0(gene, ' WT'))
      }
      
      
      prepare_fit_input = function(x, tumor_type, gene){
        
        rbind(
          mutant_samples(x = x, tumor_type = tumor_type, gene = gene),
          wt_samples(x = x, tumor_type = tumor_type, gene = gene)
        ) %>% 
          mutate(gene = gene) %>% 
          dplyr::left_join(INCOMMON::cancer_gene_census, by = 'gene') %>% 
          dplyr::mutate(gene_role = case_when(
            grepl('TSG', gene_role) ~ 'TSG',
            !grepl('TSG', gene_role) & grepl('oncogene', gene_role) ~ 'oncogene',
            TRUE ~ NA)) %>% 
          dplyr::select(sample, tumor_type, gene, gene_role, dplyr::everything())
        
      }
      
      
      #1 prepare input
      
      
      x = prepare_fit_input(data(), input$tumorType, input$gene)
      
      #2 Define palette
      
      colors = function(gene_role) {
        if(gene_role == 'TSG') 
          out = c('gainsboro', 'forestgreen','goldenrod2')
        else 
          out = c('gainsboro', 'forestgreen','purple3')
        return(out)
      }
      
      
      #3 Apply the KM fit function
      
      km_fit = survival::survfit(
        formula = survival::Surv(OS_MONTHS, OS_STATUS) ~ group,
        data = x %>% dplyr::mutate(group = factor(group, 
                                                  levels = c(
                                                    grep('WT', unique(x$group), value = T),
                                                    grep('Mutant', unique(x$group), value = T) %>% grep('with', ., invert = T, value = T),
                                                    grep('Mutant', unique(x$group), value = T) %>% grep('with', ., value = T)
                                                  ))
        )
      )
      
      
      #4 Adjust KM fit object
      
      names(km_fit$strata) = gsub(km_fit$strata %>% names(), pattern='group=',replacement='')
      
      
      #5 Plot KM fit
      
      km_plot = survminer::ggsurvplot(
        km_fit,
        censor = F,
        conf.int = F,
        data = x %>% as.data.frame(),
        ylab = "Overall Survival",
        xlab = "Time (months)",
        fontsize = 4,
        risk.table = TRUE,
        risk.table.col = "strata",
        table.fontsize = 0.1,
        ggtheme = CNAqc:::my_ggplot_theme(cex = .8),
        tables.theme = CNAqc:::my_ggplot_theme(cex = .8),
        palette = colors(unique(x$gene_role))
      ) 
      
      #6 Adjust KM plot
      
      km_plot$plot$data$tumor_type = unique(x$tumor_type)
      km_plot$data.survplot$tumor_type = unique(x$tumor_type)
      
      km_plot$plot = km_plot$plot + xlab('') + guides(color = 'none') + facet_wrap(~tumor_type)
      km_plot$table = km_plot$table + ylab('') + theme(legend.position = 'none')
      
      #7 Cox fit
      
      colnames(x) %>% sort()
      
      cox_fit = function(x, gene, tumor_type, covariates = c('age', 'sex', 'tmb')){
        
        formula = 'survival::Surv(OS_MONTHS, OS_STATUS) ~ group'
        
        for(c in covariates) {
          what = grep(c, colnames(x), ignore.case = T, value = TRUE)
          for(w in what) {
            if(is.numeric(x[[w]])){
              q =  quantile(x[w], na.rm = T)['50%']
              x[[w]] = ifelse(x[[w]] > q, paste0('>', round(q, 0)), paste0('<=', round(q, 0)))
              x[[w]] = factor(x[[w]])
              x[[w]] = relevel(x[[w]], ref = grep('<=', unique(x[[w]]), value = T))
            }
            formula = paste(formula, w, sep = ' + ')
          }
        }
        
        fit = survival::coxph(
          formula = formula %>% as.formula(),
          data = x %>%
            dplyr::mutate(group = factor(group)) %>% 
            dplyr::mutate(group = relevel(group, ref = grep('WT', unique(x$group), value = T))) %>% 
            as.data.frame()
        )
        
        return(fit)
      }
      
      
      y = cox_fit(x = x, covariates = c('AGE_AT_DEATH', 'tmb', 'MSI', 'sex'))
      
      
      forest_plot = function(x, tumor_types = FALSE){
        
        if(is.null(x)) return(NULL)
        s = summary(x)
        
        x_limits = c(s$conf.int[,'lower .95'] %>% min(),
                     s$conf.int[,'upper .95'] %>% max())
        
        what = s$conf.int %>% as_tibble()
        what$var = rownames(s$conf.int)
        what$p.value = s$coefficients[,ncol(s$coefficients)]
        
        reference_table =
          lapply(names(x$xlevels), function(n) {
            tibble(
              var = paste0(n, x$xlevels[n][[1]][1]),
              value = 1,
              low = 1,
              up = 1,
              p.value = NA
            )
          }) %>% do.call(rbind, .)
        
        toplot = tibble(
          var = what$var,
          value = what$`exp(coef)`,
          low = what$`lower .95`,
          up = what$`upper .95`,
          p.value = what$p.value
        ) %>% 
          rbind(reference_table)
        
        toplot = toplot %>% 
          mutate(var = case_when(
            grepl('group', var) ~ gsub('group', '', var),
            grepl('Sex', var, ignore.case = T) ~ gsub('Sex', 'Sex: ', var, ignore.case = T),
            TRUE ~ var
          )) 
        
        levels = c(x$xlevels$group %>% unique())
        for(c in names(x$xlevels)[-1]){
          levels = c(levels, grep(c, toplot$var, value = T, ignore.case = T) %>% rev())
        }
        
        toplot = toplot %>% 
          dplyr::mutate(var = factor(
            var,
            levels = levels))
        
        toplot$var = factor(toplot$var, levels = levels(toplot$var) %>% rev())
        
        pp = toplot %>% 
          mutate(num_label = toplot$var %>% seq_along()) %>% 
          mutate(stripe = (num_label%%2==0)) %>% 
          ggplot(aes(y = var, x = value))+
          geom_point(aes(color = p.value <= .05))+
          geom_errorbar(aes(xmin = low, xmax = up, color = p.value <= .05), width = .5)+
          geom_rect(aes(
            ymax = num_label + .5,
            ymin = num_label - .5,
            xmin = -Inf,
            xmax = Inf,
            fill = stripe
          ), alpha = .4)+
          geom_vline(xintercept = 1, linetype = 'longdash', alpha = .5)+
          scale_fill_manual(values = c('gainsboro', 'white'))+
          geom_point(aes(color = p.value <= .05))+
          geom_errorbar(aes(xmin = low, xmax = up, color = p.value <= .05), width = .1)+
          geom_text(data = toplot %>% filter(var != 'WT'), 
                    aes(x = 2.5, label = p.value %>% format_p() %>% gsub('ns', '', .)), hjust = +1)+
          # geom_errorbar(aes(xmin = low, xmax = up), width = .1)+
          scale_color_manual(values = c('indianred3','black') %>% rev())+
          scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = x_limits)+
          CNAqc:::my_ggplot_theme(cex = .8)+
          ylab('')+
          xlab('Hazard Ratio')+
          guides(fill = 'none')
        
        
        pp
      }
      
      library(patchwork)
      
      km_plot$plot/
        km_plot$table/
        forest_plot(y) +
        plot_layout(heights = c(3.5,1,2))
      
    })
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("KM_Plot-", input$tumorType, '-', input$gene, '.png', sep = "")
    },
    content = function(file) {
      
      ggsave(file, plot = last_plot(), device = "png", width = 8, height = 6, units = "in", dpi = 300)
    }
  )
  
}

shinyApp(ui = ui, server = server)
