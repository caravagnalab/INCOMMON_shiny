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

# Define UI for application that draws a histogram
fluidPage(
  titlePanel("Survival Analysis"),
  
  
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", label = HTML(paste0("Upload an INCOMMON fit object. ", 
                                                "You can also <a href='http://example.com'>download the data</a>.")), accept = c(".rds")),
      checkboxGroupInput("covariates", "Covariates:", choices = choices_covariates),
      selectizeInput("tumorType", "Tumor Type:", choices = NULL),
      selectizeInput("gene", "Gene:", choices = NULL),
      actionButton("plotButton", "Plot"),
      downloadButton("downloadPlot", "Download Plot"),
      width = 3
    ),
    
    mainPanel(
      plotOutput("survivalPlot",width = "70%", height = "800px")
    )
  )
)