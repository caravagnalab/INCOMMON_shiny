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
