#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Load modules
source("./modules/2.survival_analysis.R")

ui <- fluidPage(
  # Title
  div(style = "background-color: #A9A9A7; padding: 10px;",
      titlePanel("Welcome to INCOMMON ShinyApp", windowTitle = "INCOMMON")),
  # Subtitle asking the user what they want to do
  div(style = "background-color: #B0BEC5; padding: 10px;",
      h3("Use INCOMMON for:")),

  # # Action buttons to choose what to do
  # actionButton("classification_button", "Classification", ),
  # actionButton("survival_button", "Survival Analysis"),
  # actionButton("met_propensity_button", "Metastatic Porpensity"),
  # actionButton("met_tropism_button", "Metastatic Tropism"),

  # CSS to style the action buttons
  tags$head(
    tags$style(HTML("
      .action-btn {
        width: 200px;
        height: 100px;
        margin: 20px;
        font-size: 18px;
        border-radius: 10px;
        position: relative;
      }
      .action-btn .tooltiptext {
        visibility: hidden;
        width: 220px;
        background-color: #555;
        color: #fff;
        text-align: center;
        border-radius: 6px;
        padding: 10px;
        position: absolute;
        z-index: 1;
        bottom: 125%;
        left: 50%;
        margin-left: -110px;
        opacity: 0;
        transition: opacity 0.3s;
      }
      .action-btn:hover .tooltiptext {
        visibility: visible;
        opacity: 1;
      }

      .classification-btn { background-color: #4CAF50; }
      .survival-btn { background-color: #2196F3; }
      .met-propensity-btn { background-color: #f44336; }
      .met-tropism-btn { background-color: #ff9800; }
    "))
  ),

  # Action buttons in two columns
  fluidRow(
    column(6, align = "center",
           actionButton("classification_button", "Classification",
                        class = "action-btn classification-btn",
                        title = "Input: Data for classification.\nOutput: Predicted classes for the data.\n")),
    column(6, align = "center",
           actionButton("survival_button", "Survival Analysis", class = "action-btn survival-btn",
                        title = "Input: Survival data with covariates.\nOutput: Survival curves and hazard estimates.")))
  ,
  fluidRow(
    column(6, align = "center",
           actionButton("met_propensity_button", "Metastatic Propensity", class = "action-btn met-propensity-btn",
                        title = "Input: Clinical and biological data.\nOutput: Predicted propensity of metastasis.")),
    column(6, align = "center",
           actionButton("met_tropism_button", "Metastatic Tropism", class = "action-btn met-tropism-btn",
                        title = "Input: Tumor imaging data.\nOutput: Prediction of metastatic tropism.")))
  ,

  # Placeholder for rendering the selected shiny app module
  uiOutput("selected_module")
)
