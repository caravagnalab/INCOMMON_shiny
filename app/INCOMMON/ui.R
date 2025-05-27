
#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(waiter)
# Load modules
source("./modules/1.classification.R")
source("./modules/2.survival_analysis.R")
source("./modules/3.met_propensity.R")
source("./modules/4.met_tropism.R")
ui <- fluidPage(
  # Title
  div(style = "background-color: #A9A9A7; padding: 10px;",
      titlePanel("Welcome to the INCOMMON ShinyApp", windowTitle = "INCOMMON")),
  div(style = "margin-top: 20px;",
      p(HTML("INCOMMON is a tool for the INference of COpy number and Mutation Multiplicity
      in ONcology. INCOMMON infers the copy number and multiplicity of somatic mutations
      from tumor-only read count data, and can be applied to classify mutations from
      large-size datasets in an efficient and fast way <a href='https://doi.org/10.1101/2024.05.13.24307238'>(Calonaci N., Krisniqi E. et al.)</a>.
      <br/>
      <br/>
      INCOMMON offers a genome
      interpretation framework, in which the full inactivation of tumor suppressor genes
      (TSG) through mutations with LOH, and the enhanced activation of oncogenes through
      mutations with amplification can be detected. These events can then be used to
      perform augmented analysis of survival and metastatic patterns.
      <br/>
      <br/>
      You can use this application to perform the following tasks:
        <ul>
        <li>Classification (Inference of mutation copy number and multiplicity)</li>
        <li>Survival Analysis (Estimate hazard ratios based on INCOMMON interpreted genomes)</li>
        <li>Metastatic propensity (Estimate odds ratios of metastasis based on INCOMMON interpreted genomes)</li>
        <li>Metastatic organotropism (Estimate odds ratios of specific patterns of metastatic spread based on INCOMMON interpreted genomes)</li>
        </ul>
        <br/>
        You can run these analyses on your own data, or download our data and results from our
        <a href='https://zenodo.org/records/10927218'>Zenodo repository</a>.
        <br/>
        For more information about INCOMMON requirements, parameters and output, please visit
             <a href='https://caravagnalab.github.io/INCOMMON'>the INCOMMON website</a>."
             )
        # tags$a(href = "https://caravagnalab.github.io/INCOMMON", "INCOMMON website", target = "_blank")
        )
      ),
  # Subtitle asking the user what they want to do
  div(style = "background-color: #B0BEC5; padding: 10px;",
      h3("Select a task")),

  tags$style(HTML(
    ".classification-button {
      background-image: url('classification.png');
      background-size: contain;
      background-repeat: no-repeat;
      background-position-y: bottom;
      background-position: bottom;
      width: 220px;
      height: 200px;
      border: none;
      cursor: pointer;
    }

    .survival-button {
      background-image: url('survival.png');
      background-size: contain;
      background-repeat: no-repeat;
      background-position-y: bottom;
      background-position: bottom;
      width: 220px;
      height: 200px;
      border: none;
      cursor: pointer;
    }

    .metprop-button {
      background-image: url('met_prop.png');
      background-size: contain;
      background-repeat: no-repeat;
      background-position-y: bottom;
      background-position: bottom;
      width: 220px;
      height: 200px;
      border: none;
      cursor: pointer;
    }

    .mettrop-button {
      background-image: url('met_tropism.png');
      background-size: contain;
      background-repeat: no-repeat;
      background-position-y: bottom;
      background-position: bottom;
      width: 220px;
      height: 200px;
      border: none;
      cursor: pointer;
    }
    "
  )),

  # Action buttons in two columns
  fluidRow(
    column(6, align = "center",
           actionButton(
             inputId = "classification_button",
             label = HTML("<b>Classification</b>"),
             class = "classification-button",
                        )),
    column(6, align = "center",
           actionButton(
             inputId = "survival_button",
             label = HTML("<b>Survival Analysis</b>"),
             class = "survival-button")
           )
    ),
  fluidRow(
    column(6,
           align = "center",
           actionButton(
             "met_propensity_button",
             HTML("<b>Metastatic Propensity</b>"),
             class = "mettrop-button")
           ),
    column(6,
           align = "center",
           actionButton(
             "met_tropism_button",
             HTML("<b>Metastatic Tropism</b>"),
             class = "metprop-button")
           )
    ),

  # Placeholder for rendering the selected shiny app module
  uiOutput("selected_module")
)
