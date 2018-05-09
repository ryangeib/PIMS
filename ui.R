#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  fluidRow(
    column(width = 10,
           h1("Peptide Identification utilizing Mass Spectrometry (PIMS) Toolkit")
           , offset = 1
    )
  ),
  # Description
  fluidRow(
    column(width = 10,
      "This tool is designed to test if a constructed peptide is a match to an expected sequence,
      based on comparing the spectra of the peptide to a theoretical spectra based on a proposed
      sequence. Please enter the required inputs below, and then press the execute button."
      , offset = 1
    )
  ),
  hr(),
  fluidRow(
    column(width = 5,
           "Please upload your mass spectrometry output file here.",
           "Make sure your output file is comma-separated (tab-separated files not supported at this time):",
           fileInput("mass_spec", "MS Output File: CSV or TXT",
                     accept = c(
                       "text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv",".txt")
                     ),
           offset = 1
    ),
    column(width = 4,
           textInput("proposed", "Please enter your proposed sequence below:"),
           br(),
           actionButton("execute", "Execute"),
           offset = 1
    )
  ),
  hr(),
  fluidRow(
    column(width = 5,
           plotOutput("Plot"),
           offset = 1
    ),
    column(width = 4,
           textOutput("match", container = span),
           textOutput("details",container = span),
           "A plot of the data is to the left.",
           offset = 1
    )
  ),
  title = "PIMS Toolkit"
))
