#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)



# Define server logic required to execute program
shinyServer(function(input, output) {
  
  # Adjust Maximium upload size to 30MB.
  # Code derived from stackoverflow.com. Full source included in final documentation.
  options(shiny.maxRequestSize=30*1024^2)
  # Sources for other members' code
  source("dataReader.R")
  source("aminoAcidToThreeHighest.R")
  source("compareThreeHighest.R")
  source("countAminoAcids.R")
  source("matchAminoAcidsAlgorithm.R")
  source("mSDistributionBF.R")
  source("threeHighestVal.R")
  
  matching <- reactive(
    {
      matchTest <- TRUE
      
      msMatrix <- matrix(c(0,0,0,0), ncol=2, nrow=2)
      
      # First, we read the uploaded file (if one exists)
      upload <- input$mass_spec
      if (is.null(upload$datapath)) {
        msMatrix <- matrix(c(0,0,0,0), ncol=2, nrow=2)
      } else {
        # Data Reading (Zach)
        msMatrix <- dataReading(upload$datapath)
      }
      
      # Peak Extraction (Ryan)
      threeHighestDataPeaks <- threeHighestVal(msMatrix)
      
      # Proposed Sequence reading (Ryan)
      AAVec <- countAminoAcids(input$proposed)
      
      # Comparison Function (Ryan)
      matchingSpectrum <- matchAminoAcidsAlgorithm(AAVec,threeHighestDataPeaks)
      
      if (matchingSpectrum[7] == "No Match Found") {
        matchTest <- FALSE
      }
      
      # Plot Output
      output$Plot <- renderPlot({
        # Call spectplot (Alyssa)
        # spectplot(msMatrix)
        cleanedSpectrum <- c(matchingSpectrum[1],matchingSpectrum[2],matchingSpectrum[3],matchingSpectrum[4],matchingSpectrum[5],matchingSpectrum[6])
        plot(msMatrix,type = "h", main = "Data Plot")
      })
      
      
      # Text Outputs
      if(matchTest) {
        output$match <- renderText("Congratulations! The proposed sequence and the produced sequence were a match!")
      } else {
        output$match <- renderText("The proposed sequence and the produced sequence were not a match.")
      }
      
      if (matchingSpectrum[7] == "No Match Found" | matchingSpectrum[7] == "Sequence Matches") {
        output$details <- renderText("")
      } else {
        detail <- matchingSpectrum[7]
        output$details <- renderText(detail)
      }
      
    }
  )
  
  observeEvent(input$execute, matching())
  
})
