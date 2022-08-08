options(repos = BiocManager::repositories())


library(shiny)
library(MetaboCoreUtils)

ui <- navbarPage(
  "",
  tabPanel(
    "",
    sidebarLayout(
      sidebarPanel(
        textInput("formula", label = "Formula:", 
                  value = "C28H37N5O7")
      ),
      mainPanel(
        fluidRow(column(8, verbatimTextOutput("neutral"))),
        fluidRow(column(8, verbatimTextOutput("pos"))),
        fluidRow(column(8, verbatimTextOutput("neg")))
      )
    )
    ))
server <- function(input, output){
  
  output$neutral <- renderPrint({ 
    paste("Monoisotopic mass:", sprintf("%.5f", calculateMass(input$formula)))
  })
  
  output$pos <- renderPrint({ 
    mzneutral <- calculateMass(input$formula)
    unlist(mass2mz(mzneutral, 
                   adduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[2M+H]+"))) 
  })
  
  output$neg <- renderPrint({ 
    mzneutral <- calculateMass(input$formula)
    unlist(mass2mz(mzneutral, 
                   adduct = c("[M-H]-", "[M+Cl]-", "[M+CHO2]-", "[2M-H]-")))
  })
}

shinyApp(ui = ui, server = server)