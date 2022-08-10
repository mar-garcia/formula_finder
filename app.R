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
      )),
    hr(),
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(4,
                 numericInput(
                   inputId = "ppmt",
                   label = "theoretical m/z:",
                   value = 166.0862,
                   step = 0.0001)),
          column(4,
                 numericInput(
                   inputId = "ppme",
                   label = "experimental m/z value:",
                   value = 166.0872,
                   step = 0.0001))
    )),
    mainPanel(
      fluidRow(column(4, verbatimTextOutput("ppm")))
    ))
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
  
  output$ppm <- renderPrint({
    paste(sprintf("%.2f", round(((input$ppme - input$ppmt)/input$ppmt)*1000000, 
                                2)), "ppm")
  })
}

shinyApp(ui = ui, server = server)