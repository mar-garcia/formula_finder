options(repos = BiocManager::repositories())


library(shiny)
library(MetaboCoreUtils)
library(Rdisop)

ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}

adds <- rbind(
  adducts("positive"),
  "[M+CH2O2Na]-" = c("[M+CH2O2Na]-", 1, calculateMass("CH2O2Na"), "H", "H", 1, TRUE),
  adducts("negative"),
  "[M+CO2Na]-" = c("[M+CO2Na]-", 1, calculateMass("CO2Na"), "H", "H", -1, FALSE),
  "[2M+Na-2H]-" = c("[2M+Na-2H]-", 2, calculateMass("Na") - calculateMass("H")*2, "H", "H", -1, FALSE),
  "[M-H-H2O]-" = c("[M-H-H2O]-", 1, - calculateMass("HH2O"), "H", "H", -1, FALSE)
)
adds$mass_multi <- as.numeric(adds$mass_multi)
adds$mass_add <- as.numeric(adds$mass_add)
adds$charge <- as.numeric(adds$charge)
adds$positive <- as.logical(adds$positive)


ui <- navbarPage(
  "",
  tabPanel(
    "Main",
    sidebarLayout(
      sidebarPanel(
        textInput("formula", label = "Formula:", 
                  value = "C28H37N5O7")
      ),
      mainPanel(
        fluidRow(verbatimTextOutput("neutral")),
        fluidRow(verbatimTextOutput("pos")),
        fluidRow(verbatimTextOutput("neg"))
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
                   label = "experimental m/z:",
                   value = 166.0872,
                   step = 0.0001))
        )),
      mainPanel(
        fluidRow(verbatimTextOutput("ppm"))
      )),
    hr(),
    sidebarLayout(
      sidebarPanel(
        fluidRow(column(4,
                        numericInput(
                          inputId = "mzt",
                          label = "theoretical m/z:",
                          value = 166.0862,
                          step = 0.0001
                        )),
                 column(4,
                        numericInput(
                          inputId = "ppmdev",
                          label = "ppm deviations:",
                          value = 5,
                          step = 0.1
                        )))
      ),
      mainPanel(
        fluidRow(verbatimTextOutput("mzrange")),
        fluidRow(verbatimTextOutput("mzdif"))
      )
    )
  ), # close tabPanel Main
  # Formula finder ----
  tabPanel(
    "Formula Finder",
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(3,
                 numericInput(inputId = "mz1",
                              label = "m/z value",
                              value = 205.0971)),
          column(3,
                 numericInput(inputId = "i1",
                              label = "intensity",
                              value = 4854167))
        ),
        fluidRow(
          column(3,
                 numericInput(inputId = "mz2",
                              label = "",
                              value = 206.1005)),
          column(3,
                 numericInput(inputId = "i2",
                              label = "",
                              value = 631041))
        ),
        fluidRow(
          column(3,
                 numericInput(inputId = "mz3",
                              label = "",
                              value = 207.1038)),
          column(3,
                 numericInput(inputId = "i3",
                              label = "",
                              value = 48542))
        ),
        
        hr(),
        fluidRow(
          column(4,
                 radioButtons("p", label = "Polarity:",
                              choices = list("POS" = "POS", "NEG" = "NEG"), 
                              selected = "POS")),
          column(5, radioButtons("adduct", label = "Adduct:",
                                 choices = list("[M+H]+" = "[M+H]+", 
                                                "[M-H]-" = "[M-H]-"), 
                                 selected = "[M+H]+"))
        ),
        hr(),
        fluidRow(
          sliderInput("error_C", "% Error in C", min = 0, max = 100, value = 10)
        )
      ),
      mainPanel(
        h3("Formulas suggested:"),
        DT::dataTableOutput("table"),
        hr(),
        fluidRow(
          column(6,
                 h3("Isotopic  pattern:"),
                 textInput("form", label = "", 
                           value = "C11H12N2O2"),
                 plotOutput(outputId = "isopatt")
          ))
      ))), # close tabPanel Formula Finer
  # Instructions ----
  tabPanel("Instructions",
           p("This shiny app contain 2 main tools: 'Formula Finder' and 'Adduct calculation'."),
           p("If you have any question and/or any suggestion to improve this shiny app, please write me at mar.garcia@fmach.it. Your feedback will be very appreciated! :)"),
           br(),
           h3("Formula Finder"),
           ("This tool has been developed to help during the process of identifying features derived from HRMS experiments. The idea is to use this tool once we have a more or less clear hypothesis about what is the assignment for a certain "),
           em("m/z"), ("value. So, we will write its "),
           em("m/z"), ("value at the top-left of the screen and also indicate which adduct we think it refers to. Considering the deviation in ppm indicated in the corresponding box, the tool will calculate the possible formulas to which the "), 
           em("m/z"), (" value we are interrogating may correspond (using the R 'Rdisop' package). These formulas are included in the table on the right. This table can be filtered (by checking the corresponding boxes) by different analytical chemistry criteria. Following there is the description of these rules."),
           br(),
           br(),
           strong("Experimental rules"),
           p("Here we use the isotopic pattern to evaluate how our 'experimental' isotopic pattern (indicated on the left part of the page) fits with the 'theoretical' isotopic pattern of a given formula (indicated just above the plot). The main adducts of this 'theoretical' formula are printed on the right side of the plot. The sliders allow to fix an 'error' window when filtering for one or both isotopes."),
           br(),
           strong("Heuristic rules"),
           br(),
           ("The "), 
           code("Hydrogen Rule"),
           (" determines the maximum number of hydrogens that a formula can contain according to the following calculation:"),
           em("Max(H) = 2C + N + 2"),
           (". Note that this rule is never applied for formulas which contain phosphorus."),
           br(),
           ("The "),
           code("Refined Hydrogen Rule"),
           (" is based on the fact that the sum of the nominal mass plus the number of H-atoms is divisible by four."),
           br(),
           ("The "),
           code("Nitrogen Rule"),
           (" states that a compound with an odd molecular weight have an odd number of nitrogen's and compounds with an even molecular weight will have either no nitrogen or an even number of nitrogen atoms."),
           br(),
           ("The number of "),
           code("Double Bond Equivalent (DBE) Rule"),
           ( "is calculated according to the following expression: "),
           em("DBE = (C+Si) - 1/2*(H+Cl+F+I) + 1/2(N+P) + 1"),
           (". This filter is based on the fact that its answer is a positive integer value (i.e., it can never be negative or a fraction)."),
           br(),
           ("There is also the possibility to restrict formulas according to their"),
           code("element ratios"),
           (". Default ratios are set in a way that all compounds described at the HMDB fill the criteria."),
           br(),
           ("Formulas are "),
           code("ranked"),
           (" according the absolute mean deviation in the relative intensity of both isotopes."),
           br(),
           br(),
           h3("Adduct calculation"),
           ("The idea of the tool provided in the upper part of this page is that you can calculate the "), em("m/z"), 
           (" values for different adducts given an specific "), em("m/z"), ("value and its assignation."), br(),
           ("In the lower part of the page you can get the "), em("m/z"), 
           ("values for different adducts given an specific formula."),
           br(),
           br(),
           h3("References"),
           ("CHROMacademy. "),
           strong("LC-MS Interpretation:"),
           a("https://www.chromacademy.com/channels/lc-ms/technique/lc-ms-interpretation/"),
           br(),
           ("Kind T, Fiehn O. "),
           strong("Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry. "), 
           em("BMC Bioinformatics,"), (" 2007;8:105."),
           br(),
           ("Claesen J, Valkenborg D, Burzykowski T. "),
           strong("'Refined Hydrogen Rule' and a 'Refined Hydrogen and Halogen Rule' for Organic Molecules. "),
           em("J Am Soc Mass Spectrom,"), (" 2020;31(1):132-6."),
           br(),
           hr(),
           br(),
           ("This shiny app has been inspired by the Mass Decomposition tool developed by Jan Stanstrup: "),
           a("http://predret.org/tools/mass-decomposition/")
  )# close tabPanel Instructions
)
server <- function(input, output){
  
  output$neutral <- renderPrint({ 
    paste("Monoisotopic mass:", sprintf("%.5f", calculateMass(input$formula)))
  })
  
  output$pos <- renderPrint({ 
    mzneutral <- calculateMass(input$formula)
    unlist(mass2mz(mzneutral, 
                   adds[c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+CH2O2Na]-",
                              "[2M+H]+", "[2M+Na]+", 
                              "[M+H-H2O]+"),])) 
  })
  
  output$neg <- renderPrint({ 
    mzneutral <- calculateMass(input$formula)
    unlist(mass2mz(mzneutral, 
                   adds[c("[M-H]-", "[M+Cl]-", "[M+CHO2]-", "[M+CO2Na]-", 
                          "[2M-H]-", "[2M+Na-2H]-", "[M-H-H2O]-"),]))
  })
  
  output$ppm <- renderPrint({
    paste(sprintf("%.2f", ((input$ppme - input$ppmt)/input$ppmt)*1000000), 
          "ppm")
  })
  
  output$mzrange <- renderPrint({
    paste("m/z range:", 
          sprintf("%.5f", input$mzt - (input$ppmdev*input$mzt)/1e6), "-", 
          sprintf("%.5f", input$mzt + (input$ppmdev*input$mzt)/1e6))
  })
  
  output$mzdif <- renderPrint({
    paste("m/z difference:",
          sprintf("%.6f", (input$mzt + (input$ppmdev*input$mzt)/1e6) - 
                    (input$mzt - (input$ppmdev*input$mzt)/1e6)), "Da")
  })
  
  
  # Formula Finder -----
  
  fml <- reactive({
    ip <- data.frame(rbind(c(input$mz1,	input$i1),	
                           c(input$mz2,	input$i2),
                           c(input$mz3,	input$i3)
    ))
    colnames(ip) <- c("mz", "i")
    fn <- 100/ip$i[1] # normalization factor
    ip$i <- ip$i*fn # normalize intensities
    ip$d_mz <- ip$mz - ip$mz[1] # calculate deltas in m/z values
    # deduce to which element correspond each isotope:
    ip$E <- NA
    ip$E[which(abs(1.003355 - ip$d_mz) < 0.001)] <- "C"
    ip$E[which(abs(2.0043 - ip$d_mz) < 0.001)] <- "O"
    
    # get the isotope corresponding to 13C:
    idx <- which(ip$E == "C")
    nc <- ip$i[idx] / 1.1 # calculate the number of C
    # number of C +- 10% error:
    nc <- seq(round(nc - nc*input$error_C/100), round(nc + nc*input$error_C/100)) 
    # get the nominal mass of the main ion:
    ms <- as.numeric(mz2mass(input$mz1, input$adduct))
    
    # get all potential combinations of elements C-O-N-S:
    fml <- data.frame(
      C = rep(nc, each = 21*6*3*3),
      O = rep(seq(0, 20), length(nc), each = 6*3*3),
      N = rep(rep(seq(0, 5), each = 3*3), length(nc)*21),
      S = rep(rep(seq(0, 2), each = 3), length(nc)*21*6),
      P = rep(seq(0, 2), length(nc)*21*6*3)
    )# calculate the number of H considering the nominal mass:
    fml$H <- floor(ms) - (fml$C*12 + fml$O*16 + fml$N*14 + fml$S*32 + fml$P*31)
    # check the different rules:
    fml$H_rule <- fml$H <= (2*fml$C + fml$N + 2) # 2C + N + 2
    fml$N_rule <- fml$N %% 2 == floor(ms) %% 2
    fml$DBE <- fml$C - fml$H/2 + (fml$N + fml$P)/2 + 1 # (C+Si) - ?(H+cl+Fl+I) + ?(N+P) + 1
    
    # keep the formulas fullfiling the rules:
    fml <- fml[fml$H > 0 & fml$H_rule == T & fml$N_rule == T & fml$DBE >= 0 & (fml$DBE %% 1) == 0, ]
    # get the complete formula:
    fml$formula <- paste0("C", fml$C, "H", fml$H, "N", fml$N, "O", fml$O, 
                          "S", fml$S, "P", fml$P)
    fml$formula <- gsub("O0", "", fml$formula)
    idx <- which(fml$O == 1)
    fml$formula[idx] <- gsub("O1", "O", fml$formula[idx])
    fml$formula <- gsub("N0", "", fml$formula)
    fml$formula <- gsub("N1O", "NO", fml$formula)
    fml$formula <- gsub("S0", "", fml$formula)
    fml$formula <- gsub("S1", "S", fml$formula)
    fml$formula <- gsub("P0", "", fml$formula)
    fml$formula <- gsub("P1", "P", fml$formula)
    # calculate the theoretical m/z value of the main ion:
    fml$mz <- as.numeric(mass2mz(calculateMass(fml$formula), input$adduct))
    # calculate the deviations in ppm:
    fml$ppm <- round((abs(fml$mz - ip$mz[1]) / fml$mz)*1e6, 4)
    fml <- fml[order(fml$ppm), ]
    fml
  })
  
  output$table <- DT::renderDataTable(DT::datatable({
    fml()
  },
  options = list(pageLength = 5, dom = "tip" ),
  rownames= FALSE))
  
  isopt <- reactive({
    if(input$p == "POS"){
      t_fml <- addElements(input$form, "H")
    }else if(input$p == "NEG"){
      t_fml <- subtractElements(input$form, "H")
    }
    
    dt <- data.frame(t(getIsotope(getMolecule(t_fml))[,1:3]))
    dt$X3 <- (dt$X2/max(dt$X2))*100
    dt
  })
  
  output$isopatt <- renderPlot({
    plot(c(input$mz1, input$mz2, input$mz3),
         c(100, (input$i2/input$i1)*100, (input$i3/input$i1)*100),
         type = "h", xlab = "m/z value", ylab = "Relative intensity", 
         col = 2, lwd = 2, xaxt = 'n', 
         xlim = c(min(isopt()[1,"X1"], input$mz1) - 0.3, 
                  max(isopt()[3,"X1"], input$mz3) + 0.3)
    )
    axis(1, at = isopt()[,"X1"], 
         labels = c(
           paste0(round(isopt()[1,"X1"],4), " (", 
                  round((abs(isopt()[1,"X1"] - input$mz1)/isopt()[1,"X1"])*1e6, 1), " ppm)"),
           paste0(round(isopt()[2,"X1"],4), " (", 
                  round((abs(isopt()[2,"X1"] - input$mz2)/isopt()[2,"X1"])*1e6, 1), " ppm)"),
           paste0(round(isopt()[3,"X1"],4), " (", 
                  round((abs(isopt()[3,"X1"] - input$mz3)/isopt()[3,"X1"])*1e6, 1), " ppm)")
         ))
    rect(xleft = isopt()[,"X1"] - (ppm(isopt()[,"X1"], 100)), 
         xright =isopt()[,"X1"] + (ppm(isopt()[,"X1"], 100)), 
         ybottom = 0, ytop = isopt()[,"X3"], border = 3)
    legend("topright", legend = c("Theoretical", "Experimental"),
           col = c("#B2DF8A", "#E31A1C"), lty = 1, lwd = 3)
    
  })
  
}

shinyApp(ui = ui, server = server)