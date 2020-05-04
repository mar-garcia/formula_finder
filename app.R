options(repos = BiocManager::repositories())

library(Rdisop)
library(MetaboCoreUtils)
library(shiny)
library(CompoundDb)

# FUNCTIONS -----------
getmzneut <- function(
  mz = numeric(0),
  adduct = c("[M+H]+", "[M+Na]+", "[M+K]+", 
             "[M-H]-", "[M-H+HCOOH]-", "[M+Cl]-", "[M-H+HCOONa]-",
             "[2M-H]-", "[2M+H]+",
             "[M+H-H2O]+")){
  addtb <- rbind(
    c("[M+H]+", 1.007276),
    c("[M+Na]+", 22.98980),
    c("[M+K]+", 38.96371),
    c("[M-H]-", -1.007276),
    c("[M-H+HCOOH]-", 44.99820),
    c("[M+Cl]-", 34.96885),
    c("[M-H+HCOONa]-", 66.98017),
    c("[2M+H]+", 1.007276),
    c("[2M-H]-", -1.007276),
    c("[M+H-H2O]+", -17.00328)
  )
  if(grepl("2M", adduct)){
    (mz - as.numeric(addtb[addtb[,1] == adduct,2]))/2
  } else {
    mz - as.numeric(addtb[addtb[,1] == adduct,2])
  }
}


getform <- function(mzval = numeric(0),
                    elements = "CHNOPS",
                    ppm = 10){
  Rdisop::decomposeMass(mzval, 
                        ppm = ppm,
                        elements = initializeElements(
                          unlist(strsplit(elements, 
                                          "(?<=.)(?=[[:upper:]])", 
                                          perl=TRUE)))
  )$formula
}

update_form <- function(
  formulas = character(0),
  adduct = c("[M+H]+", "[M+Na]+", "[M+K]+", 
             "[M-H]-", "[M-H+HCOOH]-", "[M+Cl]-", "[M-H+HCOONa]-",
             "[2M-H]-", "[2M+H]+",
             "[M+H-H2O]+"),
  action = c("add", "remove")){
  addtb <- data.frame(rbind(
    c("[M+H]+", "H"),
    c("[M+Na]+", "Na"),
    c("[M+K]+", "K"),
    c("[M-H]-", "H"),
    c("[M-H+HCOOH]-", "CH2O2"),
    c("[M+Cl]-", "Cl"),
    c("[M-H+HCOONa]-", "CHO2Na"),
    c("[2M+H]+", "H"),
    c("[2M-H]-", "H"),
    c("[M+H-H2O]+", "H2O")
  ))
  colnames(addtb) <- c("adduct", "formulas")
  addtb <- cbind(
    addtb, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(addtb$formula), 
                   MetaboCoreUtils::countElements)))
  addtb$H[grep("M-H", addtb$adduct)] <- 
    addtb$H[grep("M-H", addtb$adduct)] - 1
  addtb$H[addtb$adduct == "[M-H]-"] <- -1
  addtb$H[addtb$adduct == "[2M-H]-"] <- -1
  addtb$H[addtb$adduct == "[M+H-H2O]+"] <- -1
  addtb$O[addtb$adduct == "[M+H-H2O]+"] <- -1
  
  tmp.frm1 <- data.frame(formulas)
  tmp.frm1$formulas <- as.character(tmp.frm1$formulas)
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  
  addtb[setdiff(names(tmp.frm1), names(addtb))] <- NA
  tmp.frm1[setdiff(names(addtb), names(tmp.frm1))] <- NA
  tmp.frm1 <- rbind(addtb[addtb$adduct == adduct,], tmp.frm1)
  tmp.frm1 <- tmp.frm1[,-1]
  tmp.frm1[is.na(tmp.frm1)] <- 0
  idx <- ncol(tmp.frm1)
  tmp.frm1$formok <- NA
  if(action == "add"){
    if(adduct %in% c("[2M+H]+", "[2M-H]-")){
      for(i in 2:nrow(tmp.frm1)){
        tmp.frm1[i,2:idx] <- tmp.frm1[i,2:idx]*2
      }
    }
    for(i in 2:nrow(tmp.frm1)){
      tmp.frm1[i, 2:idx] <- tmp.frm1[i, 2:idx] + tmp.frm1[1, 2:idx]
      tmp.frm1$formok[i] <- pasteElements(tmp.frm1[i, 2:idx])
    }
  } else if(action == "remove"){
    for(i in 2:nrow(tmp.frm1)){
      tmp.frm1[i, 2:idx] <- tmp.frm1[i, 2:idx] - tmp.frm1[1, 2:idx]
    }
    if(adduct %in% c("[2M+H]+", "[2M-H]-")){
      for(i in 2:nrow(tmp.frm1)){
        tmp.frm1[i,2:idx] <- tmp.frm1[i,2:idx]/2
      }
    }
    for(i in 2:nrow(tmp.frm1)){
      tmp.frm1$formok[i] <- pasteElements(tmp.frm1[i, 2:idx])
    }
  }
  tmp.frm1$formok[-1]
}


H_rule <- function(formulas = character(0)){
  
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "CHNOPS")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  tmp.frm1$rule <- tmp.frm1$H <= tmp.frm1$C*2 + tmp.frm1$N + 2
  tmp.frm1$rule[tmp.frm1$P > 0] <- TRUE
  tmp.frm1$rule
}


H_rule2 <- function(mzval = numeric(0),
                    formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "CHNOPS")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  tmp.frm1$rule <- (trunc(mzval) + tmp.frm1$H) %% 4 == 0
  tmp.frm1$rule[tmp.frm1$P > 0] <- TRUE
  tmp.frm1$rule
}


N_rule <- function(mzval = numeric(0),
                   formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "N")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  tmp.frm1$N %% 2 == round(mzval) %% 2
}


RPU_rule <- function(formulas = character(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1 <- rbind(tmp.frm1, "CHNOPSSiClFlI")
  tmp.frm1 <- cbind(
    tmp.frm1, 
    do.call(dplyr::bind_rows, 
            lapply(as.character(tmp.frm1$formula), 
                   MetaboCoreUtils::countElements)))
  tmp.frm1[is.na(tmp.frm1)] <- 0
  tmp.frm1 <- tmp.frm1[-nrow(tmp.frm1),]
  
  # RPU <- (C+Si) - (H+Cl+Fl+I)/2 + (N+P)/2 + 1
  (
    (tmp.frm1$C + tmp.frm1$Si) - 
      ((tmp.frm1$H + tmp.frm1$Cl + tmp.frm1$Fl + tmp.frm1$I)/2) + 
      ((tmp.frm1$N + tmp.frm1$P)/2) + 
      1
  )  >= 0
}

isotope_rule <- function(formulas = character(0),
                         intensity = numeric(0),
                         range = 0.4,
                         isotope = 1){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1$A1 <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmpx <- Rdisop::getMolecule(tmp.frm1$formula[i])$isotopes[[1]]
    tmp.frm1$A1[i] <- (tmpx[2,isotope+1]/tmpx[2,1])*100
  }
  A1_range <- intensity + (intensity * range * c(-1, 1))
  tmp.frm1$A1 >= A1_range[1] & tmp.frm1$A1 <= A1_range[2]
}

rank_form <- function(formulas = character(0),
                      i1 = numeric(0),
                      i2 = numeric(0)){
  tmp.frm1 <- data.frame(formula = formulas)
  tmp.frm1$formula <- as.character(tmp.frm1$formula)
  tmp.frm1$A1 <- NA
  tmp.frm1$A2 <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmpx <- Rdisop::getMolecule(tmp.frm1$formula[i])$isotopes[[1]]
    tmp.frm1$A1[i] <- (tmpx[2, 2] / tmpx[2, 1])*100
    tmp.frm1$A2[i] <- (tmpx[2, 3] / tmpx[2, 1])*100
  }
  tmp.frm1$A1d <- abs(tmp.frm1$A1 - i1)
  tmp.frm1$A2d <- abs(tmp.frm1$A2 - i2)
  tmp.frm1$Ad <- NA
  for(i in 1:nrow(tmp.frm1)){
    tmp.frm1$Ad[i] <- mean(c(tmp.frm1$A1d[i], tmp.frm1$A2d[i]))
  }
  tmp.frm1 <- tmp.frm1[order(tmp.frm1$Ad), ]
  seq(nrow(tmp.frm1))
}



# UI ------------------------------
ui <- navbarPage(
  "",
  
  # Formula finder -----------
  tabPanel("Formula Finder",
           
           # Options ------------------------------------
           sidebarLayout(
             
             sidebarPanel(
               
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz",
                                     label = "m/z value:",
                                     value = 205.0971,
                                     step = 0.0001)),
                 
                 column(4, 
                        numericInput(inputId = "ppm",
                                     label = "ppm:",
                                     value = 10,
                                     step = 1))
               ),
               
               radioButtons(
                 inputId = "adduct", 
                 label = "Adduct:",
                 choices = list("[M+H]+" = 1, "[M+K]+" = 2, 
                                "[M-H]-" = 3, "[M-H+HCOOH]-" = 4), 
                 selected = 1),
               
               hr(),
               
               h3("Heuristic rules:"),
               
               checkboxInput("H_rule", 
                             label = "H rule", 
                             value = TRUE),
               
               checkboxInput("H_rule2", 
                             label = "H rule (II)", 
                             value = TRUE),
               
               checkboxInput("N_rule", 
                             label = "N rule", 
                             value = TRUE),
               
               checkboxInput("RPU_rule", 
                             label = "RPU rule", 
                             value = TRUE),
               
               hr(),
               
               h3("Experimental rules:"),
               fluidRow(
                 column(3,
                        checkboxInput("iso_1", 
                                      label = "Isotope 1", 
                                      value = TRUE)),
                 column(9,
                        sliderInput("error1", label = "", min = 0, 
                                    max = 10, value = 0.4, step = 0.1))
               ),
               fluidRow(column(3,
                               checkboxInput("iso_2", 
                                             label = "Isotope 2", 
                                             value = TRUE)),
                        column(9,
                               sliderInput(
                                 "error2", label = "", min = 0, 
                                 max = 10, value = 2, step = 0.1))
               ),
               
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz1",
                                     label = "m/z value",
                                     value = 205.0971)),
                 column(6,
                        numericInput(inputId = "i1",
                                     label = "intensity",
                                     value = 100))
               ),
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz2",
                                     label = "",
                                     value = 206.1005)),
                 column(6,
                        numericInput(inputId = "i2",
                                     label = "",
                                     value = 13))
               ),
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz3",
                                     label = "",
                                     value = 207.1038)),
                 column(6,
                        numericInput(inputId = "i3",
                                     label = "",
                                     value = 1))
               )
               
             ),
             
             # Output ------------------------------------
             mainPanel(
               h3("Formulas suggested:"),
               DT::dataTableOutput("table"),
               hr(),
               fluidRow(
                 column(6,
                        h3("Isotopic  pattern:"),
                        textInput("formula", label = "", 
                                  value = "C11H12N2O2"),
                        plotOutput(outputId = "isopatt")
                 ),
                 column(1),
                 column(3,
                        h3("Adducts:"),
                        fluidRow(verbatimTextOutput("pos")),
                        fluidRow(verbatimTextOutput("neg"))
                 )
               )
             ) # close main panel formula finder
           ) # close sidebar layout
  ),
  
  # Adduct calculator -----
  tabPanel("Adduct calculation",
           sidebarLayout(
             
             sidebarPanel(
               numericInput(
                 inputId = "mzX",
                 label = "m/z value:",
                 value = 166.0862,
                 step = 0.0001),
               radioButtons(
                 inputId = "adductX", 
                 label = "Adduct:",
                 choices = list("[M+H]+" = 1, "[M+K]+" = 2, 
                                "[M-H]-" = 3, "[M-H+HCOOH]-" = 4), 
                 selected = 1)
             ),
             mainPanel(
               fluidRow(column(8, verbatimTextOutput("posX"))),
               fluidRow(column(8, verbatimTextOutput("negX")))
             )),
           hr(),
           sidebarLayout(
             sidebarPanel(
               textInput("formulaX", label = "Formula:", 
                         value = "C9H11NO2")
             ),
             mainPanel(
               fluidRow(column(8, verbatimTextOutput("posX2"))),
               fluidRow(column(8, verbatimTextOutput("negX2")))
             )
           )
  ),
  
  # Instructions ----
  tabPanel("Instructions",
           h3("Under construction"),
           br(),
           br(),
           br(),
           br(),
           br(),
           hr(),
           p("Inspired by Mass Decomposition tool developed by Jan Stanstrup"),
           a("http://predret.org/tools/mass-decomposition/")
           )
)



adducts <- seq(4)
names(adducts) <- c("[M+H]+", "[M+K]+", 
                    "[M-H]-", "[M-H+HCOOH]-")



# SERVER ------
server <- function(input, output){
  data <- reactive({
    mzneutral <- getmzneut(input$mz, 
                           names(which(adducts == input$adduct)))
    
    myform <- data.frame(
      formula = getform(mzval = mzneutral, 
                        ppm = input$ppm,
                        elements = "CHNOPS"))
    
    myform$ppm <- round(abs((
      decomposeMass(mzneutral, ppm = input$ppm)$exactmass - mzneutral) / 
        decomposeMass(mzneutral, ppm = input$ppm)$exactmass) * 1e6, 1)
    
    # Rules -----
    myform$H_rule <- H_rule(myform$formula)
    myform$H_rule2 <- H_rule2(mzneutral, myform$formula)
    myform$N_rule <- N_rule(mzneutral, myform$formula)
    myform$RPU_rule <- RPU_rule(myform$formula)
    myform$iso_1 <- isotope_rule(myform$formula, 
                                 (input$i2/input$i1)*100, 
                                 range = input$error1, isotope = 1)
    myform$iso_2 <- isotope_rule(myform$formula, 
                                 (input$i3/input$i1)*100, 
                                 range = input$error2, isotope = 2)
    myform$rank <- rank_form(
      myform$formula, (input$i2/input$i1)*100, (input$i3/input$i1)*100
    )
    
    if(input$H_rule){myform <- myform[myform$H_rule, ]} 
    if(input$H_rule2){myform <- myform[myform$H_rule2, ]} 
    if(input$N_rule){myform <- myform[myform$N_rule, ]} 
    if(input$RPU_rule){myform <- myform[myform$RPU_rule, ]} 
    if(input$iso_1){myform <- myform[myform$iso_1, ]} 
    if(input$iso_2){myform <- myform[myform$iso_2, ]} 
    
    colnames(myform) <- c("Formula", "ppm", "H rule", "H rule (II)", 
                          "N rule", "RPU rule", "Isotope 1", "Isotope 2",
                          "Rank")
    myform
  })
  
  isopt <- reactive({
    formneutral <- update_form(input$formula, 
                               names(which(adducts == input$adduct)),
                               action = "add")
    
    dt <- data.frame(t(getIsotope(getMolecule(formneutral))[,1:3]))
    dt$X3 <- (dt$X2/max(dt$X2))*100
    dt
  })
  
  
  
  output$table <- DT::renderDataTable(DT::datatable({
    data()
  }, 
  options = list(pageLength = 5, dom = "tip" ),
  rownames= FALSE))
  
  output$isopatt <- renderPlot({
    plot(c(input$mz1, input$mz2, input$mz3), 
         c(100, (input$i2/input$i1)*100, (input$i3/input$i1)*100), 
         type = "h", xlab = "m/z value", ylab = "Relative intensity", 
         col = "#E31A1C", lwd = 1, ylim = c(0, 100), 
         xlim = c(min(isopt()[1,"X1"], input$mz1) - 0.3, 
                  max(isopt()[3,"X1"], input$mz3) + 0.3))
    rect(xleft = isopt()[1,"X1"] - 0.05, xright =isopt()[1,"X1"] + 0.05, 
         ybottom = 0, ytop = isopt()[1,"X3"], border = "#B2DF8A")
    rect(xleft = isopt()[2,"X1"] - 0.05, xright =isopt()[2,"X1"] + 0.05, 
         ybottom = 0, ytop = isopt()[2,"X3"], border = "#B2DF8A")
    rect(xleft = isopt()[3,"X1"] - 0.05, xright =isopt()[3,"X1"] + 0.05, 
         ybottom = 0, ytop = isopt()[3,"X3"], border = "#B2DF8A")
    legend("topright", legend = c("Theoretical", "Experimental"),
           col = c("#B2DF8A", "#E31A1C"), lty = 1, lwd = 3)
  })
  
  
  output$pos <- renderPrint({ 
    mzneutral <- getmzneut(input$mz, 
                           names(which(adducts == input$adduct)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M+H]+", "[M+Na]+"))) 
  })
  
  output$neg <- renderPrint({ 
    mzneutral <- getmzneut(input$mz, 
                           names(which(adducts == input$adduct)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M-H]-", "[M-H+HCOOH]-")))
  })
  
  output$posX <- renderPrint({ 
    mzneutral <- getmzneut(input$mzX, 
                           names(which(adducts == input$adductX)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M+H]+", "[M+Na]+", "[M+K]+"))) 
  })
  
  output$negX <- renderPrint({ 
    mzneutral <- getmzneut(input$mzX, 
                           names(which(adducts == input$adductX)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M-H]-", "[M-H+HCOOH]-")))
  })
  
  output$posX2 <- renderPrint({ 
    mzneutral <- getMolecule(input$formulaX)$exactmass
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M+H]+", "[M+Na]+", "[M+K]+"))) 
  })
  
  output$negX2 <- renderPrint({ 
    mzneutral <- getMolecule(input$formulaX)$exactmass
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M-H]-", "[M-H+HCOOH]-")))
  })
  
  
}


shinyApp(ui = ui, server = server)