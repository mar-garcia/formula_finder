options(repos = BiocManager::repositories())

library(Rdisop)
library(MetaboCoreUtils)
library(shiny)
library(CompoundDb)
library(MsCoreUtils)

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


DBE_rule <- function(formulas = character(0)){
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
  
  # DBE <- (C+Si) - (H+Cl+Fl+I)/2 + (N+P)/2 + 1
  (
    (tmp.frm1$C + tmp.frm1$Si) - 
      ((tmp.frm1$H + tmp.frm1$Cl + tmp.frm1$Fl + tmp.frm1$I)/2) + 
      ((tmp.frm1$N + tmp.frm1$P)/2) + 
      1
  ) # >= 0
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


HC_rule <- function(formulas = character(0)){
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
  
  tmp.frm1$H / tmp.frm1$C
}


NC_rule <- function(formulas = character(0)){
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
  
  tmp.frm1$N / tmp.frm1$C
}

OC_rule <- function(formulas = character(0)){
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
  
  tmp.frm1$O / tmp.frm1$C
}

PC_rule <- function(formulas = character(0)){
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
  
  tmp.frm1$P / tmp.frm1$C
}

SC_rule <- function(formulas = character(0)){
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
  
  tmp.frm1$S / tmp.frm1$C
}

mycolnames <- c("Formula", "ppm", "Isotope 1", "Isotope 2", 
                "H rule", "H rule (II)", 
                "N rule", "DBE", "DBE rule", 
                "H/C", "H/C ratio", "N/C", "N/C ratio", 
                "O/C", "O/C ratio", "P/C", "P/C ratio",
                "S/C", "S/C ratio",
                "Rank")

# UI ------------------------------
ui <- navbarPage(
  "",
  
  # Formula finder -----------
  tabPanel("Formula Finder",
           
           # Options ------------------------------------
           sidebarLayout(
             
             sidebarPanel(
               
               fluidRow(
                 column(3,
                        numericInput(inputId = "mz",
                                     label = "m/z value:",
                                     value = 205.0971,
                                     step = 0.0001)),
                 
                 column(3, 
                        numericInput(inputId = "ppm",
                                     label = "ppm:",
                                     value = 10,
                                     step = 1)),
                 column(4,
                        selectInput(
                          inputId = "adduct", 
                          label = "Adduct:",
                          choices = list("[M+H]+" = 1, "[M+Na]+" = 2, "[M+K]+" = 3, 
                                         "[M-H]-" = 4, "[M+Cl]-" = 5, "[M-H+HCOOH]-" = 6), 
                          selected = 1))),
               
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
               helpText("The sliders allows to fix an specific 'error' window when filtering for one or both isotopes."),
               
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz1",
                                     label = "m/z value",
                                     value = 205.0971)),
                 column(6,
                        numericInput(inputId = "i1",
                                     label = "intensity",
                                     value = 4854167))
               ),
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz2",
                                     label = "",
                                     value = 206.1005)),
                 column(6,
                        numericInput(inputId = "i2",
                                     label = "",
                                     value = 631041))
               ),
               fluidRow(
                 column(6,
                        numericInput(inputId = "mz3",
                                     label = "",
                                     value = 207.1038)),
                 column(6,
                        numericInput(inputId = "i3",
                                     label = "",
                                     value = 48542))
               ),
               
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
               
               checkboxInput("DBE_rule", 
                             label = "DBE rule", 
                             value = TRUE),
               fluidRow(
                 column(3,
                        checkboxInput("HC_ratio", 
                                      label = "H/C ratio", 
                                      value = TRUE)),
                 column(6,
                        sliderInput("HC_ratio_range", label = "", 
                                    min = 0, max = 10, value = c(0.2, 5), step = 0.1))
               ),
               fluidRow(
                 column(3,
                        checkboxInput("NC_ratio", 
                                      label = "N/C ratio", 
                                      value = TRUE)),
                 column(6,
                        sliderInput("NC_ratio_range", label = "", 
                                    min = 0, max = 10, value = c(0, 3), step = 0.1))
               ),
               fluidRow(
                 column(3,
                        checkboxInput("OC_ratio", 
                                      label = "O/C ratio", 
                                      value = TRUE)),
                 column(6,
                        sliderInput("OC_ratio_range", label = "", 
                                    min = 0, max = 10, value = c(0, 5), step = 0.1))
               ),
               fluidRow(
                 column(3,
                        checkboxInput("PC_ratio", 
                                      label = "P/C ratio", 
                                      value = TRUE)),
                 column(6,
                        sliderInput("PC_ratio_range", label = "", 
                                    min = 0, max = 10, value = c(0, 1.5), step = 0.1))
               ),
               fluidRow(
                 column(3,
                        checkboxInput("SC_ratio", 
                                      label = "S/C ratio", 
                                      value = TRUE)),
                 column(6,
                        sliderInput("SC_ratio_range", label = "", 
                                    min = 0, max = 10, value = c(0, 3), step = 0.1))
               ),
               
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
                 column(4,
                        h3("Adducts:"),
                        fluidRow(verbatimTextOutput("pos")),
                        fluidRow(verbatimTextOutput("neg"))
                 )
               ),
               hr(),
               fluidRow(
                 checkboxGroupInput(
                   "show_vars", "Columns to show in the table:", mycolnames, 
                   selected = c("Formula", "ppm", "Isotope 1", "Isotope 2", 
                                "H rule", "H rule (II)", "N rule", 
                                "DBE", "H/C", "N/C", "O/C", "P/C", "S/C", 
                                "Rank"), 
                   inline = TRUE)
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
                 choices = list("[M+H]+" = 1, "[M+Na]+" = 2, "[M+K]+" = 3, 
                                "[M-H]-" = 4, "[M+Cl]-" = 5, "[M-H+HCOOH]-" = 6), 
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
               fluidRow(column(8, verbatimTextOutput("neutral"))),
               fluidRow(column(8, verbatimTextOutput("posX2"))),
               fluidRow(column(8, verbatimTextOutput("negX2")))
             )
           ),
           hr(),
           sidebarLayout(
             sidebarPanel(
               numericInput(
                 inputId = "mzX1",
                 label = "m/z value 1:",
                 value = 166.0862,
                 step = 0.0001),
               numericInput(
                 inputId = "mzX2",
                 label = "m/z value 2:",
                 value = 120.0807,
                 step = 0.0001),
               numericInput(inputId = "ppmX",
                            label = "ppm:",
                            value = 10,
                            step = 1)
               ),
             mainPanel(
               fluidRow(column(8, verbatimTextOutput("adductZ")))
             )
           )
  ),
  
  # Instructions ----
  tabPanel("Instructions",
           p("This shiny app contain 2 main tools: 'Formula Finder' and 'Adduct calculation'."),
           p("If you have any question and/or any suggestion to improve this shiny app, please write me at mar.garcia@eurac.edu. Your feedback will be very appreciated! :)"),
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
  )
)



adducts <- seq(6)
names(adducts) <- c("[M+H]+", "[M+Na]+", "[M+K]+", 
                    "[M-H]-", "[M+Cl]-", "[M-H+HCOOH]-")



# SERVER ------
server <- function(input, output){
  
  # reactive ----------
  
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
    
    myform$iso_1 <- isotope_rule(myform$formula, 
                                 (input$i2/input$i1)*100, 
                                 range = input$error1, isotope = 1)
    myform$iso_2 <- isotope_rule(myform$formula, 
                                 (input$i3/input$i1)*100, 
                                 range = input$error2, isotope = 2)
    myform$H_rule <- H_rule(myform$formula)
    myform$H_rule2 <- H_rule2(mzneutral, myform$formula)
    myform$N_rule <- N_rule(mzneutral, myform$formula)
    myform$DBE <- DBE_rule(myform$formula)
    myform$DBE_rule <- (myform$DBE %%1 == 0) & (myform$DBE >= 0) 
    
    myform$HC <- round(HC_rule(myform$formula), 1)
    myform$HC_ratio <- (myform$HC >= input$HC_ratio_range[1]) &
      (myform$HC < input$HC_ratio_range[2])
    myform$NC <- round(NC_rule(myform$formula), 1)
    myform$NC_ratio <- (myform$NC >= input$NC_ratio_range[1]) &
      (myform$NC <= input$NC_ratio_range[2])
    myform$OC <- round(OC_rule(myform$formula), 1)
    myform$OC_ratio <- (myform$OC >= input$OC_ratio_range[1]) &
      (myform$OC <= input$OC_ratio_range[2])
    myform$PC <- round(PC_rule(myform$formula), 1)
    myform$PC_ratio <- (myform$PC >= input$PC_ratio_range[1]) &
      (myform$PC <= input$PC_ratio_range[2])
    myform$SC <- round(SC_rule(myform$formula), 1)
    myform$SC_ratio <- (myform$SC >= input$SC_ratio_range[1]) &
      (myform$SC <= input$SC_ratio_range[2])
    
    if(input$iso_1){myform <- myform[myform$iso_1, ]} 
    if(input$iso_2){myform <- myform[myform$iso_2, ]} 
    if(input$H_rule){myform <- myform[myform$H_rule, ]} 
    if(input$H_rule2){myform <- myform[myform$H_rule2, ]} 
    if(input$N_rule){myform <- myform[myform$N_rule, ]} 
    if(input$DBE_rule){myform <- myform[myform$DBE_rule, ]} 
    if(input$HC_ratio){myform <- myform[myform$HC_ratio, ]} 
    if(input$NC_ratio){myform <- myform[myform$NC_ratio, ]} 
    if(input$OC_ratio){myform <- myform[myform$OC_ratio, ]} 
    if(input$PC_ratio){myform <- myform[myform$PC_ratio, ]} 
    if(input$SC_ratio){myform <- myform[myform$SC_ratio, ]} 
    
    myform$rank <- rank_form(
      myform$formula, (input$i2/input$i1)*100, (input$i3/input$i1)*100
    )
    
    colnames(myform) <- mycolnames
    myform[, input$show_vars]
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
         col = "#E31A1C", lwd = 1, xaxt = 'n', ylim = c(0, 100), 
         xlim = c(min(isopt()[1,"X1"], input$mz1) - 0.3, 
                  max(isopt()[3,"X1"], input$mz3) + 0.3))
    axis(1, at = isopt()[,"X1"], 
         labels = c(paste0(round(isopt()[1,"X1"],4), " (", 
                           round((abs(isopt()[1,"X1"] - input$mz1)/isopt()[1,"X1"])*1e6, 1), " ppm)"),
                    paste0(round(isopt()[2,"X1"],4), " (", 
                           round((abs(isopt()[2,"X1"] - input$mz2)/isopt()[2,"X1"])*1e6, 1), " ppm)"),
                    paste0(round(isopt()[3,"X1"],4), " (", 
                           round((abs(isopt()[3,"X1"] - input$mz3)/isopt()[3,"X1"])*1e6, 1), " ppm)")
         ))
    rect(xleft = isopt()[,"X1"] - (ppm(isopt()[,"X1"], 100)), 
         xright =isopt()[,"X1"] + (ppm(isopt()[,"X1"], 100)), 
         ybottom = 0, ytop = isopt()[,"X3"], border = "#B2DF8A")
    #rect(xleft = isopt()[2,"X1"] - 0.05, xright =isopt()[2,"X1"] + 0.05, 
     #    ybottom = 0, ytop = isopt()[2,"X3"], border = "#B2DF8A")
    #rect(xleft = isopt()[3,"X1"] - 0.05, xright =isopt()[3,"X1"] + 0.05, 
     #    ybottom = 0, ytop = isopt()[3,"X3"], border = "#B2DF8A")
    legend("topright", legend = c("Theoretical", "Experimental"),
           col = c("#B2DF8A", "#E31A1C"), lty = 1, lwd = 3)
  })
  
  
  output$pos <- renderPrint({ 
    mzneutral <- getmzneut(input$mz, 
                           names(which(adducts == input$adduct)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[2M+H]+"))) 
  })
  
  output$neg <- renderPrint({ 
    mzneutral <- getmzneut(input$mz, 
                           names(which(adducts == input$adduct)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M-H]-", "[M+Cl]-", "[M-H+HCOOH]-", "[2M-H]-")))
  })
  
  output$posX <- renderPrint({ 
    mzneutral <- getmzneut(input$mzX, 
                           names(which(adducts == input$adductX)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[2M+H]+"))) 
  })
  
  output$negX <- renderPrint({ 
    mzneutral <- getmzneut(input$mzX, 
                           names(which(adducts == input$adductX)))
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M-H]-", "[M+Cl]-", "[M-H+HCOOH]-", "[2M-H]-")))
  })
  
  output$neutral <- renderPrint({ 
    paste("Monoisotopic mass:", getMolecule(input$formulaX)$exactmass)
  })
  
  output$posX2 <- renderPrint({ 
    mzneutral <- getMolecule(input$formulaX)$exactmass
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[2M+H]+"))) 
  })
  
  output$negX2 <- renderPrint({ 
    mzneutral <- getMolecule(input$formulaX)$exactmass
    unlist(CompoundDb::mass2mz(mzneutral, 
                               adduct = c("[M-H]-", "[M+Cl]-", "[M-H+HCOOH]-", "[2M-H]-")))
  })
  
  output$adductZ <- renderPrint({ 
    tmp <- unlist(CompoundDb::mass2mz(-1.007276, adduct = adducts()))
    paste0(names(unlist(matchWithPpm(abs(input$mzX1 - input$mzX2), 
                                     abs(tmp), ppm = input$ppmX))),
           ": ",
           round((abs(
             abs(input$mzX1 - input$mzX2) - 
               abs(tmp[unlist(matchWithPpm(abs(input$mzX1 - input$mzX2), 
                                           abs(tmp), ppm = input$ppmX))])) / 
               abs(tmp[unlist(matchWithPpm(abs(input$mzX1 - input$mzX2), 
                                           abs(tmp), ppm = input$ppmX))])
             )*1e6,1),
           " ppm"
    )
  })
  
  
}


shinyApp(ui = ui, server = server)