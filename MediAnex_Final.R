# Establishing connection between stored data and network template

# ------> Loading relevant packages: ----

if (!require('shiny')) install.packages('shiny'); library('shiny')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('plotly')) install.packages('plotly'); library('plotly')
if (!require('corrplot')) install.packages('plotly'); library('corrplot')
if (!require("BiocManager")) install.packages("BiocManager"); library('BiocManager')
if (!require('mixOmics')) BiocManager::install('mixOmics'); library('mixOmics')
if (!require('ropls')) BiocManager::install('ropls'); library('ropls')
if (!require('DT')) install.packages('DT'); library('DT')
if (!require('MVN')) install.packages('MVN'); library('MVN')
if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('snowfall')) install.packages('snowfall')
if (!require('neldermead')) install.packages('neldermead')
if (!require('optimbase')) install.packages('optimbase')
if (!require('classyfire')) install.packages("classyfire_0.1-2.tar.gz", repos = NULL, type = "source"); library('classyfire')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('arm')) install.packages('arm'); library('arm')
if (!require('xgboost')) install.packages('xgboost'); library('xgboost')
if (!require('DiagrammeR')) install.packages('DiagrammeR'); library('DiagrammeR')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('Ckmeans.1d.dp')) install.packages('Ckmeans.1d.dp'); library('Ckmeans.1d.dp')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('tibble')) install.packages('tibble'); library('tibble')
if (!require('pROC')) install.packages('pROC'); library('pROC')
if (!require('shinythemes')) install.packages('shinythemes');library('shinythemes')
if (!require('bslib')) install.packages('bslib'); library('bslib')
if (!require('shinyWidgets')) install.packages('shinyWidgets'); library('shinyWidgets')
if (!require('shinydashboard')) install.packages('shinydashboard'); library('shinydashboard')
if (!require('colourpicker')) install.packages('colourpicker'); library('colourpicker')
if (!require('shinyalert')) install.packages('shinyalert'); library('shinyalert')
if (!require('reporter')) install.packages('reporter'); library('reporter')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('tidygraph')) install.packages('tidygraph'); library('tidygraph')
if (!require('visNetwork')) install.packages('visNetwork'); library('visNetwork')
if (!require('rintrojs')) install.packages('rintrojs'); library('rintrojs')
if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('pryr')) install.packages('pryr'); library('pryr')
options(scipen = 999)


# Loading the database data: ----

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#standardized naming for all files 
###Correct the name to be in the correct format 
Rename_Met_Table<-read.csv(file = "Coefficients/Rename_Met2.csv", header = TRUE, sep = ",")
colnames(Rename_Met_Table) <- c("Final_Name", "Other_Name")
#Create a dictionary with values as the old names and the rename as the values

Rename_Met<-setNames(Rename_Met_Table$Final_Name, Rename_Met_Table$Other_Name)

# #dictionary for SPMS
SPM_FA_Table<- read.csv(file = "Coefficients/SPM_FA.csv", header = TRUE, sep = ",")
colnames(SPM_FA_Table) <- c("Metabolite", "Fatty_Acid")

# Example Files (PLS-DA, Machine Learning & Network Analysis):
Test_data <- read.csv(file = "TestFiles/TestData.csv", header = TRUE, sep = ",")
Train_Data <- read.csv(file = "TestFiles/Train_RA.csv", header = TRUE, sep = ",")
Train_RML <- read.csv(file = "TestFiles/Train_RA2.csv", header = TRUE, sep = ",")
Network_Data <- read.csv(file = "TestFiles/Net_Data.csv", header = TRUE, sep = ",")

# ML Model time taken:
format_time <- function(time_diff) {
  minutes <- floor(as.numeric(time_diff) / 60)
  seconds <- round(as.numeric(time_diff) %% 60, 2)
  paste0(minutes, " minute(s) ", seconds, " second(s)")
}

SPMs_FA<-setNames(SPM_FA_Table$Fatty_Acid, SPM_FA_Table$Metabolite)

# Loading EPA Network template: ----
EPA_data <- read.csv(file = "Network_templates/EPA_network.csv", header = TRUE, sep = ",")
EPA_direction <- read.csv(file = "Network_templates/EPA_edges.csv", header = TRUE, sep = ",")

# Loading AA Network template: ----
AA_data <- read.csv(file = "Network_templates/AA_network.csv", header = TRUE, sep = ",")
AA_direction <- read.csv(file = "Network_templates/AA_edges.csv", header = TRUE, sep = ",")

# Loading the DHA Network template: ----
DHA_data <- read.csv(file = "Network_templates/DHA_network.csv", header = TRUE, sep = ",")
DHA_direction <- read.csv(file = "Network_templates/DHA_edges.csv", header = TRUE, sep = ",")

# Loading the DPA Network template: ----
DPA_data <- read.csv(file = "Network_templates/DPA_network.csv", header = TRUE, sep = ",")
DPA_direction <- read.csv(file = "Network_templates/DPA_edges.csv", header = TRUE, sep = ",")

# Building a dictionary which will be used to add subscripts to the SPM names: ----
SPM_sub <- c("LTC4S/GSTM4" = "LTC" %p% subsc('4') %p% "S/GSTM4", "PD1n-3DPA" = "PD1" %p% subsc('n-3') %p% "DPA", 
             "PD2n-3DPA" = "PD2" %p% subsc('n-3') %p% "DPA", "RvD1n-3DPA" = "RvD1" %p% subsc('n-3') %p% "DPA", 
             "RvD2n-3DPA" = "RvD2" %p% subsc('n-3') %p% "DPA", "RvD5n-3DPA" = "RvD5" %p% subsc('n-3') %p% "DPA", 
             "MaR1n-3DPA" = "MaR1" %p% subsc('n-3') %p% "DPA", "MaR2n-3DPA" = "MaR2" %p% subsc('n-3') %p% "DPA", 
             "TXB2" = "TXB" %p% subsc('2'), "PGD2" = "PGD" %p% subsc('2'), "PGE2" = "PGE" %p% subsc('2'), 
             "PGF2a" = "PGF" %p% subsc('2a'), "15-epi-LXA4" = "15-epi-LXA" %p% subsc('4'), 
             "15-epi-LXB4" = "15-epi-LXB" %p% subsc('4'), "LXA4" = "LXA" %p% subsc('4'), 
             "LXB4" = "LXB" %p% subsc('4'), "15-oxo-LXA4" = "15-oxo-LXA" %p% subsc('4'), 
             "13,14-dehydro-15-oxo LXA4" = "13,14-dehydro-15-oxo LXA" %p% subsc('4'), 
             "6-trans-12-epi-LTB4" = "6-trans-12-epi-LTB" %p% subsc('4'), "6-trans-LTB4" = "6-trans-LTB" %p% subsc('4'),
             "LTA4H" = "LTA" %p% subsc('4') %p% "H", "LTB4" = "LTB" %p% subsc('4'), "LTD4" = "LTD" %p% subsc('4'), 
             "20-OH-LTB4" = "20-OH-LTB" %p% subsc('4'), "20-COOH-LTB4" = "20-COOH-LTB" %p% subsc('4'), 
             "LTC4S" = "LTC" %p% subsc('4') %p% "S", "LTC4" = "LTC" %p% subsc('4'), "LTE4" = "LTE" %p% subsc('4'))

# Using dictionary to add subscripts to SPM names:
EPA_data <- EPA_data %>% mutate(label = recode(label, !!!SPM_sub))

AA_data <- AA_data %>% mutate(label = recode(label, !!!SPM_sub))

DPA_data <- DPA_data %>% mutate(label = recode(label, !!!SPM_sub))

DHA_data <- DHA_data %>% mutate(label = recode(label, !!!SPM_sub))

# Ui section: ----

ui <- navbarPage("MediANEX",
                 
                 theme = bslib::bs_theme(
                   primary = "#008CBA",
                   secondary = "#C6C6C6",
                   heading_font = "Helvetica",
                   base_font = "Helvetica", 
                   'enable-gradients' = TRUE,
                   'enable-shadows' = TRUE,
                   preset = "yeti"
                 ),
                 
                 tags$head(
                   tags$style(HTML("
                    .navbar-brand {
                      font-family: 'Helvetica', sans-serif; /* Choose a futuristic font */
                      font-size: 2em; /* Adjust font size as needed */
                      color: #FFFFFF ; /* Adjust text color as needed */
                    }
                  "))),
                 
                 navbarMenu(
                   title = tags$span("Multivariate Statistics", id = "multivar_stats_nav"),
                   tabPanel(
                     title = "PCA/PLS-DA",
                     id = "PCA/PLS-DA",
                     titlePanel("PCA/PLS-DA Analysis"),
                     page_sidebar(
                       sidebar = sidebar(
                       #size of sidebar
                       width = 250,
                       
                       actionButton("Help_PCA", "Tutorial"),
                       
                       tags$hr(),
                       
                       div(id = "PCA",
                           #File upload code
                           fileInput("PCAPLSDA_file", "Choose upload file",
                                     multiple = FALSE,
                                     accept = c("text/csv", 
                                                "text/comma-separated-values,text/plain",
                                                ".csv")),
                           #Check box to see if file has headers 
                           checkboxInput("header_PCA", "Header", TRUE),
                           
                           #input--Separator for the file 
                           radioButtons("sep_PCA", "Separator", 
                                        choices = c(Comma = ",",
                                                    Tab = "\t")),
                           
                           #input---Location of samples 
                           radioButtons("samp_PCA", "Sample location",
                                        choices = c(Row = "row_PCA", 
                                                    Column = "col_PCA"),
                                        selected = "row_PCA")
                       ),
                       
                       div(id = "PCA_col",
                         #input -- where the sample names column is
                         #numericInput("SampCol_PCA", "Sample name Column/Row Number", 1, min = 1, max = NA),
                         selectInput("SampleCol_PCA", "Select Sample Column Name", "Please upload file"),
                         helpText("If samples were the rownames, please select x."),
                         #input -- where the sample group info column is
                         selectInput("Group_PCA", "Select Group Name", "Please upload file"),
                         downloadButton("downTest", "Download Example File"),
                       ),
                       
                       
                       
                       #line separator 
                       tags$hr(),
                       
                       #Input --- Selecting the test
                       radioButtons("test_PCAPLSDA", "Select test",
                                    choices = c(PCA = "PCA",
                                                PLSDA = "PLS-DA"),
                                    selected = "PCA"),
                       
                       #line separator
                       tags$hr(),
                       
                       #calculate button 
                       actionButton("Go_PCAPLSDA", "Calculate")
                     ),
                     navset_bar(
                       fillable = TRUE,
                       fluid = TRUE,
                       gap = "0px",
                       padding = "0px",
                       nav_panel(
                         title = "2D Scores Plot",
                         value = "",
                                 uiOutput("error_PLSDA"),
                                 uiOutput("error_PLSDA2"),
                                 uiOutput("error_PLSDA3"),
                                 uiOutput("error_PLSDA4"),
                                 uiOutput("error_PLSDA_group_select"),
                                 plotlyOutput("Xscoresplot"),
                                 #Check box to show samples 
                                 checkboxInput("ViewSamples_PCA", "View Samples", TRUE),
                                 #Update information 
                                 actionButton("Go_PCAPLSDA2", "Update Plot"),
                                 downloadButton("DownloadScore", "Download X-scores Data"),
                                     
                                 
                                 bslib::accordion(id = "R2Q2",
                                                  open = FALSE,
                                                  
                                                  bslib::accordion_panel(title = "R2Q2 Plot",
                                                                         value = "",
                                                                         uiOutput("R2Q2_Title"),
                                                                         dataTableOutput("R2Q2_Table"),
                                                                         plotlyOutput("R2Q2_Plot")))
                       ),
                       nav_panel("Loading Plot", 
                                 plotOutput("Loadingplot"),
                                 downloadButton("downloadLoading_Plot", "Download Loading Plot"),
                                 downloadButton("DownloadLoading", "Download Loading Score Data")),
                       nav_panel("Variance Plot", 
                                 plotlyOutput("PerVarplot", height = "600px"),
                                 downloadButton("DownloadVar", "Download Variance Data")),
                       nav_panel(
                         title = "VIP Scores", 
                         id = "VIPScore",
                         
                         textOutput("VIPtext"),
                         plotOutput("VIPplot"), 
                         downloadButton("downloadPLSDA_VIP_Plot", "Download VIP plot"),
                         downloadButton("DownloadVIP", "Download VIP Table"),
                         
                                 bslib::accordion(id = "VIP_path", 
                                                  open = FALSE,
                                                  bslib::accordion_panel(
                                                    title = "VIP Pathway Plot",
                                                    value = "",
                                                    uiOutput("VIP_Pathway_Title"),
                                                    plotlyOutput("VIP_PathwayPlot")))
                         
                     ))
                     )),
                   
                   tabPanel(title ="Differential Analysis",
                            id = "Diff",
                            titlePanel("Differential Analysis"),
                            
                            page_sidebar(
                              sidebar = sidebar(
                                introjsUI(),
                                
                                #size of sidebar
                                width = 250,
                                
                                actionButton("Help_Diff", "Tutorial"),
                                
                                div(id = "Diff_tut",
                                    #File upload code 
                                    fileInput("lm_profiles", "Choose file to upload",
                                              multiple = FALSE,
                                              accept = c("text/csv", 
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    
                                    #Check box to see if file has headers 
                                    checkboxInput("header_Diff", "Header", TRUE),
                                    
                                    #input--Separator for the file 
                                    radioButtons("sep_Diff", "Separator", 
                                                 choices = c(Comma = ",",
                                                             Tab = "\t")),
                                    
                                    
                                    
                                    tags$hr(),
                                    
                                    downloadButton("downDiff", "Download Example File"),
                                    #line separator 
                                    tags$hr()
                                ),
                                
                                div(id = "normality",
                                    # Input: Select separator ----
                                    radioButtons("mvn", "Multivariate Normality Test",
                                                 choices = c(Mardia = "mardia",
                                                             `Henze-Zirkler` = "hz",
                                                             Royston = "royston"),
                                                 selected = "royston"),
                                    
                                    #line separator
                                    tags$hr(),
                                    
                                    # Group selection ----
                                    selectInput("group_A", "Group A",""),
                                    selectInput("group_B", "Group B",""),
                                    
                                    
                                    #calculate button 
                                    actionButton("Go_Diff", "Calculate")
                                    
                                    # #refresh button
                                    # actionButton("reset_diff", "Reset inputs")
                                ),
                                
                              
                            ),
                            
                            navset_bar(
                              
                              fillable = TRUE,
                              fluid = TRUE,
                              gap = "0px",
                              padding = "0px",
                              
                              nav_panel("Results Table",
                                       dataTableOutput("DiffTable"),
                                      uiOutput("error_Diff_format"),
                                      uiOutput("error_Diff_groups"),
                                      uiOutput("error_Diff_numeric"),
                                      uiOutput("error_Diff_Snames"),
                                       downloadButton("DownloadDiffTable", "Download Data")),
                              nav_panel("Volcano Plots",
                                        bslib::accordion(
                                          id = "Volcano_options",
                                          open = FALSE,
                                          bslib::accordion_panel(
                                            title = "Volcano plot options",
                                            value = "",
                                            numericInput("Sig_DiffPlot", "Adjusted P-Value cut-off", 0.05, min = 0, max = NA),
                                            numericInput("FC_DiffPlot", "Fold Change cut-off", 1, min = 0, max = NA),
                                            numericInput("TopN_DiffPlot", "Number of Top points to Label", 10, min = 0, max = NA),
                                            actionButton("Go_DiffPlot", "Update Plot")
                                          )
                                        ),
                                       plotOutput("DiffPlot", height = "600px"),
                                       downloadButton("downloadDiffPlot", "Download Volcano Plot")),
                              nav_panel("Pathways",
                                        
                                        bslib::accordion(id = "Path_analysis",
                                                         open = TRUE,
                                                         bslib::accordion_panel(
                                                           title = "Pathways showing Upregulation",
                                                           value = "",
                                                           uiOutput("Diff_UpRegPath_Title"),
                                                           plotlyOutput("Diff_UpRegPath_Plot")
                                                           )),
                                        bslib::accordion(id = "Path_analysis2",
                                                         open = FALSE,
                                                         bslib::accordion_panel(
                                                           title = "Pathways showing downregulation",
                                                           value = "",
                                                           uiOutput("Diff_DownRegPath_Title"),
                                                           plotlyOutput("Diff_DownRegPath_Plot")
                                                         ),
                                                         bslib::accordion_panel(
                                                           title = "Number of Significant SPMs per Pathway",
                                                           value = "",
                                                           uiOutput("Diff_PValPath_Title"),
                                                           plotlyOutput("Diff_PValPath_Plot")
                                                         )
                                                         ),
                              nav_panel("Network Analysis"))         
                            ))),
                   
                   tabPanel(
                     "Network analysis",
                     
                       titlePanel("Network analysis"),
                       
                       page_sidebar(
                         sidebar = sidebar(
                           introjsUI(),
                           width = 250,
                           actionButton("Help_Net", "Tutorial"),
                           br(),
                           "Upload SPM data file (tsv or txt) or csv file.",
                           
                           div(id = "SPM",
                               fileInput("SPM_val", "SPM values",
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv", ".tsv")),
                               tags$hr(),
                               
                               # Input: Checkbox if file has header ----
                               checkboxInput("header", "Header", TRUE),
                               
                               # Input: Select separator ----
                               radioButtons("sep", "Separator",
                                            choices = c(Comma = ",",
                                                        Semicolon = ";",
                                                        Tab = "\t"),
                                            selected = ","),
                               
                               # Horizontal line ----
                               tags$hr(),
                               
                               # Input: Select number of rows to display ----
                               radioButtons("disp", "Display",
                                            choices = c("Adjusted P values" = "pval",
                                                        "VIP scores" = "vip.score"),
                                            selected = "pval")
                           ),
                           
                           downloadButton("Net_test", "Download Example File"),
                           tags$br(),
                           actionButton("info_button", "More information", class = "btn btn-primary"),
                           
                           tags$script(HTML("
                                                            $(document).ready(function(){
                                                            $('#info_button').tooltip({
                                                            title: 'In this part of the options section you may select the parameter of significance used in your file.' +
                                                            ' If SPM regulation significance is assesed using adjusted p values, and your file contains this column,' + 
                                                            ' then the adjusted p value should be selected.', 
                                                            placement: 'top', 
                                                            trigger: 'hover'
                                                            });
                                                            });
                                                                             ")),
                           style = "info",
                           
                           # Horizontal line ----
                           tags$hr(),
                           
                           div(id = "Net_select",
                               selectInput("selection",
                                           label = "Network Selection",
                                           choices = c("No Selection", 
                                                       "EPA Network",
                                                       "AA Network",
                                                       "DHA Network",
                                                       "n-3 DPA Network"),
                                           selected = "No Selection"),
                               
                               actionButton("Go","Generate Network"),
                           ),
                           
                           bslib::accordion(id = "A1",
                                            
                                            open = TRUE,
                                            
                                            bslib::accordion_panel(
                                              
                                              title = "Network Colour Options",
                                              div(id = "Color_options",
                                                  "Select Node Color",
                                                  colourInput("col1", "Select Colour", value = "lightblue"),
                                                  actionButton("update_node_color", "Update Node Color"),
                                                  tags$hr(),
                                                  "Select Edge color",
                                                  colourInput("col2", "Select Colour", value = "lightblue"),
                                                  actionButton("update_edge_color", "Update Edge Color")
                                              )
                                              )
                                            )
                           ),
                         
                         navset_bar(
                           
                           fillable = TRUE,
                           fluid = TRUE,
                           gap = "0px",
                           padding = "0px",
                           
                           nav_panel(title = "Human Model",
                                     value = "",
                                     # Main panel for displaying outputs ----
                                     bslib::accordion(
                                       id = "A2",
                                       bslib::accordion_panel(
                                         title = "Pathway analysis",
                                         card(fluidPage(visNetworkOutput(
                                           "network", 
                                           width = "100%", 
                                           height = "1000px"),
                                           height = "1400px", width = "100%")))),
                                     
                                     bslib::accordion(
                                       id = "A3",
                                       open = FALSE,
                                       bslib::accordion_panel(
                                         title = "Unrecognised SPMs",
                                         HTML("This panel is showing any SPM names from the uploaded file which have not been recognised.
                                         Please check the names which have not been matched from the databasse and ensure that the names match on your file. 
                                         (<u><b>Please Note:</b></u> Some SPMs may have different names if the nomenclature is different. For example, PDX was previously known as 10S,17S-diHDHA.)"),
                                         tags$hr(),
                                         dataTableOutput("unused_database"),
                                         tags$hr(),
                                         dataTableOutput("unused_file")))
                                     ),
                           nav_panel(title = "Mouse model",
                                     value = "",
                                     # Main panel for displaying outputs ----
                                     bslib::accordion(
                                       id = "A2",
                                       bslib::accordion_panel(
                                         title = "Pathway analysis",
                                         card(fluidPage(visNetworkOutput(
                                           "network_mouse", 
                                           width = "100%", 
                                           height = "1000px"),
                                           height = "1400px", width = "100%")))),
                                     
                                     bslib::accordion(
                                       id = "A3",
                                       open = FALSE,
                                       bslib::accordion_panel(
                                         title = "Unrecognised SPMs",
                                         HTML("This panel is showing any SPM names from the uploaded file which have not been recognised.
                                         Please check the names which have not been matched from the databasse and ensure that the names match on your file. 
                                         (<u><b>Please Note:</b></u> Some SPMs may have different names if the nomenclature is different. For example, PDX was previously known as 10S,17S-diHDHA.)"),
                                         tags$hr(),
                                         dataTableOutput("unused_database2"),
                                         tags$hr(),
                                         dataTableOutput("unused_file2")))
                             
                           )
                         
                         
                         )
                         )
                       
                     )
                   ),
                 
                 navbarMenu(
                   "Machine Learning Models",
                   tabPanel(
                     title ="Machine Learning ",
                     id = "ML",
                     titlePanel("Machine Learning"),
                     
                     # Sidebar panel for inputs ----
                     page_sidebar(
                       sidebar = sidebar(
                         introjsUI(),
                         width = 270,
                         
                         actionButton("Help_ML", "Tutorial"),
                         
                         div(id = "ML_input",
                             fileInput(
                               "LM_ML_File", 
                               "Lipid Mediator profiling file:",
                               multiple = FALSE,
                               accept = c("text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv", ".tsv", ".txt")),
                             
                             # Input: Select separator ----
                             radioButtons(
                               "sep_LM_ML", "Separator",
                               choices = c(Comma = ",",
                                           Semicolon = ";",
                                           Tab = "\t"),
                               selected = ",")
                           
                         ),
                         
                         div(id = "Select_ML",
                             #input -- where the sample group info column is
                             selectInput("Group_ML", "Select Group Column Name", "Please upload file"),
                         ),       
                                
                         # Download button for the example file ----
                         downloadButton(outputId = "TrainTest", label = "Example ML file"),
                                
                         actionButton("info_button2", "More information", class = "btn btn-primary"),
                         
                         tags$script(HTML("
                         $(document).ready(function(){
                         $('#info_button2').tooltip({
                         title: 'The ML profiling file consist in a table with samples per row and each lipid mediator as columns.' + 
                         ' An extra column is added called <b><u> groups </u></b>, specifying to which group every sample belongs to;' + 
                         ' and an extra row, added after the lipid mediators names, meaning the <b><u> second row </u></b>' +
                         ' which contains the fatty acid substrates to which every lipid mediator comes from.' +
                         ' You can see the format dowloading the example file.',
                         placement: 'top', 
                         trigger: 'hover',
                         html: true,
                         });
                         });
                                          ")),
                         style = "info",
                         
                         # Horizontal line ----
                         tags$hr(),
                            
                         # Horizontal line ----
                         tags$hr(),
                         
                         div(id = "Gen_ML",
                             checkboxGroupInput("Model_ML", "Select Machine Learning Model",
                                                c("Random Forests" = "RF_ML", 
                                                  "Extreme Gradient Boosting" = "XGB_ML",
                                                  "Support Vector Machine" = "SVM_ML",
                                                  "Elastic Net Regresion" = "LA_ML",
                                                  "Bayesian Linear Model" = "BC_ML")),
                             
                             # Action bottom to create and run the ML models ---
                             actionButton("Go_ML","Create Machine Learning Models")
                           
                         )
                         
                         
                         ),
                       
                       # Main panel for displaying outputs ----
                       
                       navset_bar(
                         fillable = TRUE,
                         fluid = TRUE,
                         gap = "0px",
                         padding = "0px",
                         
                         nav_panel("Accuracy",
                                   
                                   bslib::accordion(
                                     id = "Acc_1 ",
                                     open = TRUE,
                                     bslib::accordion_panel(
                                       title = "Accuracy scores plot",
                                       value = "",
                                       uiOutput("AccPlot_Title"),
                                       tags$hr(),
                                       plotOutput("Accuracy_ML"),
                                       uiOutput("error_ML_coln"),
                                       uiOutput("error_ML_Snames"),
                                       uiOutput("error_ML_groups"),
                                       htmlOutput("system_info"),
                                       downloadButton("downloadAcc_ML_Plot", "Download Accuracy Plot"))),
                                   
                                   bslib::accordion(
                                     id = "Acc_2",
                                     open = FALSE,
                                     bslib::accordion_panel(
                                       title = "Machine learning model table",
                                       value = "",
                                       uiOutput("AccTable_Title"),
                                       tags$hr(),
                                       dataTableOutput("ML_Table"),
                                       downloadButton("downloadAcc_ML_Table", "Download Accuracy Table")))),
                         
                         nav_panel("Random Forests",
                                   navset_card_tab(
                                     id = "RF_card1",
                                     nav_panel("Parameters",
                                               uiOutput("RF_Plot_Title"),
                                               tags$hr(),
                                               bslib::accordion(
                                                 id = "RF_1",
                                                 open = TRUE,
                                                 bslib::accordion_panel(
                                                   title = "Average Percent Accuracy Plot of All ML Models",
                                                   value = "",
                                                   plotOutput("RF_Plot1", height = "500px"),
                                                   tags$hr()
                                                 )
                                               ),
                                               bslib::accordion(
                                                 id = "RF_2",
                                                 open = FALSE,
                                                 bslib::accordion_panel(
                                                   title = "Average Percent Accuracy Plot of DHA ML Model",
                                                   value = "",
                                                   plotOutput("RF_Plot2", height = "500px"),
                                                   tags$hr()
                                                 ), 
                                                 
                                                 bslib::accordion_panel(
                                                   title = "Average Percent Accuracy Plot of n-3 DPA ML Model",
                                                   value = "",
                                                   plotOutput("RF_Plot3", height = "500px"),
                                                   tags$hr()
                                                 ),
                                                 
                                                 bslib::accordion_panel(
                                                   title = "Average Percent Accuracy Plot of EPA ML Model",
                                                   value = "",
                                                   plotOutput("RF_Plot4", height = "500px"),
                                                   tags$hr()
                                                 ),
                                                 
                                                 bslib::accordion_panel(
                                                   title = "Average Percent Accuracy Plot of AA ML Model",
                                                   value = "",
                                                   plotOutput("RF_Plot5", height = "500px"),
                                                   tags$hr()
                                                 )
                                               ),
                                               
                                               selectInput("RF_Mod_Num1", "Select which model you want to download", ""),
                                               
                                               downloadButton("download_RF_Plot", "Download Random Forest Plot"),
                                               downloadButton("download_RF_Mod1", "Download Random Forest Model"),
                                               tags$hr()),
                                     
                                     nav_panel("Importance",
                                               uiOutput("RF_VIPPlot_Title"),
                                               tags$hr(),
                                               
                                               bslib::accordion(
                                                 id = "RF_VIP1",
                                                 open = TRUE,
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for All ML Model",
                                                   value = "",
                                                   plotOutput("RF_VIPPlot1", height = "700")
                                                 )
                                               ),
                                               
                                               bslib::accordion(
                                                 id = "RF_VIP2",
                                                 open = FALSE,
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for DHA ML Model",
                                                   value = "",
                                                   plotOutput("RF_VIPPlot2", height = "800px")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for n-3 DPA ML Model",
                                                   value = "",
                                                   plotOutput("RF_VIPPlot3", height = "800px")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for EPA ML Model",
                                                   value = "",
                                                   plotOutput("RF_VIPPlot4", height = "800px")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for AA ML Model",
                                                   value = "",
                                                   plotOutput("RF_VIPPlot5", height = "800px")
                                                 )
                                               ),
                                               
                                               selectInput("RF_Mod_Num", "Select which plot/model you want to download", ""),
                                               downloadButton("download_RF_VIPPlot", "Download Random Forest VIP Plot"),
                                               downloadButton("download_RF_Mod", "Download Random Forest Model"),
                                               tags$hr()))
                                   
                         ),
                         
                         nav_panel("Extreme Gradient Boosting",
                                   navset_card_underline(
                                     id = "EGB1",
                                     nav_panel("Parameters", uiOutput("XGB_Plot_Title"),
                                               bslib::accordion(
                                                 id = "EGB2",
                                                 open = TRUE,
                                                 bslib::accordion_panel(
                                                   title = "Model Table",
                                                   value = "",
                                                   dataTableOutput("XGB_Models"),
                                                   tags$hr(),
                                                   downloadButton("download_XGB_Table", "Download XGBoost Table")
                                                 )
                                               ),
                                               
                                               bslib::accordion(
                                                 id = "EGB3",
                                                 open = FALSE,
                                                 bslib::accordion_panel(
                                                   title = "Test Error Rate Plot for All ML Models",
                                                   plotlyOutput("XGB_Plot1")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "Test Error Rate Plot for DHA ML Model",
                                                   plotlyOutput("XGB_Plot2")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "Test Error Rate Plot for n-3 DPA  ML Model",
                                                   plotlyOutput("XGB_Plot3")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "Test Error Rate Plot for EPA ML Model",
                                                   plotlyOutput("XGB_Plot4")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "Test Error Rate Plot for AA ML Model",
                                                   plotlyOutput("XGB_Plot5")
                                                 )
                                               ),
                                               
                                               selectInput("XGB_Mod_Num", "Select which Model you want to download", ""),
                                               downloadButton("download_XGB_Plot", "Download XGBoost Plot"),
                                               downloadButton("download_XGB_Mod", "Download Extreme Gradien Boosting Model"),
                                               tags$hr()),
                                     nav_panel("Importance Plots", uiOutput("XGB_VIPPlot_Title"),
                                               
                                               bslib::accordion(
                                                 id = "XGB_VIP1",
                                                 open = TRUE,
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for ALL ML Models",
                                                   value = "",
                                                   plotOutput("XGB_VIPPlot1", height = "800px")
                                                 ),
                                               ),
                                               
                                               bslib::accordion(
                                                 id = "XGB_VIP2",
                                                 open = FALSE,
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for DHA ML Model",
                                                   value = "",
                                                   plotOutput("XGB_VIPPlot2", height = "800px")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for n-3 DPA ML Model",
                                                   value = "",
                                                   plotOutput("XGB_VIPPlot3", height = "800px")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for EPA Model",
                                                   value = "",
                                                   plotOutput("XGB_VIPPlot4", height = "800px")
                                                 ),
                                                 bslib::accordion_panel(
                                                   title = "VIP Plot for AA ML Model",
                                                   value = "",
                                                   plotOutput("XGB_VIPPlot5", height = "800px")
                                                 ),
                                               ),
                                               
                                               tags$hr(),
                                               selectInput("XGB_Mod_Num1", "Select which plot/Model you want to download", ""),
                                               downloadButton("download_XGB_VIPPlot", "Download XGBoost VIP Plot"),
                                               downloadButton("download_XGB_Mod1", "Download Extreme Gradien Boosting Model"),
                                               tags$hr()))
                         ),
                         
                         nav_panel("Support Vector Machine",
                                   uiOutput("SVM_Plot_Title"),
                                   tags$hr(),
                                   
                                   bslib::accordion(
                                     id = "SVM1",
                                     open = TRUE,
                                     bslib::accordion_panel(
                                       title = "Average Percent accuracy Plot of All ML Models",
                                       value = "",
                                       plotOutput("SVM_Plot1", height = "500px"),
                                     )
                                   ),
                                   bslib::accordion(
                                     id = "SVM2",
                                     open = FALSE,
                                     bslib::accordion_panel(
                                       title = "Average Percent accuracy Plot of DHA ML Model",
                                       value = "",
                                       plotOutput("SVM_Plot2", height = "500px"),
                                     ),
                                     bslib::accordion_panel(
                                       title = "Average Percent accuracy Plot of n-3 DPA ML Model",
                                       value = "",
                                       plotOutput("SVM_Plot3", height = "500px"),
                                     ),
                                     bslib::accordion_panel(
                                       title = "Average Percent accuracy Plot of EPA ML Model",
                                       value = "",
                                       plotOutput("SVM_Plot4", height = "500px"),
                                     ),
                                     bslib::accordion_panel(
                                       title = "Average Percent accuracy Plot of AA ML Model",
                                       value = "",
                                       plotOutput("SVM_Plot5", height = "500px"),
                                     ),
                                   ),
                                   
                                   
                                   selectInput("SVM_Mod_Num", "Select which Model you want to download", ""),
                                   downloadButton("download_SVM_Plot", "Download SVM Plot"),
                                   downloadButton("download_SVM_Mod", "Download SVM Model"),
                                   tags$hr()),
                         
                         nav_panel("Elastic Net Regression",
                                   uiOutput("LA_Plot_Title"),
                                   tags$hr(),
                                   bslib::accordion(
                                     id = "ENR1",
                                     open = TRUE,
                                     bslib::accordion_panel(
                                       title = "Average Accuracy Plot for ALL ML Models",
                                       value = "",
                                       plotOutput("LA_Plot1", height = "500px"),
                                     )
                                     ),
                                   bslib::accordion(
                                     id = "ENR2",
                                     open = FALSE,
                                     bslib::accordion_panel(
                                       title = "Average Accuracy Plot for DHA ML Model",
                                       value = "",
                                       plotOutput("LA_Plot2", height = "500px"),
                                     ),
                                     bslib::accordion_panel(
                                       title = "Average Accuracy Plot for n-3 DPA ML Model",
                                       value = "",
                                       plotOutput("LA_Plot3", height = "500px"),
                                     ),
                                     bslib::accordion_panel(
                                       title = "Average Accuracy Plot for EPA ML Model",
                                       value = "",
                                       plotOutput("LA_Plot4", height = "500px"),
                                     ),
                                     bslib::accordion_panel(
                                       title = "Average Accuracy Plot for AA ML Model",
                                       value = "",
                                       plotOutput("LA_Plot5", height = "500px"),
                                     )
                                   ),
                                   
                                   
                                   selectInput("LA_Mod_Num", "Select which Model you want to download", ""),
                                   downloadButton("download_LA_Plot", "Download Elastic Net Plot"),
                                   downloadButton("download_LA_Mod", "Download Elastic Net Model"),
                                   tags$hr()),
                         
                         nav_panel("Bayesian Linear Model",
                                   selectInput("BC_Mod_Num", "Select which Model you want to download", ""),
                                   downloadButton("download_BC_Mod", "Download Bayes GLM Model")),
                         
                         )
                       )
                     ),
                   
                   tabPanel(title ="Optimize ML Model",
                            id = "ML2",
                            titlePanel("Optimize Machine Learning Model"),
                            page_sidebar(
                              
                              sidebar = sidebar(
                                # Input: Select a file ----
                                width = 270,
                                
                                
                                fileInput("LM_ML_File2", "Lipid Mediator profiling file:",
                                          multiple = FALSE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv", ".tsv", ".txt")),
                                
                                # Input: Select separator ----
                                radioButtons("sep_LM_ML2", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                # Horizontal line ----
                                tags$hr(style = "border-top: 1px solid #000000;"),
                                
                                #option to build a random forest model 
                                checkboxInput("RF_Build", "Do you want to build a Random Forests Model?", FALSE),
                                
                                conditionalPanel(
                                  condition = "input.RF_Build == true",
                                  #select number of trees
                                  numericInput("RF_ntrees", "Select Number of Trees", 10000, min = 1, max = NA)),
                                
                                # Horizontal line ----
                                tags$hr(style = "border-top: 1px solid #000000;"),
                                
                                #option to build a XGBoost model 
                                checkboxInput("XGB_Build", "Do you want to build a Extreme Gradient Boosting Model?", FALSE),
                                conditionalPanel(
                                  condition = "input.XGB_Build == true",
                                  #select number of rounds 
                                  numericInput("XGB_nrounds", "Select Number of Rounds", 10000, min = 1, max = NA),
                                  
                                  #select model for All LM
                                  numericInput("XGB_ALL_LM_Model", "Select Model for All LM.", 1, min = 1, max = 5),
                                  
                                  #select model for DHA
                                  checkboxInput("XGB_Build_DHA", "Does Data contain DHA SPMs?", FALSE),
                                  conditionalPanel( 
                                    condition = "input.XGB_Build_DHA == true",
                                    numericInput("XGB_DHA_Model", "Select Model for DHA.", 1, min = 1, max = 5)),
                                  
                                  #Select model for n3-DPA
                                  checkboxInput("XGB_Build_n3DPA", "Does Data contain n3-DPA SPMs?", FALSE),
                                  conditionalPanel( 
                                    condition = "input.XGB_Build_n3DPA == true",
                                    numericInput("XGB_n3DPA_Model", "Select Model for n3-DPA.", 1, min = 1, max = 5)),
                                  
                                  #Select model for EPA
                                  checkboxInput("XGB_Build_EPA", "Does Data contain EPA SPMs?", FALSE),
                                  conditionalPanel( 
                                    condition = "input.XGB_Build_EPA == true",
                                    numericInput("XGB_EPA_Model", "Select Model for EPA.", 1, min = 1, max = 5)),
                                  
                                  #Select model for n3-DPA
                                  checkboxInput("XGB_Build_AA", "Does Data contain AA SPMs?", FALSE),
                                  conditionalPanel( 
                                    condition = "input.XGB_Build_AA == true",
                                    numericInput("XGB_AA_Model", "Select Model for AA.", 1, min = 1, max = 5))
                                ),
                                # Horizontal line ----
                                tags$hr(style = "border-top: 1px solid #000000;"),
                                
                                #option to build a support vector mechanism model 
                                checkboxInput("SVM_Build", "Do you want to build a Support Vector Machine Model?", FALSE),
                                conditionalPanel(
                                  condition = "input.SVM_Build == true",
                                  #select number of trees
                                  numericInput("SVM_Ensemble", "Select Number of Ensembles", 70, min = 1, max = NA),
                                  numericInput("SVM_BootNum", "Select Number of Boostraps", 70, min = 1, max = NA)),
                                # Horizontal line ----
                                tags$hr(style = "border-top: 1px solid #000000;"),
                                
                                #option to build a Lasso model 
                                checkboxInput("LA_Build", "Do you want to build a Elastic Net Regression Model?", FALSE),
                                conditionalPanel(
                                  condition = "input.LA_Build == true",
                                  #select number of trees
                                  numericInput("LA_BootNum", "Select Number of Bootstraps", 70, min = 1, max = NA)),
                                # Horizontal line ----
                                tags$hr(style = "border-top: 1px solid #000000;"),
                                
                                #option to build a Bayseian model 
                                checkboxInput("BC_Build", "Do you want to build a Bayesian Linear Model?", FALSE),
                                conditionalPanel(
                                  condition = "input.BC_Build == true",
                                  #select number of trees
                                  numericInput("BC_BootNum", "Select Number of Bootstraps", 70, min = 1, max = NA)),
                                # Horizontal line ----
                                tags$hr(style = "border-top: 1px solid #000000;"),
                                
                                # Action bottom to create and run the ML models ---
                                actionButton("Build_ML","Build Machine Learning Models")
                              
                            ),
                            
                            navset_bar(
                              nav_panel("Accuracy",
                                        
                                        bslib::accordion(
                                          id = "AC_Opt1",
                                          open = TRUE,
                                          bslib::accordion_panel(
                                            title = "Accuracy ML Plot",
                                            value = "",
                                            plotOutput("Build_Accuracy_ML"),
                                            uiOutput("error_ML_Opt"),
                                            downloadButton("downloadBuild_Acc_ML_Plot", "Download Accuracy Plot")
                                          )),
                                        
                                        bslib::accordion(
                                          id = "AC_Opt2",
                                          open = FALSE,
                                          bslib::accordion_panel(
                                            title = "Accuracy ML Table",
                                            value = "",
                                            dataTableOutput("Build_ML_Table"),
                                            downloadButton("downloadBuild_Acc_ML_Table", "Download Accuracy Table")
                                          ))),
                              
                              nav_panel("Random Forests", 
                                       navset_card_underline(id = "BML1",
                                                   nav_panel("Parameters",
                                                             bslib::accordion(
                                                               id = "RF_Par_Opt1",
                                                               open = TRUE,
                                                               bslib::accordion_panel(
                                                                 title = "Average Percent Accuracy Plot for All ML Models",
                                                                 value = "",
                                                                 plotOutput("Build_RF_Plot1", height = "500px")
                                                               )
                                                             ),
                                                             bslib::accordion(
                                                               id = "RF_Par_Opt2",
                                                               open = FALSE,
                                                               bslib::accordion_panel(
                                                                 title = "Average Percent Accuracy Plot for DHA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_Plot2", height = "500px")
                                                               ),
                                                               bslib::accordion_panel(
                                                                 title = "Average Percent Accuracy Plot for n-3 DPA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_Plot3", height = "500px")
                                                               ),
                                                               bslib::accordion_panel(
                                                                 title = "Average Percent Accuracy Plot for EPA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_Plot4", height = "500px")
                                                               ),
                                                               bslib::accordion_panel(
                                                                 title = "Average Percent Accuracy Plot for AA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_Plot5", height = "500px")
                                                               )
                                                             ),
                                                            downloadButton("downloadBuild_RF_Plot", "Download Random Forest Plot"),
                                                            ),
                                                   nav_panel("Importance",
                                                             bslib::accordion(
                                                               id = "RF_Imp_Opt1",
                                                               open = TRUE,
                                                               bslib::accordion_panel(
                                                                 title = "VIP Plot for All ML Models",
                                                                 value = "",
                                                                 plotOutput("Build_RF_VIPPlot1", height = "800px"),
                                                               )
                                                             ),
                                                             
                                                             bslib::accordion(
                                                               id = "RF_Imp_Opt2",
                                                               open = FALSE,
                                                               bslib::accordion_panel(
                                                                 title = "VIP Plot for DHA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_VIPPlot2", height = "800px"),
                                                               ),
                                                               bslib::accordion_panel(
                                                                 title = "VIP Plot for n-3 DPA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_VIPPlot3", height = "800px"),
                                                               ),
                                                               bslib::accordion_panel(
                                                                 title = "VIP Plot for EPA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_VIPPlot4", height = "800px"),
                                                               ),
                                                               bslib::accordion_panel(
                                                                 title = "VIP Plot for AA ML Model",
                                                                 value = "",
                                                                 plotOutput("Build_RF_VIPPlot5", height = "800px"),
                                                               )
                                                             ),
                                                            
                                                            numericInput("Build_RF_VIPPlot_Num", "Select which plot you want to download", 1, min = 1, max = 5),
                                                            downloadButton("downloadBuild_RF_VIPPlot", "Download Random Forest VIP Plot"),
                                                            
                                                            ))),
                              
                              nav_panel("Extreme Gradient Boosting", 
                                        navset_card_underline(id = "BML2",
                                                              nav_panel("Parameters",
                                                                        bslib::accordion(
                                                                          id = "XGB_Par_Opt1",
                                                                          open = TRUE,
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate Table for All ML Models",
                                                                            value = "",
                                                                            dataTableOutput("XGB_Build_Models"),
                                                                            downloadButton("downloadBuild_XGB_Table", "Download XGBoost Table"),
                                                                          )
                                                                        ),
                                                                        
                                                                        bslib::accordion(
                                                                          id = "XGB_Par_Opt2",
                                                                          open = FALSE,
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate for All ML Models",
                                                                            value = "",
                                                                            plotlyOutput("Build_XGB_Plot1", height = "500px")
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate for DHA ML Model",
                                                                            value = "",
                                                                            plotlyOutput("Build_XGB_Plot2", height = "500px")
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate for n-3 DPA ML Model",
                                                                            value = "",
                                                                            plotlyOutput("Build_XGB_Plot3", height = "500px")
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate for EPA ML Model",
                                                                            value = "",
                                                                            plotlyOutput("Build_XGB_Plot4", height = "500px")
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate for AA ML Model",
                                                                            value = "",
                                                                            plotlyOutput("Build_XGB_Plot5", height = "500px")
                                                                          )
                                                                        )
                                                              ),
                                                              nav_panel("Importance",
                                                                        bslib::accordion(
                                                                          id = "XGB_Vip_Opt1",
                                                                          open = TRUE,
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate Vip Plot for All Ml Models",
                                                                            value = "",
                                                                            plotOutput("Build_XGB_VIPPlot1", height = "800px")
                                                                          )
                                                                        ),
                                                                        bslib::accordion(
                                                                          id = "XGB_Vip_Opt2",
                                                                          open = FALSE,
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate Vip Plot for DHA ML Model",
                                                                            value = "",
                                                                            plotOutput("Build_XGB_VIPPlot2", height = "800px"),
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate Vip Plot for n-3 DPA ML Model",
                                                                            value = "",
                                                                            plotOutput("Build_XGB_VIPPlot3", height = "800px"),
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate Vip Plot for EPA ML Model",
                                                                            value = "",
                                                                            plotOutput("Build_XGB_VIPPlot4", height = "800px"),
                                                                          ),
                                                                          bslib::accordion_panel(
                                                                            title = "Test Error Rate Vip Plot for AA ML Model",
                                                                            value = "",
                                                                            plotOutput("Build_XGB_VIPPlot5", height = "800px"),
                                                                          ),
                                                                        ),
                                                                        
                                                                        tags$hr(),
                                                                        numericInput("Build_XGB_VIPPlot_Num", "Select which plot you want to download", 1, min = 1, max = 5),
                                                                        downloadButton("downloadBuild_XGB_VIPPlot", "Download XGBoost VIP Plot")
                                                                        
                                                              ))),
                              nav_panel("Support Vector Machine",
                                        bslib::accordion(
                                          id = "SVM_Opt1",
                                          open = TRUE,
                                          bslib::accordion_panel(
                                            title = "Average Percent Accuracy Plot for All ML Models",
                                            value = "",
                                            plotOutput("Build_SVM_Plot1", height = "500px"),
                                          )
                                        ),
                                        bslib::accordion(
                                          id = "SVM_Opt2",
                                          open = FALSE,
                                          bslib::accordion_panel(
                                            title = "Average Percent Accuracy Plot for DHA ML Model",
                                            value = "",
                                            plotOutput("Build_SVM_Plot2", height = "500px"),
                                          ),
                                          bslib::accordion_panel(
                                            title = "Average Percent Accuracy Plot for n-3 DPA ML Model",
                                            value = "",
                                            plotOutput("Build_SVM_Plot3", height = "500px"),
                                          ),
                                          bslib::accordion_panel(
                                            title = "Average Percent Accuracy Plot for EPA ML Model",
                                            value = "",
                                            plotOutput("Build_SVM_Plot4", height = "500px"),
                                          ),
                                          bslib::accordion_panel(
                                            title = "Average Percent Accuracy Plot for AA ML Model",
                                            value = "",
                                            plotOutput("Build_SVM_Plot5", height = "500px"),
                                          )
                                        ),
                                        downloadButton("downloadBuild_SVM_Plot", "Download SVM Plot")),
                              
                              nav_panel("Elastic Net Regression",
                                       bslib::accordion(
                                         id = "ENR_Opt1",
                                         open = TRUE,
                                         bslib::accordion_panel(
                                           title = "Average Accuracy Plot for All ML Models",
                                           value = "",
                                           plotOutput("Build_LA_Plot1", height = "500px")
                                         )
                                       ),
                                       bslib::accordion(
                                         id = "ENR_Opt2",
                                         open = FALSE,
                                         bslib::accordion_panel(
                                           title = "Average Accuracy Plot for DHA ML Models",
                                           value = "",
                                           plotOutput("Build_LA_Plot2", height = "500px")
                                         ),
                                         bslib::accordion_panel(
                                           title = "Average Accuracy Plot for n-3 DPA ML Models",
                                           value = "",
                                           plotOutput("Build_LA_Plot3", height = "500px")
                                         ),
                                         bslib::accordion_panel(
                                           title = "Average Accuracy Plot for EPA ML Models",
                                           value = "",
                                           plotOutput("Build_LA_Plot4", height = "500px")
                                         ),
                                         bslib::accordion_panel(
                                           title = "Average Accuracy Plot for AA ML Models",
                                           value = "",
                                           plotOutput("Build_LA_Plot5", height = "500px")
                                         )
                                       ),
                                       downloadButton("downloadBuild_LA_Plot", "Download Elastic Net Regression Plot")),
                              nav_panel("Bayesian Linear Model",
                                        selectInput("Build_BC_Mod_Num", "Select which model you want to download", ""),
                                        downloadButton("downloadBuild_BC_Mod", "Download Bayes GLM Model")),
                              ))),
                   
                   
                   tabPanel(title ="Build ML Model",
                            id = "ML2_5",
                            titlePanel("Build Machine Learning Model"),
                            page_sidebar(
                              sidebar = sidebar(
                              
                              width = 270,
                              
                              # Input: Select a file ----
                              
                              
                              
                              fileInput("LM_ML_File2_5", "Lipid Mediator profiling file:",
                                        multiple = FALSE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv", ".tsv", ".txt")),
                              
                              # Input: Select separator ----
                              radioButtons("sep_LM_ML2_5", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              #input -- where the sample group info column is
                              selectInput("Group_ML2_5", "Select Group Column Name", "Please upload file"),
                              
                              #selecting metabolites of interest 
                              selectInput("MetName_ML", "Select all Metabolites of intrest", "Please upload file", multiple = TRUE, selectize = TRUE),
                              
                              tags$hr(),
                              
                              #option to build a random forest model 
                              checkboxInput("RF_BuildMet", "Do you want to build a Random Forests Model?", FALSE),
                              
                              #option to build a XGBoost model 
                              checkboxInput("XGB_BuildMet", "Do you want to build a Extreme Gradient Boosting Model?", FALSE),
                              
                              #option to build a support vector mechanism model 
                              checkboxInput("SVM_BuildMet", "Do you want to build a Support Vector Machine Model?", FALSE),
                              
                              #option to build a Lasso model 
                              checkboxInput("LA_BuildMet", "Do you want to build a Elastic Net Regression Model?", FALSE),
                              
                              #option to build a Bayseian model 
                              checkboxInput("BC_BuildMet", "Do you want to build a Bayesian Linear Model?", FALSE),
                              
                              tags$hr(),
                              
                              # Action bottom to create and run the ML models ---
                              actionButton("Build_ML2_5","Build Machine Learning Models")
                            ),
                            navset_bar(
                              id = "Build_ML1",
                              nav_panel("Accuracy",
                                        bslib::accordion(
                                          id = "Acc_BLD1",
                                          open = TRUE,
                                          bslib::accordion_panel(
                                            title = "ML Model Accuracy Plot",
                                            value = "",
                                            plotOutput("BuildMet_AccPlot"),
                                            uiOutput("error_ML_Build"),
                                            downloadButton("downloadBuildMet_AccPlot", "Download Accuracy Plot")
                                          )
                                        ),
                                        bslib::accordion(
                                          id = "Acc_BLD2",
                                          open = FALSE,
                                          bslib::accordion_panel(
                                            title = "ML Model Accuracy Table",
                                            value = "",
                                            dataTableOutput("BuildMet_ML_Table"),
                                            downloadButton("downloadBuildMet_Acc_ML_Table", "Download Accuracy Table")
                                          )
                                        )),
                              nav_panel("Random Forests", 
                                       navset_card_underline(id = "RF_build1",
                                                   nav_panel("Parameters",
                                                            plotOutput("BuildMet_RF_Plot", height = "500px"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_RF_Plot", "Download Random Forest Plot"),
                                                            downloadButton("downloadBuildMet_RF_Mod", "Download Random Forest Model")),
                                                   nav_panel("Importance",
                                                            plotOutput("BuildMet_RF_VIPPlot", height = "500px"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_RF_VIPPlot", "Download Random Forest VIP Plot"),
                                                            downloadButton("downloadBuildMet_RF_Mod1", "Download Random Forest Model")))),
                              nav_panel("Extreme Gradient Boosting", 
                                       navset_card_underline(id = "EGB_Build1",
                                                   nav_panel("Parameters",
                                                            plotlyOutput("BuildMet_XGB_Plot"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_XGB_Mod", "Download Extreme Gradient Boosting Model")),
                                                   nav_panel("Importance",
                                                            plotOutput("BuildMet_XGB_VIPPlot", height = "500px"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_XGB_VIPPlot", "Download XGBoost VIP Plot"),
                                                            downloadButton("downloadBuildMet_XGB_Mod1", "Download Extreme Gradient Boosting Model"),
                                                   ))),
                              nav_panel("Support Vector Machine",
                                        plotOutput("BuildMet_SVM_Plot", height = "500px"),
                                        tags$hr(),
                                        downloadButton("downloadBuildMet_SVM_Plot", "Download SVM Plot"),
                                        downloadButton("downloadBuildMet_SVM_Mod", "Download SVM Model")),
                              nav_panel("Elastic Net Regression",
                                       plotOutput("BuildMet_LA_Plot", height = "500px"),
                                       tags$hr(),
                                       downloadButton("downloadBuildMet_LA_Plot", "Download Elastic Net Regression Plot"),
                                       downloadButton("downloadBuildMet_LA_Mod", "Download Elastic Net Model")),
                              nav_panel("Bayesian Linear Model",
                                        tags$hr(),
                                        downloadButton("downloadBuildMet_BC_Mod", "Download Bayes GLM Model"))
                                          
                              ))),
                   
                   tabPanel(title ="Run ML Model",
                            id = "ML3",
                            titlePanel("Run Machine Learning Model"),
                            page_sidebar(
                              sidebar = sidebar(
                                introjsUI(),
                                width = 270,
                              
                              # Input: Select a file ----
                              actionButton("Help_ML2", "Tutorial"),
                              
                              
                              div(id = "ML_in2",
                                  fileInput("LM_ML_File3", "Lipid Mediator profiling file:",
                                            multiple = FALSE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv", ".tsv", ".txt")),
                                  
                                  # actionButton("filechoose=",label = "Pick a file")
                                  
                                  # Input: Select separator ----
                                  radioButtons("sep_LM_ML3", "Separator",
                                               choices = c(Comma = ",",
                                                           Semicolon = ";",
                                                           Tab = "\t"),
                                               selected = ",")
                              ),
                              
                              downloadButton("Train_ML", "Example File"),
                              
                              div(id = "ML_file_up",
                                  #input -- where the sample group info column is
                                  selectInput("Group_ML3", "Select Group Column Name", "Please upload file"),
                                  
                                  
                                  #checkbox input
                                  checkboxInput("Choose_MetName3", "Do you want to filter the metabolites included?", FALSE),
                              ),
                              
                              
                              #conditional panel
                              conditionalPanel(
                                condition = "input.Choose_MetName3 == true",
                                #selecting metabolites of interest 
                                selectInput("MetName_ML3", "Select all Metabolites of intrest", "Please upload file", multiple = TRUE, selectize = TRUE)
                                
                              ),
                              
                              
                              # Horizontal line ----
                              tags$hr(style = "border-top: 1px solid #000000;"),
                              
                              #option to build a random forest model 
                              checkboxInput("RF_Run", "Do you want to build a Random Forests Model?", FALSE),
                              
                              conditionalPanel(
                                condition = "input.RF_Run == true",
                                fileInput("RF_Mod_File", "Random Forest Model:",
                                          multiple = FALSE,
                                          accept = ".rds")),
                              
                              # Horizontal line ----
                              tags$hr(style = "border-top: 1px solid #000000;"),
                              
                              #option to build a XGBoost model 
                              checkboxInput("XGB_Run", "Do you want to Run an Extreme Gradient Boosting Model?", FALSE),
                              conditionalPanel(
                                condition = "input.XGB_Run == true",
                                
                                fileInput("XGB_Mod_File", "XGBoost Model:",
                                          multiple = FALSE,
                                          accept = ".rds")),
                              # Horizontal line ----
                              tags$hr(style = "border-top: 1px solid #000000;"),
                              
                              #option to build a support vector mechanism model 
                              checkboxInput("SVM_Run", "Do you want to Run a Support Vector Machine Model?", FALSE),
                              conditionalPanel(
                                condition = "input.SVM_Run == true",
                                fileInput("SVM_Mod_File", "SVM Model:",
                                          multiple = FALSE,
                                          accept = ".rds")),
                              # Horizontal line ----
                              tags$hr(style = "border-top: 1px solid #000000;"),
                              
                              #option to build a Lasso model 
                              checkboxInput("LA_Run", "Do you want to Run a Elastic Net Regression Model?", FALSE),
                              conditionalPanel(
                                condition = "input.LA_Run == true",
                                
                                fileInput("LA_Mod_File", "Elastic Net Regression Model:",
                                          multiple = FALSE,
                                          accept = ".rds")),
                              
                              # Horizontal line ----
                              tags$hr(style = "border-top: 1px solid #000000;"),
                              
                              #option to build a Bayseian model 
                              checkboxInput("BC_Run", "Do you want to run a Bayesian Linear Model?", FALSE),
                              conditionalPanel(
                                condition = "input.BC_Run == true",
                                
                                #upload BC model 
                                fileInput("BC_Mod_File", "Bayes Regression Model:",
                                          multiple = FALSE,
                                          accept = ".rds")),
                              # Horizontal line ----
                              tags$hr(style = "border-top: 1px solid #000000;"),
                              
                              # Action bottom to create and run the ML models ---
                              actionButton("Run_ML","Run Machine Learning Models")
                            ),
                            navset_bar(
                              nav_panel("Accuracy",
                                        bslib::accordion(
                                          id = "Acc_RML1",
                                          open = TRUE,
                                          bslib::accordion_panel(
                                            title = "Accuracy plot of ML Model",
                                            value = "",
                                            plotOutput("ROC_Plot", height = "600px"),
                                            uiOutput("error_ML_Run"),
                                            tags$hr(),
                                            downloadButton("downloadRun_ROCPlot", "Download Accuracy Plot")
                                          )
                                        ),
                                        bslib::accordion(
                                          id = "Acc_RML2",
                                          open = FALSE,
                                          bslib::accordion_panel(
                                            title = "Accuracy table of ML Model",
                                            value = "",
                                            dataTableOutput("ROC_Table"),
                                            tags$hr(),
                                            downloadButton("downloadRun_ROCTable", "Download Accuracy Table")
                                          ))),
                              nav_panel("Random Forest ROC",
                                       plotOutput("Run_RF_ROCPlot", height = "500px"),
                                       tags$hr(),
                                       downloadButton("downloadRun_RF_ROCPlot", "Download Random Forest ROCPlot")),
                              nav_panel("Extreme Gradient Boosting ROC",
                                       plotOutput("Run_XGB_ROCPlot", height = "500px"),
                                       tags$hr(),
                                       downloadButton("downloadRun_XGB_ROCPlot", "Download XGboost ROC Plot")),
                              nav_panel("Support Vector Machine ROC",
                                       plotOutput("Run_SVM_ROCPlot", height = "500px"),
                                       tags$hr(),
                                       downloadButton("downloadRun_SVM_ROCPlot", "Download SVM ROC Plot")),
                              nav_panel("Elastic Net Regression ROC",
                                       plotOutput("Run_LA_ROCPlot", height = "500px"),
                                       tags$hr(),
                                       downloadButton("downloadRun_LA_ROCPlot", "Download Elastic Net ROC Plot")),
                              nav_panel("Bayesian Regression ROC",
                                       plotOutput("Run_BC_ROCPlot", height = "500px"),
                                       tags$hr(),
                                       downloadButton("downloadRun_BC_ROCPlot", "Download Bayes Classifier ROC Plot"))
                                          
                              )))
                   
                 ),
)

server <- shinyServer(function(input, output, session) {
  
  output$downTest <- downloadHandler(filename = function(){"Example_Data.csv"}, 
                                          content = function(fname){
                                            write.csv(Test_data, fname, row.names = FALSE)})
  output$TrainTest <- downloadHandler(filename = function(){"Example_Data.csv"}, 
                                     content = function(fname){
                                       write.csv(Train_Data, fname, row.names = FALSE)})
  output$Net_test <- downloadHandler(filename = function(){"Example_Data.csv"}, 
                                      content = function(fname){
                                        write.csv(Network_Data, fname, row.names = FALSE)})
  
  output$Train_ML <- downloadHandler(filename = function(){"Example_Data.csv"}, 
                                     content = function(fname){
                                       write.csv(Train_RML, fname, row.names = FALSE)})
    
      # '<img src=/Random_Forest_Plotn3DPA />' 
  observeEvent(input$Help_PCA, { 
    introjs(session, options = list(
      "showBullets" = FALSE,
      "showProgress" = TRUE,
      "showStepNumbers" = FALSE,
      "nextLabel" = "Next",
      "prevLabel" = "Prev",
      "skipLabel" = "Skip",
      steps = list(
        list(element = "#multivar_stats_nav",
             intro = "Welcome to the Mediator Analysis Nexus, or MediAnex for sort. Through this and various other tutorials
            you can learn to use and navigate this application. This is the primary navigation bar which will allow you to move 
             between different sections and perform different analysis. ",
             position = "bottom"),
        list(element = "#multivar_stats_nav",
             intro = "To move to a different section all you need to do is simply press the Multivariate Statistics 
             or Machine Learning Models panels and then select one of the relevant options available.",
             position = "bottom"),
        list(element = "#PCA",
             intro = "This is the input field where you have the options to upload a specific file. To do so, simply press the browse button and then upload a relevant 
             file from your files. Most sections of this application accept comma-separated, semicolon-separated or tab-separated files.",
             position = "right"),
        list(element = "#PCA",
             intro = "Please note: This application accepts specific format and content for certain files so please feel free to familiarise yourself with this format. 
             To do so, simply press on the download exmample file button and this will allow you to download a csv file so you can become familiar with the format accepted. 
             Please be aware that sometimes column names need to be the same as in the example file.",
             position = "right"),
        list(element = "#PCA_col",
             intro = "The next step once you upload your file is to select from the drop-down menu the relevant data on which analysis will be performed. For the first 
             input please select the column containing the names of your samples. If you do not wish to indicate the sample 
             names on your graph, you can do so by un-ticking the box named sample names 
             and press the update plot button located in the main panel.",
             position = "right"),
        list(element = "#PCA_col",
             intro = "The last input on the second drop down menu should be the column containing your group names. This will allow separation into 2 (or more) groups
             based on the component created in the 2D-scores plot.",
             position = "right"),
        list(element = "#Go_PCAPLSDA",
             intro = "Once you have selected all the available parameters, the next step is to simply select whether you wish to perform PCA or PLS-DA and press the button. 
             Once analysis is complete, output will appear on the main panel. ",
             position = "right"),
        list(element = "#Go_PCAPLSDA",
             intro = "From there you can click on the panels of the second navigation bar to open the plots and table of results developed.",
             position = "right")
      )
    ))
  })
  
  observeEvent(input$Help_Diff, { 
    introjs(session, options = list(
      "showBullets" = FALSE,
      "showProgress" = TRUE,
      "showStepNumbers" = FALSE,
      "nextLabel" = "Next",
      "prevLabel" = "Prev",
      "skipLabel" = "Skip",
      steps = list(
        list(element = "#Diff_tut",
             intro = "This is the differential analysis section where you can perform statistical analysis tests 
             on the data used as input on the PCA/PLS-DA analysis section. This section also utilises the same 
             file format as PCA/PLS-DA.",
             position = "right"),
        list(element = "#Diff_tut",
             intro = "However, if you missed the previous tutorial or would like to make yourself 
             familiar with the format and file content accepted please use the download example file button.",
             position = "right"),
        list(element = "#normality",
             intro = "Once you have uploaded the relevant file containing your samples, groups and different 
             concentrations of SPMs, the next step will be to select the normality test which should be utilised.",
             position = "right"),
        list(element = "#normality",
             intro = "The application automatically reads the uploaded file and updates the drop-down menus
             with the relevant groups on which statistical analysis is performed. However, if there are more than 2 groups available,
             you can specify which grouops to compare. Once the normality test has been selected, the analysis can commence 
             by pressing the calculate button.",
             position = "right")
      )
    ))
  })
  
  observeEvent(input$Help_Net, { 
    introjs(session, options = list(
      "showBullets" = FALSE,
      "showProgress" = TRUE,
      "showStepNumbers" = FALSE,
      "nextLabel" = "Next",
      "prevLabel" = "Prev",
      "skipLabel" = "Skip",
      steps = list(
        list(element = "#SPM",
             intro = "This is the Network pathway analysis section, which utilises the output of Differential analysis
             to develop dynamic and interactive visulisation networks of 4 fatty acids and their constituent SPMs. These 
             fatty acids include the Arachidonic Acid (AA), Eicosaptaenoic Acid (EPA), Docosahexaenoic Acid (DHA) 
             and Docosapentaenoic Acid (n-3 DPA). ",
             position = "right"),
        list(element = "#SPM",
             intro = "Once the input file is uploaded, it automatically reads the adjusted p-values or VIP scores available and 
             based on your selection, it will update the node colours, sizes and shapes of SPMs to reflect on regulation patterns 
             and if those are significant.",
             position = "right"),
        list(element = "#Net_select",
             intro = "In this section, you have to select from the drop-down menu the name of one of the fatty acids and by pressing
             the generate network button, the application will output the relevant network with the updated nodes to reflect on regulation 
             changes.",
             position = "right"),
        list(element = "#Color_options",
             intro = "In this section you can make changes to the color of the nodes. 
             Simply by selecting a node, then selecting the preferred color from the ones available in the color wheel and clicking the update button will allow you to change the color. 
             The same can be done below for the edge colors.",
             position = "right")
      )
    ))
  })
  
  observeEvent(input$Help_ML, { 
    introjs(session, options = list(
      "showBullets" = FALSE,
      "showProgress" = TRUE,
      "showStepNumbers" = FALSE,
      "nextLabel" = "Next",
      "prevLabel" = "Prev",
      "skipLabel" = "Skip",
      steps = list(
        list(element = "#ML_input",
             intro = "This is the Machine learning section where you can upload files, train models and then 
             download the ML models to use in the last ML panel, named Run ML. All the ML sections have the same premise
             so a tutorial has been provided only for this one.",
             position = "right"),
        list(element = "#ML_input",
             intro = "Similar to previous sections there is another example file button to help understand the format and
             content of the file which should be used as input. ",
             position = "right"),
        list(element = "#Select_ML",
             intro = "In this drop-down menu you have to specify which column of your file contains the group names. 
             The application will use the group names to train the ML model into differentiating between phenotypes 
             and the trained ML model can then be downloaded and used on datasets in the Run ML section.",
             position = "right"),
        list(element = "#Gen_ML",
             intro = "The next step is to select which machine learning method to use and then by pressing the create 
             machine learning models button, you can start building your ML model.",
             position = "right"),
        list(element = "#Gen_ML",
             intro = "Please note: Some ML methods may take more time to build than others ranging from a few seconds to a maximum
             of one hour and 40 minutes. However, this will depend on your system and componentns as well as how many models you choose
             to run simultaneously.",
             position = "right"),
        list(element = "#Gen_ML",
             intro = "Furthermore, this section tends to take more time as a very large number of iterations is used. In the other 
                  two sections (Optimise and build ML) you will have the chance to change some of the parameters and as such Models generally
                  take less time to run.",
             position = "right")
             
      )
    ))
  })
  observeEvent(input$Help_ML2, { 
    introjs(session, options = list(
      "showBullets" = FALSE,
      "showProgress" = TRUE,
      "showStepNumbers" = FALSE,
      "nextLabel" = "Next",
      "prevLabel" = "Prev",
      "skipLabel" = "Skip",
      steps = list(
        list(element = "#ML_in2",
             intro = "This is the final page of the ML models section. This part is the same as the other ML sections
             where you can upload your file and then if needed select specific metabolites of interest to train the ML models on.",
             position = "right"),
        list(element = "#ML_file_up",
             intro = "This part is the new element found only in this ML page. Once the user selects which ML method to use, a new input box 
             opens up. In this input field you are able to upload a file of a ML model trained on one of the previous ML pages and use it here 
             on actual data.",
             position = "right"),
        list(element = "#ML_file_up",
             intro = "With that you have reached the end of the application and the final tutorial. Each one of these sections can run as a work flow
             from beginning to end or you can start in each one for specific analysis. If you wish to complete any of the tutorials again,
             you can do do simply by pressing the tutorial button.",
             position = "right")
    )))
  })
  
  observeEvent(input$Test_intro, { introjs(session, options = list("showBullets"="false", 
                                                                   "showProgress"="true", 
                                                                   "showStepNumbers"="false",
                                                                   "nextLabel"="Next",
                                                                   "prevLabel"="Prev",
                                                                   "skipLabel"="Skip", 
                                                                   steps = data.frame(element = c("#sep_Diff"),
                                                                   intro = c("Testing that the box works"),
                                                                   position = c("right")))) })
  
  observe({
    tryCatch({
      
      
      req(input$PCAPLSDA_file)
      #form dataframe with file 
      df<- read.csv(input$PCAPLSDA_file$datapath,
                    header = input$header_PCA,
                    sep = input$sep_PCA)
      
      #transpose table if samples are at columns
      if(input$samp_PCA == "col_PCA"){
        df<-data.frame(t(df))
        df_samp<-rownames(df)
        df<-data.frame(cbind(df_samp, df))
        colnames(df)<- df[1,]
        df<- df[-1, ]}
      
      #get list of conditions as options 
      Group_Names<-unique(as.character(unname(unlist(colnames(df)))))
      updateSelectInput(session, "Group_PCA", choices = Group_Names, selected = Group_Names[[1]])
      updateSelectInput(session, "SampleCol_PCA", choices = Group_Names, selected = Group_Names[[1]])
      
      output$error_PC_PLS <- renderUI(NULL)
      },
      error = function(e) {
        output$error_PC_PLS <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          For more information, please look at the example file available through the 
          download example file button.", e$message
          )
        })
      }
    )
  })
  
  observeEvent(input$Go_PCAPLSDA|input$Go_PCAPLSDA2, {
    tryCatch({
      
      
      
      req(input$PCAPLSDA_file)
      #form dataframe with file 
      df<- read.csv(input$PCAPLSDA_file$datapath,
                    header = input$header_PCA,
                    sep = input$sep_PCA)
      
      #transpose table if samples are at columns
      if(input$samp_PCA == "col_PCA"){
        df<-data.frame(t(df))
        df_samp<-rownames(df)
        df<-data.frame(cbind(df_samp, df))
        colnames(df)<- df[1,]
        df<- df[-1, ]}
      
      Met_Rename<-c()
      
      #remove spaces from column names 
      colnames(df)<-gsub("\\.","", colnames(df)) 
      
      for (x in 1:ncol(df)){
        if (colnames(df)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(df)[x]]]
        }else{
          Met_Rename[x] <-colnames(df)[x]
        }
      }
      colnames(df)<-Met_Rename
      
      #script for PCA analysis
      if(input$test_PCAPLSDA == "PCA"){
        withProgress(message = 'Initiating PCA Analysis:', detail = "This may take a while.", value = 0, {  
          GroupName<- as.character(input$Group_PCA)
          SampleName<-as.character(input$SampleCol_PCA)
          
          #save group columns
          Groups_PCA<-df[,GroupName]
          Sample_PCA<-df[,SampleName]
          
          #remove group column from dataframe
          #df_PCA<-df[,-c(input$SampCol_PCA)]
          df_PCA<-df[ , !(names(df) %in% c(GroupName, SampleName))]
          
          #ensure dataframe is numeric
          df_PCA<-data.frame(sapply(df_PCA, function(x) as.numeric(as.character(x))))
          
          #zero handling data
          cols<- 1:(ncol(df_PCA))
          #replace zeros for each column to 1/5 the smallest value for each column
          df_PCA[cols] <- lapply(df_PCA[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
          
          #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
          df_PCA<-as.matrix(df_PCA)
          df_PCA[is.infinite(df_PCA)] <- NA
          min_index <- which.min(df_PCA)
          zero_replace <- (df_PCA[min_index]/5)
          df_PCA <- as.data.frame(df_PCA)
          df_PCA[is.na(df_PCA)] <- zero_replace
          
          #exclude columns with the same value for all rows
          df_PCA<-Filter(var, df_PCA)
          
          
          #PCA test with normalization 
          PCA<-prcomp(df_PCA, scale = TRUE)
          
          #----------PCA Score (PC1vsPC2)
          
          #get PCA Scores into dataframe and change it to numeric
          PCA_Scores<-data.frame(PCA$x)
          PCA_Scores<-data.frame(sapply(PCA_Scores, function(x) as.numeric(as.character(x))))
          PCA_Scores<-(cbind(Sample_PCA, Groups_PCA, PCA_Scores))
          
          #PC1 vs PC2 scores plot
          if(input$ViewSamples_PCA == TRUE){
            output$Xscoresplot<-renderPlotly({ggplotly(ggplot(PCA_Scores, aes(x=PC1, y=PC2, color = Groups_PCA)) + geom_point() + 
                                                         stat_ellipse(level = 0.95, geom = "polygon", aes(fill = Groups_PCA), alpha = 0.25) + theme_bw() + 
                                                         geom_text(label=Sample_PCA, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T))})
          }
          else{
            output$Xscoresplot<-renderPlotly({ggplotly(ggplot(PCA_Scores, aes(x=PC1, y=PC2, color = Groups_PCA)) + geom_point() + 
                                                         stat_ellipse(level = 0.95, geom = "polygon", aes(fill = Groups_PCA), alpha = 0.25) + theme_bw())})
          }
          
          #download PC scores 
          output$DownloadScore <- downloadHandler(filename = function(){"PCA_Score.csv"}, 
                                                  content = function(fname){
                                                    write.csv(PCA_Scores, fname, row.names = FALSE)})
          
          #----------PCA Loadings 
          
          #get loading scores, convert to numeric, get row names as column 
          loading<-data.frame(PCA$rotation)
          Variable<-rownames(loading)
          
          loading<-data.frame(sapply(loading, function(x) as.numeric(as.character(x))))
          loading<-cbind(Variable, loading)
          
          #change the names to be in the proper format 
          loading$Variable[grep("^X", loading$Variable)] <- paste(gsub("^X", "'", loading$Variable[grep("^X", loading$Variable)]), "'", sep = "")
          
          loading$Variable <- gsub("A4'", "A'[4]", loading$Variable) #converts names ending with 4 to subscript 5
          loading$Variable <- gsub("B4'", "B'[4]", loading$Variable) #B4 converted to subscript 4
          
          loading$Variable <- gsub("TXB2", "TXB[2]", loading$Variable)
          loading$Variable <- gsub("PGD2", "PGD[2]", loading$Variable)
          loading$Variable <- gsub("PGE2", "PGE[2]", loading$Variable)
          loading$Variable <- gsub("LXB4", "LXB[4]", loading$Variable)
          loading$Variable <- gsub("LXA4", "LXA[4]", loading$Variable)
          loading$Variable <- gsub("LTB4", "LTB[4]", loading$Variable)
          loading$Variable <- gsub("LTC4", "LTC[4]", loading$Variable)
          loading$Variable <- gsub("LTE4", "LTE[4]", loading$Variable)
          loading$Variable <- gsub("LTD4", "LTD[4]", loading$Variable)
          
          
          loading$Variable <- gsub("\\.", "-", loading$Variable) # Replace "." with "-" when required
          loading$Variable <- gsub("PGF2a", "PGF[2~a]", loading$Variable)
          loading$Variable <- gsub("n-3-dpa", "[n-3]~DPA", loading$Variable)
          
          lipids <- as.character(loading$Variable)
          
          loading<-loading[order(-loading$PC1),] #order dataframe based on PC1 score
          loading$Variable<-factor(loading$Variable, levels = loading$Variable)
          
          Loading_plot<-ggplot(loading, aes(x=reorder(Variable, PC1), y=PC1)) + geom_bar(stat = "identity") + 
            scale_x_discrete(labels = parse(text = lipids)) +
            coord_flip() + labs(x = "Lipid Mediator", y= "Loading Value") +
            theme(axis.title = element_text(size = 20),
                  axis.text.x  = element_text(size = 15, hjust = 0.5),
                  axis.text.y  = element_text(size = 15, hjust = 1),
                  panel.background = element_rect(fill = "white")) 
          
          
          
          #plot the loading plot
          output$Loadingplot<-renderPlot(Loading_plot)
          
          
          #download loading plot
          output$downloadLoading_Plot <- downloadHandler(filename = function(){paste("Loading_Plot",'.png',sep='')},
                                                         content = function(file){
                                                           ggsave(file,plot=Loading_plot, width = 12, height = 12)})
          
          #download Loading Data
          output$DownloadLoading <- downloadHandler(filename = function(){"PCA_LoadingScore.csv"}, 
                                                    content = function(fname){
                                                      write.csv(loading, fname, row.names = FALSE)})
          
          #----------Percentage of Variance 
          
          #percentage of variance
          Importance<-summary(PCA) #summary of PCA without normalization
          Importance<-t(Importance$importance)#save only the importance matrix
          Component<-rownames(Importance)
          Importance<-data.frame(cbind(Component, Importance[,2]))
          colnames(Importance)[2]<-"Percentage_of_Variance"
          Importance$Percentage_of_Variance<-(as.numeric(Importance$Percentage_of_Variance)*100)
          Importance$Component<-factor(Importance$Component, levels = Importance$Component)
          
          #plot the percentage variance 
          output$PerVarplot<- renderPlotly(ggplotly(ggplot(Importance, aes(Component, Percentage_of_Variance)) + 
                                                      geom_bar(stat = "identity", fill= "skyblue") + labs(y="Percentage of Variance") +
                                                      theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
          
          #download Percent Variance Data
          output$DownloadVar <- downloadHandler(filename = function(){"PCA_PercentVariance.csv"}, 
                                                content = function(fname){
                                                  write.csv(Importance, fname, row.names = FALSE)})
          
          #ensures that the VIP plot is 
          output$VIPplot<-NULL
          output$VIPtable<-NULL
          output$DownloadVIP <-NULL
          output$VIPtext <- renderText({ 
            "Please select PLS-DA Analysis to obtain VIP Scores "
          })
        }) 
      }
      #script for PLS-DA analysis
      else{
        withProgress(message = 'Initiating PLS-DA Analysis:', detail = "This may take a while.", value = 0, {
          if (nrow(df) < 7) {
            output$error_PLSDA <- renderUI({
              div(
                class = "alert alert-danger", role = "alert",
                strong("Error: "), "The uploaded file must contain 4 samples per group to ensure PLS-DA analysis runs."
              )
            })
            output$Xscoresplot <- renderPlotly(NULL)
          } else if (ncol(df) < 4) {
            output$error_PLSDA2 <- renderUI({
              div(
                class = "alert alert-danger", role = "alert",
                strong("Error: "), "The uploaded file must contain at least 3 SPMs to ensure 2 components are generated for PLS-DA analysis to run."
              )
            })
            output$Xscoresplot <- renderPlotly(NULL)
          } else if (!all(c("Group_names", "Sample_names") %in% colnames(df))) {
            output$error_PLSDA4 <- renderUI({
              div(
                class = "alert alert-danger", role = "alert",
                strong("Error: "), "There appears to be a problem with the format of the file.  
            Please ensure that the format of the uploaded file follows the example file. 
            Please ensure that the names of the columns in your file are named as in the example file. 
            For Example, your first column should contain sample with the column being named Sample_names and the second column should contain at least 2 groups and be named Group_names.
            For more information on the file format you can download an example through the Download Example Button of the sidebar section."
              )
            })
            output$Xscoresplot <- renderPlotly(NULL)
          } else if (ncol(df) >= 2 && (length(unique(df[[1]])) < 2 || length(unique(df[[2]])) < 2)) {
            output$error_PLSDA3 <- renderUI({
              div(
                class = "alert alert-danger", role = "alert",
                strong("Error: "), "The uploaded file must contain at least 2 groups to ensure PLS-DA analysis to runs."
              )
            })
            output$Xscoresplot <- renderPlotly(NULL)
          } else {
            
            output$error_PLSDA <- renderUI(NULL)
            output$error_PLSDA2 <- renderUI(NULL)
            output$error_PLSDA3 <- renderUI(NULL)
            output$error_PLSDA4 <- renderUI(NULL)
            
            GroupName<- as.character(input$Group_PCA)
            SampleName<-as.character(input$SampleCol_PCA)
            
            Groups_PLSDA <- df[[GroupName]]
            
            if (length(unique(Groups_PLSDA)) == length(Groups_PLSDA)) {
              output$error_PLSDA_group_select <- renderUI({
                div(
                  class = "alert alert-danger", role= "alert", 
                  strong("Error: "), "The column which you have selected, appears to have different values for each row. 
              Please ensure that the column you have selected contains group names and there are multiple samples per group in your file."
                )
              })
            } else {
              output$error_PLSDA_group_select <- renderUI(NULL)
              
            }
            
            #save group columns
            Groups_PLSDA = df[,GroupName]
            Sample_PLSDA = df[,SampleName]
            
            #remove group column from dataframe
            df_PLSDA<-df[ , !(names(df) %in% c(GroupName, SampleName))]
            Groups_PLSDA1<-as.factor(Groups_PLSDA)
            
            #ensure dataframe is numeric
            df_PLSDA<-data.frame(sapply(df_PLSDA, function(x) as.numeric(as.character(x))))
            
            #zero handling data
            cols<- 1:(ncol(df_PLSDA))
            #replace zeros for each column to 1/5 the smallest value for each column
            df_PLSDA[cols] <- lapply(df_PLSDA[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
            
            #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
            df_PLSDA<-as.matrix(df_PLSDA)
            df_PLSDA[is.infinite(df_PLSDA)] <- NA
            min_index <- which.min(df_PLSDA)
            zero_replace <- (df_PLSDA[min_index]/5)
            df_PLSDA <- as.data.frame(df_PLSDA)
            df_PLSDA[is.na(df_PLSDA)] <- zero_replace
            
            #exclude columns with the same value for all rows
            df_PLSDA<-Filter(var, df_PLSDA)
            
            #get the number of components using PCA test to get the length of the scores table = number of components 
            PCA<-prcomp(df_PLSDA, scale = TRUE)
            PCA_Scores<-data.frame(PCA$x)
            
            #run plsda() with data values and class
            PLSDA<- mixOmics::plsda(df_PLSDA, Groups_PLSDA, scale = TRUE, ncomp = (ncol(PCA_Scores)-1))
            
            
            #-------------PLS-DA X-scores
            
            #getting the scores and convert to numeric
            XScores<-data.frame(PLSDA$variates$X)
            XScores<-data.frame(sapply(XScores, function(x) as.numeric(as.character(x))))
            XScores<- cbind(Sample_PLSDA, XScores)
            
            #ploting comp1 vs comp 2 to show sample names 
            if(input$ViewSamples_PCA == TRUE){
              output$Xscoresplot<-renderPlotly({ggplotly(ggplot(XScores, aes(x=comp1, y=comp2, color = Groups_PLSDA)) + geom_point() + 
                                                           stat_ellipse(level = 0.95, geom = "polygon", aes(fill = Groups_PLSDA), alpha = 0.25) + 
                                                           theme_bw() + labs(x = "Component 1", y= "Component 2") + 
                                                           geom_text(label=Sample_PLSDA, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T))})
            }
            #plotting comp1 vs comp 2 wihtout sample names 
            else{
              output$Xscoresplot<-renderPlotly({ggplotly(ggplot(XScores, aes(x=comp1, y=comp2, color = Groups_PLSDA)) + geom_point() + 
                                                           stat_ellipse(level = 0.95, geom = "polygon", aes(fill = Groups_PLSDA), alpha = 0.25) + 
                                                           labs(x = "Component 1", y= "Component 2") + theme_bw())})
            }
            
            
            #download PC scores 
            output$DownloadScore <- downloadHandler(filename = function(){"PLSDA_Score.csv"}, 
                                                    content = function(fname){
                                                      write.csv(XScores, fname, row.names = FALSE)})
            
            ##get the p-value and the Q2 R2 plot 
            
            #get the R2 and Q2 values for the first component 
            PLSDA_R2Q2<-ropls::opls(df_PLSDA, Groups_PLSDA, predI = 2, fig.pdfC = "none", info.txtC = "none")
            Summ_R2Q2<-ropls::getSummaryDF(PLSDA_R2Q2)
            
            #output datatable
            R2Q2_Table<- data.frame(R2Y = Summ_R2Q2$`R2Y(cum)`, Q2Y = Summ_R2Q2$`Q2(cum)`, P.Value_R2Y = Summ_R2Q2$pR2Y , P.Value_Q2Y = Summ_R2Q2$pQ2)
            output$R2Q2_Table<-renderDataTable(R2Q2_Table, rownames = FALSE)
            output$R2Q2_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("R2, Q2, and P Values for PLSDA", align = "center") })
            
            #plot for R2Q2 values 
            R2Q2<-as.data.frame(t(R2Q2_Table[,1:2])) #get the R2 and Q2 values
            R2Q2<-cbind(rownames(R2Q2), data.frame(R2Q2, row.names=NULL)) #change rownames to the first column
            colnames(R2Q2)<-c("PLSDA_R2Q2", "Values") #change the column names 
            
            output$R2Q2_Plot<- renderPlotly(ggplotly(ggplot(R2Q2, aes(PLSDA_R2Q2, Values)) + 
                                                       geom_bar(stat = "identity", fill = c("gray", "black")) +
                                                       theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
            
            #-------------PLS-DA Loading Scores
            
            #obtan loading dataframes and convert to numeric 
            Xloading<-data.frame(PLSDA$loadings$X) #get x loadings as dataframe 
            Variable<-rownames(Xloading)
            Xloading<-as.data.frame(sapply(Xloading, function (x) as.numeric(as.character(x))))
            Xloading<-cbind(Variable, Xloading)
            
            #change the names to be in the proper format 
            Xloading$Variable[grep("^X", Xloading$Variable)] <- paste(gsub("^X", "'", Xloading$Variable[grep("^X", Xloading$Variable)]), "'", sep = "")
            Xloading$Variable <- gsub("A4'", "A'[4]", Xloading$Variable) #converts names ending with 4 to subscript 5
            Xloading$Variable <- gsub("B4'", "B'[4]", Xloading$Variable) #B4 converted to subscript 4
            
            Xloading$Variable <- gsub("TXB2", "TXB[2]", Xloading$Variable)
            Xloading$Variable <- gsub("PGD2", "PGD[2]", Xloading$Variable)
            Xloading$Variable <- gsub("PGE2", "PGE[2]", Xloading$Variable)
            Xloading$Variable <- gsub("LXB4", "LXB[4]", Xloading$Variable)
            Xloading$Variable <- gsub("LXA4", "LXA[4]", Xloading$Variable)
            Xloading$Variable <- gsub("LTB4", "LTB[4]", Xloading$Variable)
            Xloading$Variable <- gsub("LTC4", "LTC[4]", Xloading$Variable)
            Xloading$Variable <- gsub("LTE4", "LTE[4]", Xloading$Variable)
            Xloading$Variable <- gsub("LTD4", "LTD[4]", Xloading$Variable)
            
            Xloading$Variable <- gsub("\\.", "-", Xloading$Variable) # Replace "." with "-" when required
            Xloading$Variable <- gsub("PGF2a", "PGF[2~a]", Xloading$Variable)
            Xloading$Variable <- gsub("n-3-dpa", "[n-3]~DPA", Xloading$Variable) # Transform n-3 DPA as a subscript: 
            lipids <- as.character(Xloading$Variable)
            
            Xloading<-Xloading[order(-Xloading$comp1),] #order dataframe based on PC1 score
            Xloading$Variable<-factor(Xloading$Variable, levels = Xloading$Variable)
            
            XLoading_plot<-ggplot(Xloading, aes(x=reorder(Variable, comp1), y=comp1)) + geom_bar(stat = "identity") + 
              scale_x_discrete(labels = parse(text = lipids)) +
              coord_flip() + labs(x = "Lipid Mediator", y= "Loading Value") +
              theme(axis.title = element_text(size = 20),
                    axis.text.x  = element_text(size = 15, hjust = 0.5),
                    axis.text.y  = element_text(size = 15, hjust = 1),
                    panel.background = element_rect(fill = "white")) 
            
            #plot the loading plot
            output$Loadingplot<-renderPlot(XLoading_plot)
            
            
            #Download loading plot
            output$downloadLoading_Plot <- downloadHandler(filename = function(){paste("Loading_Plot",'.png',sep='')},
                                                           content = function(file){
                                                             ggsave(file,plot=XLoading_plot, width = 12, height = 12)})
            #download Loading Data
            output$DownloadLoading <- downloadHandler(filename = function(){"PLSDA_LoadingScore.csv"}, 
                                                      content = function(fname){
                                                        write.csv(Xloading, fname, row.names = FALSE)})
            
            #-------------PLS-DA Percentage Variance
            
            #Variance dataframes 
            XVariance<-data.frame(PLSDA$prop_expl_var$X)
            colnames(XVariance)[1]<-"Percentage_of_Var"
            Component<-rownames(XVariance)
            XVariance<-cbind(Component, XVariance)
            XVariance$Percentage_of_Var<-as.numeric(XVariance$Percentage_of_Var)*100
            XVariance$Component<- factor(XVariance$Component, levels = XVariance$Component)
            
            #percentage variance plots
            output$PerVarplot<- renderPlotly(ggplotly(ggplot(XVariance, aes(Component, Percentage_of_Var)) + 
                                                        geom_bar(stat = "identity", fill = "skyblue") + labs(y = "Percentage of Variance") +
                                                        theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
            #download Percent Variance Data
            output$DownloadVar <- downloadHandler(filename = function(){"PLSDA_PercentVariance.csv"}, 
                                                  content = function(fname){
                                                    write.csv(XVariance, fname, row.names = FALSE)})
            
            #-------------PLS-DA VIP Scores 
            
            #VIP scores plot
            VIP<-data.frame(vip(PLSDA)) #get VIP list
            Variable<-rownames(VIP) #get row names 
            VIP<-data.frame(sapply(VIP, function(x) as.numeric(as.character(x))))
            VIP<-data.frame(cbind(Variable, VIP)) #add variables as column
            
            
            Fatty_Acid <- 1:length(VIP$Variable)
            cbind(VIP, Fatty_Acid)
            
            #fatty acid points 
            for (i in 1:length(VIP$Variable)){
              VIP$Fatty_Acid[i]<-SPMs_FA[VIP$Variable[i]]
            }
            
            VIP$Fatty_Acid<-as.character(VIP$Fatty_Acid)
            
            #VIP["lm_names"]<-VIP$Variable
            #format table to be downloaded and used for pathway analysis 
            VIP_Table <- VIP
            VIP_Table$Variable<- gsub("\\.", "-", VIP_Table$Variable) # Replace "." with "-" when required 
            VIP_Table$Variable[grep("^X", VIP_Table$Variable)] <- paste(gsub("^X", "", VIP_Table$Variable[grep("^X", VIP_Table$Variable)]), "", sep = "") #remove x at the beginning
            
            VIP<-VIP[VIP["comp1"]>1,]
            
            #Change the naming format to be for a scientific paper
            VIP$Variable[grep("^X", VIP$Variable)] <- paste(gsub("^X", "'", VIP$Variable[grep("^X", VIP$Variable)]), "'", sep = "")
            
            VIP$Variable <- gsub("A4'", "A'[4]", VIP$Variable) #converts names ending with 4 to subscript 5
            VIP$Variable <- gsub("B4'", "B'[4]", VIP$Variable) #B4 converted to subscript 4
            
            VIP$Variable <- gsub("TXB2", "TXB[2]", VIP$Variable)
            VIP$Variable <- gsub("PGD2", "PGD[2]", VIP$Variable)
            VIP$Variable <- gsub("PGE2", "PGE[2]", VIP$Variable)
            VIP$Variable <- gsub("LXB4", "LXB[4]", VIP$Variable)
            VIP$Variable <- gsub("LXA4", "LXA[4]", VIP$Variable)
            VIP$Variable <- gsub("LTB4", "LTB[4]", VIP$Variable)
            VIP$Variable <- gsub("LTC4", "LTC[4]", VIP$Variable)
            VIP$Variable <- gsub("LTE4", "LTE[4]", VIP$Variable)
            VIP$Variable <- gsub("LTD4", "LTD[4]", VIP$Variable)
            
            VIP$Variable <- gsub("\\.", "-", VIP$Variable) # Replace "." with "-" when required
            VIP$Variable <- gsub("PGF2a", "PGF[2~a]", VIP$Variable)
            VIP$Variable <- gsub("n-3-dpa", "[n-3]~DPA", VIP$Variable)
            
            
            
            VIP<-VIP[order(VIP$comp1),] #order dataframe based on -comp1 score
            VIP$Variable <- factor(VIP$Variable, levels = VIP$Variable[order(VIP$comp1)]) #keeps variable column as factor
            #VIP$Fatty_Acid2 <- factor(VIP$Fatty_Acid, levels = c(1:length(VIP$Fatty_Acid)))
            #keeps variable column as factor
            
            lipids <- as.character(VIP$Variable)
            
            #add colors to SPM names
            lm_classes <- as.factor(VIP$Fatty_Acid)
            names(lm_classes) <- VIP$Variable
            
            dha_index<-which((lm_classes=="DHA")== TRUE)
            n_three_index<-which((lm_classes=="n3DPA")== TRUE)
            epa_index<-which((lm_classes=="EPA")== TRUE)
            aa_index<-which((lm_classes=="AA")== TRUE)
            
            lm_colors <- NULL
            lm_colors[dha_index] <- "blue"
            lm_colors[n_three_index] <- "brown"
            lm_colors[epa_index] <- "darkgoldenrod1"
            lm_colors[aa_index] <- "darkslategray"
            
            
            #ggplot2
            VIP_plot <- ggplot(data = VIP, aes(x = Variable, y = comp1, color = Fatty_Acid)) + geom_point(size =3) +
              scale_y_continuous(name = "VIP Score") +
              labs(x = "Lipid Mediators") +
              scale_x_discrete(labels = parse(text = lipids)) + 
              scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue", "n3DPA" = "brown", 
                                                                                   "EPA"="darkgoldenrod1",  "AA" ="darkslategray")) +
              coord_flip() + 
              theme(axis.title = element_text(size = 25),
                    axis.text.x  = element_text(size = 15, hjust = 0.5),
                    axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                    legend.position = "top",
                    aspect.ratio = 2/1,
                    legend.title = element_text(size = 20),
                    legend.text  = element_text(size = 15),
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
            
            
            
            #outputs 
            output$VIPtext<-NULL
            output$VIPplot<-renderPlot(VIP_plot)
            
            output$VIPtable<-renderDataTable(VIP, options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE), 
                                             selection = "single", rownames = FALSE)
            
            output$downloadPLSDA_VIP_Plot <- downloadHandler(filename = function(){paste("PLSDA_VIP_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot=VIP_plot, width = 12, height = 8)})
            output$DownloadVIP <- downloadHandler(filename = function(){"PLSDA_VIPScore.csv"}, 
                                                  content = function(fname){
                                                    write.csv(VIP_Table, fname, row.names = FALSE)})
            
            #Pathway Analysis PLSDA VIP score
            VIP_Table_Comp1<-VIP_Table[,1:2]
            #open csv with the pathways for each metabolite
            Pathway <-read.csv(file = "coefficients/Pathways2.csv", header = TRUE, sep = ",")
            #set it up like a dictionary 
            MetPath<-setNames(Pathway$Pathway, Pathway$Metabolite)
            VIP_Table_Comp1["Pathways"]<-"N/A"#column for pathways with N/A as default
            
            #loop to put the pathway for the metabolite 
            for (x in 1:length(VIP_Table_Comp1$Variable)){
              if (VIP_Table_Comp1$Variable[x] %in% Pathway$Metabolite){
                VIP_Table_Comp1$Pathways[x] <- MetPath[[VIP_Table_Comp1$Variable[x]]]
              }else{
                VIP_Table_Comp1$Pathways[x]<-VIP_Table_Comp1$Pathways[x]
              }
            }
            
            #Keep rows that have a VIP Score greater than 1 and pathway != N/A
            VIP_Table_Comp1<-VIP_Table_Comp1[VIP_Table_Comp1["comp1"]>1,]
            VIP_Table_Comp1<-VIP_Table_Comp1[VIP_Table_Comp1["Pathways"] != "N/A",]
            
            #get list of unique pathways
            VIP_PathList<-unique(VIP_Table_Comp1$Pathways)
            #empty table
            VIPPath_Table<-data.frame(Pathway = "del", Number_SPM = 1)
            #loop to calculate the number of SPMs for each pathway
            for (x in VIP_PathList){
              x_VIP <- VIP_Table_Comp1[VIP_Table_Comp1["Pathways"] == x, ]
              Table<-data.frame(Pathway = x, Number_SPM = length(x_VIP$Pathways))
              VIPPath_Table<-rbind(VIPPath_Table, Table)
            }
            #remove plot 
            VIPPath_Table<-VIPPath_Table[-1,]
            VIPPath_Table<-VIPPath_Table[order(-VIPPath_Table$Number_SPM),] #reorder based on Number SPM
            
            #get columns in correct data format
            VIPPath_Table$Number_SPM<-as.numeric(VIPPath_Table$Number_SPM)
            VIPPath_Table$Pathway<-as.character(VIPPath_Table$Pathway)
            #save as factor to keep the order
            VIPPath_Table$Pathway<-factor(VIPPath_Table$Pathway, levels = VIPPath_Table$Pathway)
            
            #plot title and plot output
            output$VIP_Pathway_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("VIP Pathways Plot", align = "center") })
            output$VIP_PathwayPlot<- renderPlotly(ggplotly(ggplot(VIPPath_Table, aes(x=reorder(Pathway, -Number_SPM), Number_SPM)) + 
                                                             geom_bar(stat = "identity", fill = c("gray")) +
                                                             labs(x = "Pathway", y = "Number of SPMs")+
                                                             theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
            
          }
          
          
          
        }) }
      
      output$error_PC_PLS <- renderUI(NULL)
    },
    error = function(e) {
      output$error_PC_PLS <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          For more information, please look at the example file available through the 
          download example file button.", e$message
        )
      })
    }
    )
  })
  
  observe({
    tryCatch({
      inFile <- input$lm_profiles
      
      print(inFile)
      
      if(is.null(inFile))
        
        return(NULL)
      
      file_data <- read.csv(inFile$datapath, header = input$header_Diff, sep = input$sep_Diff)
      
      if (!all(c("Group_names", "Sample_names") %in% colnames(file_data))) {
        output$error_Diff_format <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the format of the file.  
            Please ensure that the format of the uploaded file follows the example file. 
            Please ensure that the names of the columns in your file are named as in the example file. 
            For Example, your first column should contain sample with the column being named Sample_names and the second column should contain at least 2 groups and be named Group_names.
            For more information on the file format you can download an example through the Download Example Button of the sidebar section."
          )
        })
      } 
      
      else if (anyDuplicated(file_data[, 1])) {
        output$error_Diff_Snames <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the content of the file.  
            Please ensure that the first column is named Sample_names and that there are no duplicate values in the rows of that column."
          )
        })
        return()
      } 
      
      else {
        
        output$error_Diff_format <- renderUI(NULL)
        lm_profiles = read.csv(inFile$datapath,
                               header=input$header_Diff,
                               sep=input$sep_Diff,
                               row.names = 1)
        lm_profiles<-data.frame(t(lm_profiles))}
      
      # This option allows to chose the groups to make the comparison:
      
      options <- unique(as.character(unname(unlist(lm_profiles[1, ]))))
      
      if (length(options) < 2) {
        output$error_Diff_groups <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the content of the file.  
            Please ensure that there are at least two groups available for differential analysis to take place."
          )
        })
      } else {
        output$error_Diff_format <- renderUI(NULL)
        output$error_Diff_groups <- renderUI(NULL)
        updateSelectInput(session, "group_A", choices = options, selected = options[[2]])
        updateSelectInput(session, "group_B", choices = options)
      }
      error_Diff <- renderUI(NULL)
    },
    error = function(e) {
      output$error_Diff <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          For more information, please look at the example file available through the 
          download example file button.", e$message
        )
      })
    }
    )
  })
  
  observeEvent(input$Go_Diff|input$Go_DiffPlot, {
    tryCatch({
      
      
      inFile <- input$lm_profiles
      
      print(inFile)
      
      if(is.null(inFile))
        
        return(NULL)
      
      file_data <- read.csv(inFile$datapath, header = input$header_Diff, sep = input$sep_Diff)
      
      # In the event that the file does not contain the specified columns by the name required in the application, 
      # the following error meassage will take effect.
      
      if (!all(c("Group_names", "Sample_names") %in% colnames(file_data))) {
        output$error_Diff_format <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the format of the file.  
            Please ensure that the format of the uploaded file follows the example file. 
            Please ensure that the names of the columns in your file are named as in the example file. 
            For Example, your first column should contain sample with the column being named Sample_names and the second column should contain at least 2 groups and be named Group_names.
            For more information on the file format you can download an example through the Download Example Button of the sidebar section."
          )
        })
        # The below output have been change to return NULL. This is done to prevent confusion
        # If the user finishes their analysis and uploads a new file which does not meet the requirements then 
        # the previous analysis is removed and the error message pops up.
        
        output$DiffTable <- renderDataTable(NULL)
        output$DiffPlot <- renderPlot(NULL)
        output$Diff_UpRegPath_Plot <- renderPlotly(NULL)
        output$Diff_DownRegPath_Plot <- renderPlotly(NULL)
        output$Diff_PValPath_Plot <- renderPlotly(NULL)
        
      }
      
      else if (anyDuplicated(file_data[, 1])) {
        output$error_Diff_Snames <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the content of the file.  
            Please ensure that the first column is named Sample_names and that there are no duplicate values in the rows of that column."
          )
        })
        return()
      } 
      
      else if (length(unique(file_data$Group_names)) < 2) {
        output$error_Diff_groups <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the content of the file.  
            Please ensure that there are at least two groups available for differential analysis to take place."
          )
        })
        
        output$DiffTable <- renderDataTable(NULL)
        output$DiffPlot <- renderPlot(NULL)
        output$Diff_UpRegPath_Plot <- renderPlotly(NULL)
        output$Diff_DownRegPath_Plot <- renderPlotly(NULL)
        output$Diff_PValPath_Plot <- renderPlotly(NULL)
        
      } 
      
      else {
        
        data_to_check <- file_data[2:nrow(file_data), 3:ncol(file_data)]
        
        # Define a function to check if a value can be converted to numeric
        is_numeric_value <- function(x) {
          # Try converting to numeric, return FALSE if NA or any conversion issue
          !is.na(suppressWarnings(as.numeric(x)))
        }
        
        # Initialize a list to store non-numeric details
        non_numeric_details <- list()
        
        # Check each column and row for non-numeric values
        for (col_index in 1:ncol(data_to_check)) {
          column <- data_to_check[, col_index]
          for (row_index in 1:length(column)) {
            value <- column[row_index]
            if (!is_numeric_value(value)) {
              non_numeric_details <- append(non_numeric_details, list(
                list(
                  sample = file_data[row_index + 1, 1], # Adjust for the row offset
                  group = file_data[row_index + 1, 2],
                  mediator = colnames(file_data)[col_index + 2], # Adjust for the column offset
                  row = row_index + 1, 
                  value = value
                )
              ))
            }
          }
        }
        
        # Check if there are any non-numeric values
        if (nrow(file_data) > 0 && length(non_numeric_details) > 0) {
          # Construct the error message
          error_message <- "There appears to be a problem with the content of the file. Some cells in the contain non-numeric values:\n\n"
          for (detail in non_numeric_details) {
            error_message <- paste(
              error_message,
              sprintf(
                "Sample: %s, Group: %s, Mediator: %s, Row: %d, Value: %s\n",
                detail$sample, detail$group, detail$mediator, detail$row, detail$value
              )
            )
          }
          
          output$error_Diff_numeric <- renderUI({
            div(
              class = "alert alert-danger", role = "alert",
              strong("Error: "), error_message
            )
          })
          output$DiffTable <- renderDataTable(NULL)
          output$DiffPlot <- renderPlot(NULL)
          output$Diff_UpRegPath_Plot <- renderPlotly(NULL)
          output$Diff_DownRegPath_Plot <- renderPlotly(NULL)
          output$Diff_PValPath_Plot <- renderPlotly(NULL)
        } 
        
        else { 
          output$error_Diff_numeric <- renderUI(NULL)
          
          lm_profiles = read.csv(inFile$datapath,
                                 header=input$header_Diff,
                                 sep=input$sep_Diff,
                                 row.names = 1)
          lm_profiles<-data.frame(t(lm_profiles))
          
          # Get the classification from the main file:
          
          classification <- data.frame(samples = colnames(lm_profiles),
                                       group = as.character(unname(unlist(lm_profiles[1, ]))))
          
          # Transform the lipid mediator concentrations to numeric values so we can do the calculations: 
          
          lipid_m_profile <- as.data.frame(sapply(lm_profiles[-1, ], function (x) as.numeric(as.character(x))))
          rownames(lipid_m_profile) <- rownames(lm_profiles[-1, ])
          
          #transpose dataframe
          lm_profiles_t <- data.frame(t(lipid_m_profile))
          
          # Replace all the zero values by NA (This is for the sake of replacing zero/missing values)
          lm_profiles_t[lm_profiles_t == 0] <- NA
          
          
          #zero handling data
          cols<- 1:(ncol(lm_profiles_t))
          #replace zeros for each column to 1/5 the smallest value for each column
          lm_profiles_t[cols] <- lapply(lm_profiles_t[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
          
          #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
          lm_profiles_t<-as.matrix(lm_profiles_t)
          lm_profiles_t[is.infinite(lm_profiles_t)] <- NA
          min_index <- which.min(lm_profiles_t)
          zero_replace <- (lm_profiles_t[min_index]/5)
          lm_profiles_t <- as.data.frame(lm_profiles_t)
          lm_profiles_t[is.na(lm_profiles_t)] <- zero_replace
          
          
          # MULTINORMALITY CHECK:
          
          # For normality it takes a transpose dataframe as the one we made before.
          lm_profiles_mvn <- as.data.frame(t(lipid_m_profile))
          
          # Since normality does not work with variables with zero values in all its samples, we delete this variables. 
          lm_profiles_mvn <- lm_profiles_mvn[, colSums(lm_profiles_mvn != 0) > 0]
          
          # Also, the function doesn't work well when variables has the same values (duplicates), so we delete them. 
          lm_profiles_mvn <- lm_profiles_mvn[!duplicated(as.list(lm_profiles_mvn))]
          
          # This function calculates multivariate normality:
          testing_mult_norm <- mvn(data = lm_profiles_mvn, mvnTest = input$mvn, covariance = TRUE, tol = 1e-100)
          
          if (input$mvn == "mardia") {mnv_result <- as.character(testing_mult_norm$multivariateNormality$Result[3])}
          else {mnv_result <- as.character(testing_mult_norm$multivariateNormality$MVN)}
          
          
          # Comparison:
          
          # Please especify the names of the groups that you want to compare as it is in the classification file. Here,
          # the comparison is going to be made as "Group A vs Group B". Meaning, that group B, generally, is going to be
          # your control group. 
          
          group_a <- input$group_A
          group_b <- input$group_B  
          
          
          # Specify classification table:
          
          classification <- classification[(classification$group == group_a) | (classification$group == group_b), ]
          
          lm_profiles_t <- lm_profiles_t[rownames(lm_profiles_t) %in% as.character(classification$samples), ]
          
          # Check if both of your group follows a normal distribution:
          
          lm_profiles_group_a <- lipid_m_profile[ ,colnames(lipid_m_profile) %in% 
                                                    as.character(classification$samples[classification$group == group_a])]
          
          lm_profiles_group_b <- lipid_m_profile[ ,colnames(lipid_m_profile) %in% 
                                                    as.character(classification$samples[classification$group == group_b])]
          
          toptable <- data.frame(name = "delete",
                                 FC = 1,
                                 logFC = 1,
                                 P.Value = 1,
                                 adj.P.Val = 1,
                                 test = "delete")
          
          
          
          for(i in 1:nrow(lipid_m_profile)) {
            
            if (mnv_result == "NO")  {
              
              # Calculates Mann_Whitney test:
              
              mw_test <- wilcox.test(lm_profiles_t[, i]~classification$group)
              
              # Calculates the Fold Change (A vs B):
              
              a <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_a]),i])
              b <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_b]),i])
              
              # Log FC:
              
              logFC <- log2(a) - log2(b)
              FC <- abs(a/b)
              
              # Put results together in a data frame:
              
              results_table <- data.frame(name = colnames(lm_profiles_t[i]),
                                          FC = format(FC, digits = 4),
                                          logFC = format(logFC, digits = 4),
                                          P.Value = format(mw_test$p.value, digits = 3),
                                          adj.P.Val = format(1, digits = 3),
                                          test = "Mann Whitney test")
              
              toptable <- rbind(toptable, results_table)
              
            }
            
            if (mnv_result == "YES") {
              
              t_test <- t.test(lm_profiles_t[, i]~classification$group)
              
              # Calculates the Fold Change (A vs B):
              
              a <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_a]),i])
              b <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_b]),i])
              
              # Log FC:
              
              logFC <- log2(a) - log2(b)
              FC <- abs(a/b)
              
              # Put results together in a data frame:
              
              results_table <- data.frame(name = colnames(lm_profiles_t[i]),
                                          FC = format(FC, digits = 4),
                                          logFC = format(logFC, digits = 4),
                                          P.Value = format(t_test$p.value, digits = 3),
                                          adj.P.Val = format(1, digits = 3),
                                          test = "t test")
              
              toptable <- rbind(toptable, results_table)
              
            }
            
          }
          
          toptable <- toptable[-1, ]
          toptable$adj.P.Val <- format(p.adjust(toptable$P.Value, "BH"), digits = 3)
          toptable$P.Value[toptable$P.Value == "NaN"] <- "1.000"
          toptable$adj.P.Val[toptable$adj.P.Val == "   NaN"] <- "1.000"
          toptable <- toptable[order(toptable$adj.P.Val), ]
          
          # #change the name to be in correct format 
          toptable$name[grep("^X", toptable$name)] <- paste(gsub("^X", "", toptable$name[grep("^X", toptable$name)]), "", sep = "")
          
          # toptable$name <- gsub("4$", "[4]", toptable$name) #converts names ending with 4 to subscript 5
          # toptable$name <- gsub("B4'", "B'[4]", toptable$name) #B4 converted to subscript 4
          toptable$name <- gsub("\\.", "-", toptable$name) # Replace "." with "-" when required 
          # toptable$name <- gsub("n.3.DPA", "[n-3~DPA]", toptable$name) # Transform n-3 DPA as a subscript: 
          # lipids <- as.character(toptable$name)
          
          #Volcano Plot 
          toptable$logFC<-as.numeric(toptable$logFC)
          toptable$adj.P.Val<-as.numeric(toptable$adj.P.Val)
          p_cutoff <- input$Sig_DiffPlot
          fc_cutoff <- input$FC_DiffPlot
          topN <- input$TopN_DiffPlot
          # VolPlot<- toptable %>% 
          #   mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
          #   mutate(Rank = 1:n(), Label = ifelse(Rank < topN, name,"")) %>% 
          #   ggplot(aes(x = logFC, y = -log(adj.P.Val), col=Significant, label=Label)) + geom_point(size = 3)
          
          VolPlot<- toptable %>% 
            mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
            mutate(Rank = 1:n(), Label = ifelse(Rank < topN, name,"")) %>% 
            ggplot(aes(x = logFC, y = -log(adj.P.Val), col=Significant, label=Label)) + geom_point(size = 3) 
          
          VolPlot<-VolPlot + geom_text_repel(col="black", size = 5) +labs(x=" LogFC", y="Adjusted p-value" ) + 
            theme_minimal(base_size = 16) +
            scale_color_manual(values = c('black', "red"))
          
          # OUTPUTS:
          
          output$DiffTable <- renderDataTable(toptable,selection = "single", rownames = FALSE)
          output$DiffPlot<-renderPlot(VolPlot)
          
          output$DownloadDiffTable <- downloadHandler(filename = function(){"Differential_Analysis.csv"}, 
                                                      content = function(fname){
                                                        write.csv(toptable, fname)})
          
          output$downloadDiffPlot <- downloadHandler(filename = function(){paste("VolcanoPlot",'.png',sep='')},
                                                     content = function(file){
                                                       ggsave(file,plot=VolPlot, width = 12, height = 8)})
          
          
          ##calculating the pathway analysis information-----
          
          #open csv with the pathways for each metabolite
          Pathway <-read.csv(file = "coefficients/Pathways2.csv", header = TRUE, sep = ",")
          #set it up like a dictionary 
          MetPath<-setNames(Pathway$Pathway, Pathway$Metabolite)
          toptable_P<-toptable
          toptable_P["Pathways"]<-"N/A"#column for pathways with N/A as default
          
          #loop to put the pathway for the metabolite 
          for (x in 1:length(toptable_P$name)){
            if (toptable_P$name[x] %in% Pathway$Metabolite){
              toptable_P$Pathways[x] <- MetPath[[toptable_P$name[x]]]
            }else{
              toptable_P$Pathways[x]<-toptable_P$Pathways[x]
            }
          }
          #remove all the rows missing the metabolite names
          toptable_Path<-toptable_P[toptable_P["Pathways"] != "N/A",]
          
          ###UPREGULATED PATHWAYS----
          toptable_UpReg<-toptable_Path[toptable_Path["logFC"] >0,]
          toptable_DownReg<-toptable_Path[toptable_Path["logFC"]<0,]
          
          UpReg_PathList<-unique(toptable_UpReg$Pathways)
          
          UpReg_Table<-data.frame(Pathway = "del", meanlogFC = 1)
          
          for (x in UpReg_PathList){
            x_UpReg <- toptable_UpReg[toptable_UpReg["Pathways"] == x, ]
            Avg_LogFC<-mean(x_UpReg$logFC)
            Table<-data.frame(Pathway = x, meanlogFC = Avg_LogFC)
            UpReg_Table<-rbind(UpReg_Table, Table)
          }
          
          UpReg_Table<-UpReg_Table[-1,]
          UpReg_Table<-UpReg_Table[order(-UpReg_Table$meanlogFC),]
          
          UpReg_Table$meanlogFC<-as.numeric(UpReg_Table$meanlogFC)
          UpReg_Table$Pathway<-as.character(UpReg_Table$Pathway)
          UpReg_Table$Pathway<-factor(UpReg_Table$Pathway, levels = UpReg_Table$Pathway)
          
          
          
          #output for title and plot 
          output$Diff_UpRegPath_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("Upregulated Pathways Plot", align = "center") })
          
          output$Diff_UpRegPath_Plot<-renderPlotly(ggplotly(ggplot(UpReg_Table, aes(x=reorder(Pathway, -meanlogFC), meanlogFC)) + 
                                                              geom_bar(stat = "identity", fill = c("gray")) +
                                                              labs(x = "Pathway", y = "Mean LogFC")+
                                                              theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
          
          ##DOWNREGULATED PATHWAYS-----
          DownReg_PathList<-unique(toptable_DownReg$Pathways)
          DownReg_Table<-data.frame(Pathway = "del", meanlogFC = 1)
          
          for (x in DownReg_PathList){
            x_DownReg <- toptable_DownReg[toptable_DownReg["Pathways"] == x, ]
            Avg_LogFC<-mean(x_DownReg$logFC)
            Table<-data.frame(Pathway = x, meanlogFC = Avg_LogFC)
            DownReg_Table<-rbind(DownReg_Table, Table)
          }
          
          DownReg_Table<-DownReg_Table[-1,]
          DownReg_Table<-DownReg_Table[order(DownReg_Table$meanlogFC),]
          
          DownReg_Table$meanlogFC<-as.numeric(DownReg_Table$meanlogFC)
          DownReg_Table$Pathway<-as.character(DownReg_Table$Pathway)
          DownReg_Table$Pathway<-factor(DownReg_Table$Pathway, levels = DownReg_Table$Pathway)
          
          output$Diff_DownRegPath_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("Downregulated Pathways Plot", align = "center") })
          output$Diff_DownRegPath_Plot<-renderPlotly(ggplotly(ggplot(DownReg_Table, aes(x=reorder(Pathway, meanlogFC), meanlogFC)) + 
                                                                geom_bar(stat = "identity", fill = c("gray")) +
                                                                labs(x = "Pathway", y = "Mean LogFC")+
                                                                theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
          
          
          
          ##P VALUES PATHWAYS
          PVal_toptable<-toptable_Path[toptable_Path["P.Value"] <0.05, ]
          PVal_PathList<-unique(PVal_toptable$Pathways)
          
          PVal_Table<-data.frame(Pathway = "del", Number_SPM = 1)
          
          for (x in PVal_PathList){
            x_PVal <- PVal_toptable[PVal_toptable["Pathways"] == x, ]
            Table<-data.frame(Pathway = x, Number_SPM = length(x_PVal$Pathways))
            PVal_Table<-rbind(PVal_Table, Table)
          }
          
          PVal_Table<-PVal_Table[-1,]
          PVal_Table<-PVal_Table[order(-PVal_Table$Number_SPM),]
          
          PVal_Table$Number_SPM<-as.numeric(PVal_Table$Number_SPM)
          PVal_Table$Pathway<-as.character(PVal_Table$Pathway)
          
          PVal_Table$Pathway<-factor(PVal_Table$Pathway, levels = PVal_Table$Pathway)
          
          output$Diff_PValPath_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("Significant SPMs Pathways Plot", align = "center") })
          output$Diff_PValPath_Plot<-renderPlotly(ggplotly(ggplot(PVal_Table, aes(x=reorder(Pathway, -Number_SPM), Number_SPM)) + 
                                                             geom_bar(stat = "identity", fill = c("gray")) +
                                                             labs(x = "Pathway", y = "Number of Signficant SPMs") +
                                                             theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        }
      }
      error_Diff <- renderUI(NULL)
    },
    error = function(e) {
      output$error_Diff <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          For more information, please look at the example file available through the 
          download example file button.", e$message
        )
      })
    }
    )
  })
  
  observe({
    
    req(input$LM_ML_File)
    
    tryCatch({
      # Attempt to read the file
      lm_profile_temp <- read.table(input$LM_ML_File$datapath, 
                                    sep = input$sep_LM_ML,
                                    header = TRUE,
                                    stringsAsFactors = FALSE)
      
      # Check if the "Group_names" column exists
      if (!"Group_names" %in% colnames(lm_profile_temp)) {
        output$error_ML_coln <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the format of the file.  
          Please ensure that the format of the uploaded file follows the example file. 
          Please ensure that the names of the columns in your file are named as in the example file. 
          For Example, your first column should contain sample with the column being named Sample_names and 
          the second column should contain at least 2 groups and be named Group_names.
          For more information on the file format you can download an example through the 
          Download Example Button of the sidebar section."
          )
        })
      } 
      # Check for duplicate sample names
      else if (anyDuplicated(lm_profile_temp[, 1])) {
        output$error_ML_Snames <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the content of the file.  
          Please ensure that the first column is named Sample_names and that there are no duplicate values in the rows of that column."
          )
        })
      } 
      else {
        # Clear any previous error messages
        output$error_ML_coln <- renderUI(NULL)
        output$error_ML_Snames <- renderUI(NULL)
        
        # Process the file
        rownames(lm_profile_temp) <- lm_profile_temp[, 1] # Setting the first column as the row names 
        lm_profile <- lm_profile_temp[, -1] # Removing the first column afterwards so it is not duplicated
        
        # Get the column names and update inputs
        Col_NamesML <- unique(as.character(unname(unlist(colnames(lm_profile)))))
        updateSelectInput(session, "Group_ML", choices = Col_NamesML, selected = Col_NamesML[[1]])
        
        # Separate the profiles data by lipid mediators types
        name_List <- unique(unname(unlist(lm_profile[1, -1])))
        name_List <- c("All_LM", name_List)
        
        # Update the names of the models available in the datafile 
        updateSelectInput(session, "RF_Mod_Num", choices = name_List)
        updateSelectInput(session, "RF_Mod_Num1", choices = name_List)
        updateSelectInput(session, "XGB_Mod_Num", choices = name_List)
        updateSelectInput(session, "XGB_Mod_Num1", choices = name_List)
        updateSelectInput(session, "SVM_Mod_Num", choices = name_List)
        updateSelectInput(session, "LA_Mod_Num", choices = name_List)
        updateSelectInput(session, "BC_Mod_Num", choices = name_List)
        
        # Check if the "Group_names" column has enough groups for analysis
        if ("Group_names" %in% colnames(lm_profile)) {
          lm_profile <- lm_profile[-1, ]
          if (length(unique(lm_profile$Group_names)) < 2) {
            output$error_ML_groups <- renderUI({
              div(
                class = "alert alert-danger", role = "alert",
                strong("Error: "), "There appears to be a problem with the content of the file.  
              Please ensure that there are at least two groups available for differential analysis to take place."
              )
            })
          } 
        }
      }
      # If re-run then the error should dissappear
      output$error_ML_coln <- renderUI(NULL)
      
    }, error = function(e) {
      # General error handling
      output$error_ML_coln <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred: ", e$message
        )
      })
    })
    
  })
  
  
  observeEvent(input$Go_ML, {
    

    tryCatch({
    # Call the lipid mediator profiling file variant:
    inFile <- input$LM_ML_File
    
    print(inFile)
    
    if(is.null(inFile))
      
      return(NULL)
    
    set.seed(415) # To get same results even with the random part.
    options(digits = 3) # To get only part of the decimals. 
    
    # Open the file:
      lm_profile_temp = read.table(inFile$datapath, 
                                   sep=input$sep_LM_ML,
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
      
      if (!"Group_names" %in% colnames(lm_profile_temp)) {
        output$error_ML_coln <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the format of the file.  
            Please ensure that the format of the uploaded file follows the example file. 
            Please ensure that the names of the columns in your file are named as in the example file. 
            For Example, your first column should contain sample with the column being named Sample_names and 
            the second column should contain at least 2 groups and be named Group_names.
            For more information on the file format you can download an example through the 
            Download Example Button of the sidebar section."
          )
        })
      } 
      

      
      if (ncol(lm_profile_temp) >= 2 && (length(unique(lm_profile_temp[[1]])) < 2 || length(unique(lm_profile_temp[[2]])) < 2)) {
        output$error_ML_groups <- renderUI({
          div(
            class = "alert alert-danger", role = "alert",
            strong("Error: "), "There appears to be a problem with the content of the file.  
            Please ensure that there are at least two groups available for differential analysis to take place."
          )
        })
      } 
      
      else {
        
        output$error_ML_coln <- renderUI(NULL)
        output$error_ML_Snames <- renderUI(NULL)
        
        rownames(lm_profile_temp) <- lm_profile_temp[, 1] # Setting the first column as the row names 
        lm_profile <- lm_profile_temp[, -1] # Removing the first column afterards so it is not duplicated
        
        
        #rename column names to fit correct name
        
        
        Met_Rename<-c()
        #remove spaces from column names 
        colnames(lm_profile)<-gsub("\\.","", colnames(lm_profile)) 
        
        for (x in 1:ncol(lm_profile)){
          if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
            Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
          }else{
            Met_Rename[x] <-colnames(lm_profile)[x]
          }
        }
        colnames(lm_profile)<-Met_Rename
        ##put coefficents into file!! 
        #list of the unique substrates for the plots 
        name_List <- unique(unname(unlist(lm_profile[1, -1])))
        name_List<-c("All_LM", name_List)
        
        # Separates the profiles data by lipid mediators types:
        
        substrates <- unique(unname(unlist(lm_profile[1, -1])))
        
        
        
        # Creates data frames for each subset:
        
        dataframes_list <- list("ALL LM" = lm_profile[, -1])
        
        for (i in 1:length(substrates)) {
          substrate <- lm_profile[ ,lm_profile[1, ]== substrates[[i]]]
          assign(substrates[[i]], substrate)
          
        }
        
        if (exists('DHA') == TRUE) {dataframes_list[["DHA"]] <- DHA}
        if (exists('n3DPA') == TRUE) {dataframes_list[["n3DPA"]] <- n3DPA}
        if (exists('EPA') == TRUE) {dataframes_list[["EPA"]] <- EPA}
        if (exists('AA') == TRUE) {dataframes_list[["AA"]] <- AA}
        
        
        # Creates the data frame that is going to collect all the info regarding the models:
        
        final_table <- data.frame(machine_learning = "del",
                                  groups = "del",
                                  percentage_accuracy = 1,
                                  sensitivity = 1,
                                  specificity = 1,
                                  TP_per = 1,
                                  FP_per = 1,
                                  TN_per = 1,
                                  FN_per = 1,
                                  stringsAsFactors = FALSE)
        
        #dataframe with top 5 models for XGBoost
        XGB_Model_Table <- data.frame(model = "del",
                                      Accuracy = 1,
                                      specificity = 1,
                                      sensitivity = 1,
                                      TP_per = 1, 
                                      FP_per = 1,
                                      TN_per = 1,
                                      FN_per = 1, 
                                      nrounds = 1,
                                      eta= 1,
                                      max_depth= 1,
                                      gamma= 1,
                                      min_child_weight= 1,
                                      subsample= 1, 
                                      colsample_bytree= 1
        ) 
        
        #selecting tests
        RF_ML <- "RF_ML" %in% input$Model_ML
        XGB_ML <-"XGB_ML" %in% input$Model_ML
        SVM_ML <- "SVM_ML" %in% input$Model_ML
        BC_ML <- "BC_ML" %in% input$Model_ML
        LA_ML <- "LA_ML" %in% input$Model_ML
        
        #progress bar
        #withProgress(message = 'Building Models', style = style, value = 0.1, {
        # Sys.sleep(0.25)})
        
        
        withProgress(message = 'Building Models:', detail = "This may take a while.", value = 0, {
          
          beginning <- Sys.time()
          
          for (j in 1:length(dataframes_list)) {
            
            lm_profiles <- dataframes_list[[j]]
            
            Fatty_Acid<-lm_profiles[1,]
            
            # Save all the values as numeric:
            lm_profile_number <- sapply(lm_profiles[-1, ], function(x) as.numeric(x))
            row.names(lm_profile_number) <- row.names(lm_profiles[-1, ])
            
            # Scale the data because is better when you are running Machine Learning models:
            lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))
            
            # If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors
            # replace the NA for zeros. 
            lm_profiles_scale[is.na(lm_profiles_scale)] <- 0
            
            #--------------------
            ### change 0 to 1/5 of lowest value 
            
            
            lm_profile_zero<-lm_profiles_scale[-1,]
            
            
            
            #get vector the same length as number of columns 
            cols<- 1:(ncol(lm_profile_zero))
            #replace zeros for each column to 1/5 the smallest value for each column 
            
            lm_profile_zero[cols] <- lapply(lm_profile_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
            
            #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
            lm_profile_zero<-as.matrix(lm_profile_zero)
            lm_profile_zero[is.infinite(lm_profile_zero)] <- NA
            min_index <- which.min(lm_profile_zero)
            zero_replace <- (lm_profile_zero[min_index]/5)
            lm_profile_zero <- as.data.frame(lm_profile_zero)
            lm_profile_zero[is.na(lm_profile_zero)] <- zero_replace  
            
            
            
            #add first row to dataframe 
            lm_profiles_scale<-rbind(lm_profiles_scale[1,] ,lm_profile_zero)
            
            #exclude columns with the same value for all rows
            #lm_profiles_scale<-Filter(var, lm_profiles_scale[-1,])
            
            
            
            #---------------------
            # Add the classification variable to the data frame (Responder and non responder):
            ResponseName<-as.character(input$Group_ML)
            lm_profiles_scale$responses <- factor(lm_profile[-1, ResponseName])
            
            # Make sure that column names do not represent a problem to randomForest making them a valid name to R.
            names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))
            
            #---> MACHINE LEARNING (randomForest R):
            if (RF_ML){
              # In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model. 
              oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.
              for (mtry in 1:(ncol(lm_profiles_scale) - 1)) { 
                rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                                      importance = TRUE, ntree = 10000)
                oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
              }
              # Define the best mtry according to the best prediction value. 
              final_mtry <- which.max(oob_error)
              
              # Run the model again with the right mtry value. 
              rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                                   importance = TRUE, ntree = 10000)
              # Get the confusion matrix of the model, sensitivity and specificity: 
              confusion_lm_profiles <- as.data.frame(rf_lm_profiles_final$confusion)
              confusion_lm_profiles[is.na(confusion_lm_profiles)] <- 0
              
              # Calculates sensitivity, specificity and AUC.
              sensitivity_lm_profiles <- confusion_lm_profiles[2, 2]/(confusion_lm_profiles[2, 2] + confusion_lm_profiles[2, 1])
              specificity_lm_profiles <- confusion_lm_profiles[1, 1]/(confusion_lm_profiles[1, 1] + confusion_lm_profiles[1, 2])
              
              # Final table for random forest:
              no_oob_error_table <- data.frame(machine_learning = "RF",
                                               groups = names(dataframes_list)[j],
                                               percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[10000])*100),
                                               sensitivity = sensitivity_lm_profiles,
                                               specificity = specificity_lm_profiles,
                                               TP_per = (1 - confusion_lm_profiles[2, 3])*100,
                                               FP_per = confusion_lm_profiles[1, 3]*100,
                                               TN_per = (1 -confusion_lm_profiles[1, 3])*100,
                                               FN_per = confusion_lm_profiles[2, 3]*100,
                                               stringsAsFactors = FALSE)
              final_table <- rbind(final_table, no_oob_error_table)
              # Number of trees plot:
              tree_table <- data.frame(Ensemble = c(1:10000),
                                       OBB = rf_lm_profiles_final$err.rate[, 1],
                                       AvgAcc = 100 - ((rf_lm_profiles_final$err.rate[, 1])*100))
              ntree_plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) +
                geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
                ggtitle(names(dataframes_list)[j]) +
                scale_x_continuous(name = "Trees") +
                scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
                theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                      axis.title.y = element_text(size = 25, colour = "black"),
                      axis.title.x = element_text(size = 25, colour = "black"),
                      axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                      axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                      legend.position = ("none"),
                      panel.background = element_rect(fill = "white"),
                      axis.ticks.length = unit(0.4, "cm"))
              assign(paste("ntree_", names(dataframes_list)[j], sep = ""), ntree_plot)
              
              #get importance dataframe 
              rf_Importance<-as.data.frame(importance(rf_lm_profiles_final, type = 1))
              rf_Importance$lipid_mediators <- rownames(rf_Importance)
              
              Fatty_Acid <- 1:length(rf_Importance$lipid_mediators)
              cbind(rf_Importance, Fatty_Acid)
              
              
              # #fatty acid points 
              for (i in 1:length(rf_Importance$lipid_mediators)){
                rf_Importance$Fatty_Acid[i]<-SPMs_FA[rf_Importance$lipid_mediators[i]]
              }
              
              #change the names to be in the proper format 
              rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)] <- paste(gsub("^X", "'", rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)]), "'", sep = "")
              
              rf_Importance$lipid_mediators <- gsub("A4'", "A'[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
              rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
              
              rf_Importance$lipid_mediators <- gsub("TXB2", "TXB[2]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("PGD2", "PGD[2]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("PGE2", "PGE[2]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("LXB4", "LXB[4]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("LXA4", "LXA[4]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("LTB4", "LTB[4]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("LTC4", "LTC[4]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("LTE4", "LTE[4]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("LTD4", "LTD[4]", rf_Importance$lipid_mediators)
              
              
              rf_Importance$lipid_mediators <- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required
              rf_Importance$lipid_mediators <- gsub("PGF2a", "PGF[2~a]", rf_Importance$lipid_mediators)
              rf_Importance$lipid_mediators <- gsub("n-3-dpa", "[n-3]~DPA", rf_Importance$lipid_mediators)
              
              #rf_Importance$lipid_mediators <- gsub("4$", "[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
              #rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
              #rf_Importance$lipid_mediators<- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required 
              #rf_Importance$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", rf_Importance$lipid_mediators) # Transform n-3 DPA as a subscript: 
              
              # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
              # (this with the purpose of get the LM in decreasing order in the figure):
              
              rf_Importance <- rf_Importance[order(rf_Importance$MeanDecreaseAccuracy), ]
              
              rf_Importance$lipid_mediators <- factor( rf_Importance$lipid_mediators, 
                                                       levels = rf_Importance$lipid_mediators[
                                                         order(rf_Importance$MeanDecreaseAccuracy)])
              
              lipids <- as.character(rf_Importance$lipid_mediators)
              #add colors to SPM names
              lm_classes <- as.factor(rf_Importance$Fatty_Acid)
              names(lm_classes) <- rf_Importance$lipid_mediators
              
              dha_index<-which((lm_classes=="DHA")== TRUE)
              n_three_index<-which((lm_classes=="n3DPA")== TRUE)
              epa_index<-which((lm_classes=="EPA")== TRUE)
              aa_index<-which((lm_classes=="AA")== TRUE)
              
              lm_colors <- NULL
              lm_colors[dha_index] <- "blue"
              lm_colors[n_three_index] <- "brown"
              lm_colors[epa_index] <- "darkgoldenrod1"
              lm_colors[aa_index] <- "darkslategray"
              
              #rf_Importance1<- subset(rf_Importance, MeanDecreaseAccuracy >1)
              rf_VIP_plot <- ggplot(data = rf_Importance, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = Fatty_Acid)) + geom_point(size = 3) +
                scale_y_continuous(name = "Mean Decrease Accuracy") +
                labs(x = "Lipid Mediators", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
                scale_x_discrete(labels = parse(text = lipids)) + 
                scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                     "EPA" ="darkgoldenrod1", "AA" = "darkslategray"
                )) +
                coord_flip() + 
                theme(axis.title = element_text(size = 20),
                      axis.text.x  = element_text(size = 15, hjust = 0.5),
                      axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                      legend.position = "top",
                      aspect.ratio = 2/1,
                      legend.title = element_text(size = 20),
                      legend.text  = element_text(size = 15),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
              
              assign(paste("rf_VIP_", names(dataframes_list)[j], sep = ""), rf_VIP_plot)
              assign(paste("rf_Mod_", names(dataframes_list)[j], sep = ""), rf_lm_profiles_final)
              
              
            }
            
            #---> MACHINE LEARNING (EXTREME GRADIENT BOOSTING)
            if (XGB_ML){
              df<-lm_profiles_scale
              #need to convert responses to integer of 0 and 1
              
              df$responses<-as.factor(df$responses)
              responses<-df$responses 
              label<-as.integer(df$responses) -1
              
              #remove labels from dataframe
              df$responses<-NULL
              df<-sapply(df, function(x) as.numeric(x))
              
              #
              #split it into training and testing dataset with 75/20
              n = nrow(df) #get number of rows for dataframe
              train.index = sample(n,floor(0.7*n)) #randomize rows to get dataset 
              train.data = as.matrix(df[train.index,])
              train.label = label[train.index]
              test.data = as.matrix(df[-train.index,])
              test.label = label[-train.index]
              
              #transform dataset to xgb matrix 
              xgb.train = xgb.DMatrix(data=train.data,label=train.label)
              xgb.test = xgb.DMatrix(data=test.data, label=test.label)
              
              #create grid with all the options for xgboost
              searchGridSubCol <- expand.grid(subsample = c(0.6, 0.8, 1), 
                                              colsample_bytree = c(0.6, 0.8, 1),
                                              max_depth = c(3, 4, 5),
                                              gamma_val = c(0.5, 1, 1.5, 2, 5), 
                                              eta = c(0.1, 0.2, 0.3),
                                              min_child = c(1, 5, 10)
              )
              
              #run through each combination of the grid 
              system.time(
                ErrorsHyperparameters <- apply(searchGridSubCol, 1, function(parameterList){
                  
                  #Extract Parameters to test
                  currentSubsampleRate <- parameterList[["subsample"]]
                  currentColsampleRate <- parameterList[["colsample_bytree"]]
                  currentDepth <- parameterList[["max_depth"]]
                  currentEta <- parameterList[["eta"]]
                  currentGamma <- parameterList[["gamma_val"]]
                  currentMinChild <- parameterList[["min_child"]]
                  
                  #run selected parameter through the model 
                  xgboostModelCV <- xgb.cv(data =  xgb.train, nrounds = 10000, nfold = 5, showsd = TRUE, 
                                           metrics = "error", verbose = TRUE, "eval_metric" = "error",
                                           "objective" = "binary:logistic", "max.depth" = currentDepth, "eta" = currentEta,       
                                           "subsample" = currentSubsampleRate, "colsample_bytree" = currentColsampleRate
                                           , print_every_n = 10, booster = "gbtree",
                                           early_stopping_rounds = 10, "gamma" = currentGamma, "min_child_weight" = currentMinChild)
                  
                  #have error evaluation score as dataframe 
                  xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
                  test_error <- tail(xvalidationScores$test_error_mean, 1)
                  train_error <- tail(xvalidationScores$train_error_mean,1)
                  output <- return(c(test_error, train_error, currentSubsampleRate, currentColsampleRate, currentDepth, currentEta, currentGamma, currentMinChild))})) #add to final table called output
              
              
              output_Error <- as.data.frame(t(ErrorsHyperparameters))
              
              varnames <- c("TestError", "TrainError", "SubSampRate", "ColSampRate", "Depth", "eta", "gamma", "currentMinChild")
              names(output_Error) <- varnames #change the variable names 
              top_output<-output_Error[order(output_Error$TestError, output_Error$TrainError), ] #order it by the lowest error rate 
              
              #empty table to put the top 10 model information into 
              Final_Acc <- data.frame(model = "test",
                                      Accuracy = 0,
                                      specificity = 0, 
                                      sensitivity = 0,
                                      TP_per = 0, 
                                      FP_per = 0,
                                      TN_per = 0,
                                      FN_per = 0, 
                                      nrounds = 0,
                                      eta= 0,
                                      max_depth= 0,
                                      gamma= 0,
                                      min_child_weight=0,
                                      subsample=0, 
                                      colsample_bytree= 0
              )
              
              top<-1:50
              
              #loop to build models for the parameters with the lowest error rate 
              for (i in top){
                #parameter list with the lowest output 
                params <- list(booster = "gbtree", objective = "binary:logistic", eta= top_output[i,6], gamma= top_output[i,7],
                               max_depth= top_output[i,5], min_child_weight=top_output[i,8], subsample=top_output[i,3], 
                               colsample_bytree= top_output[i,4])
                #Build model based on the training data 
                xgb<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
                
                #run model on test data = see how good it is at predicting
                pred<-predict(xgb, test.data)
                
                #obtain confusion matrix 
                pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1 
                pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
                pred<-as.numeric(pred) #makes it numeric
                test.label_num <-as.numeric(test.label) #makes it numeric
                response_List<-unique(responses)#get list of the 2 conditions 
                pred_y <- cut(pred, 2, label = c(response_List[1], response_List[2])) #converts to factor of responder or non responder 
                test_y <-cut(test.label_num, 2, label = c(response_List[1], response_List[2])) #convers to factor of responder or non responder 
                ConMatrix <- confusionMatrix(reference = test_y, data = pred_y) #build confusion matrix 
                TP = as.numeric(ConMatrix$table[2,2])
                FP = as.numeric(ConMatrix$table[2,1])
                TN = as.numeric(ConMatrix$table[1,1])
                FN = as.numeric(ConMatrix$table[1,2])
                
                #convert build table with accuracy for each model 
                Table <- data.frame(model = names(dataframes_list)[j],
                                    Accuracy = as.numeric(ConMatrix$overall[1])*100,
                                    specificity = as.numeric(ConMatrix$byClass[2]),
                                    sensitivity = as.numeric(ConMatrix$byClass[1]),
                                    TP_per = (TP/(TP+FN)*100), 
                                    FP_per = (FP/(FP+TN)*100),
                                    TN_per = (TN/(FP+TN)*100),
                                    FN_per = (FN/(TP+FN)*100), 
                                    nrounds = 10000,
                                    eta= top_output[i,6],
                                    max_depth= top_output[i,5],
                                    gamma= top_output[i,7],
                                    min_child_weight=top_output[i,8],
                                    subsample=top_output[i,3], 
                                    colsample_bytree= top_output[i,4]
                )
                Final_Acc<-rbind(Final_Acc, Table)
                
              }
              
              #remove the first row of random values 
              Final_Acc<-Final_Acc[-1,]
              
              #order the top 10 models based on highest accuracy, sensitivity, and specificity 
              Top_Acc<-Final_Acc[with(Final_Acc, order(-Accuracy, -sensitivity, -specificity)), ]
              
              # model<- c(paste(names(dataframes_list)[j],"_Model_1", sep = ""), 
              #            paste(names(dataframes_list)[j],"_Model_2", sep = ""),
              #            paste(names(dataframes_list)[j],"_Model_3", sep = ""),
              #            paste(names(dataframes_list)[j],"_Model_4", sep = ""),
              #            paste(names(dataframes_list)[j],"_Model_5", sep = ""))
              
              Top_Acc<-Top_Acc[1:5,]
              
              Top_Acc1 <- Top_Acc
              
              
              top5<-1:5
              for (i in top5){
                Top_Acc1[i,1] <-paste(names(dataframes_list)[j],"_Model", i,  sep = "")
              }
              
              
              # XGBoost Table:
              XGBoost_table <- data.frame(machine_learning = "XGB",
                                          groups = names(dataframes_list)[j],
                                          percentage_accuracy =Top_Acc[1,2],
                                          sensitivity = Top_Acc[1,3],
                                          specificity = Top_Acc[1,4],
                                          TP_per = Top_Acc[1,5], 
                                          FP_per = Top_Acc[1,6],
                                          TN_per = Top_Acc[1,7],
                                          FN_per = Top_Acc[1,8], 
                                          stringsAsFactors = FALSE)
              
              final_table <- rbind(final_table, XGBoost_table)
              
              XGB_Model_Table <- rbind(XGB_Model_Table, Top_Acc1)
              
              
              
              #preparing nrounds plots for top 5 models 
              iterations<-1:10000
              top_table<- data.frame(iterations)
              
              for (i in top5){
                
                #paramters for the i model 
                params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[i,10], gamma= Top_Acc[i,12], 
                               max_depth= Top_Acc[i,11], min_child_weight= Top_Acc[i,13], subsample= Top_Acc[i,14])
                
                #run the cv to get the error rate for each nround 
                xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = 10000, nfold = 5, showsd = T, 
                                 stratified = T, early_stop_round = 20, maximize = F, metrics = "error", verbose = FALSE)
                
                Table1<-data.frame( test = as.numeric(xgbcv$evaluation_log$test_error_mean))
                colnames(Table1)[1] <-paste(names(dataframes_list)[j],"_Model_", i, sep = "")
                #colnames(Table1)[2] <-paste(names(dataframes_list)[j],"Test_model_", i, sep = "")
                top_table<-cbind(top_table, Table1) #add to table
                
              }
              
              #plot using plotly 
              nround_plot<-plot_ly(data = top_table, x = ~iterations) 
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[2])), name = colnames(top_table)[2], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[3])), name = colnames(top_table)[3], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[4])), name = colnames(top_table)[4], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[5])), name = colnames(top_table)[5], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[6])), name = colnames(top_table)[6], mode = 'lines')
              nround_plot <- nround_plot %>% layout(title = paste(names(dataframes_list)[j], "Plot", sep =" "), xaxis = list(title = 'Number of Rounds'), 
                                                    yaxis = list(title = 'Test Error Rate'))
              #makes download quality of plot better
              
              nround_plot <- nround_plot%>% config(toImageButtonOptions = list(format = "jpeg", width = 1500, height = 750))
              
              assign(paste("nround_", names(dataframes_list)[j], sep = ""), nround_plot)
              
              
              
              #importance plot 
              
              #build model for the top model 
              params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[1,10], gamma= Top_Acc[1,12], 
                             max_depth= Top_Acc[1,11], min_child_weight= Top_Acc[1,13], subsample= Top_Acc[1,14], colsample_bytree= Top_Acc[1,15])
              xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
              
              
              XGB_Importance = xgb.importance(model = xgb_Mod1) #important matrix
              
              Fatty_Acid <- 1:length(XGB_Importance$Feature)
              cbind(XGB_Importance, Fatty_Acid)
              
              #fatty acid points 
              for (i in 1:length(XGB_Importance$Feature)){
                XGB_Importance$Fatty_Acid[i]<-SPMs_FA[XGB_Importance$Feature[i]]
              }
              
              # #change the names to be in the proper format
              XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)] <- paste(gsub("^X", "'", XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)]), "'", sep = "")
              
              XGB_Importance$Feature <- gsub("A4'", "A'[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 5
              XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature) #B4 converted to subscript 4
              
              XGB_Importance$Feature <- gsub("TXB2", "TXB[2]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("PGD2", "PGD[2]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("PGE2", "PGE[2]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("LXB4", "LXB[4]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("LXA4", "LXA[4]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("LTB4", "LTB[4]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("LTC4", "LTC[4]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("LTE4", "LTE[4]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("LTD4", "LTD[4]", XGB_Importance$Feature)
              
              
              XGB_Importance$Feature <- gsub("\\.", "-", XGB_Importance$Feature) # Replace "." with "-" when required
              XGB_Importance$Feature <- gsub("PGF2a", "PGF[2~a]", XGB_Importance$Feature)
              XGB_Importance$Feature <- gsub("n-3-dpa", "[n-3]~DPA", XGB_Importance$Feature)
              
              #XGB_Importance$Feature <- gsub("4$", "[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 4
              #XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature ) #B4 converted to subscript 4
              #XGB_Importance$Feature  <- gsub("\\.", "-", XGB_Importance$Feature ) # Replace "." with "-" when required
              #XGB_Importance$Feature <- gsub("n.3.DPA", "[n-3~DPA]", XGB_Importance$Feature ) # Transform n-3 DPA as a subscript:
              
              
              # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
              # (this with the purpose of get the LM in decreasing order in the figure):
              
              XGB_Importance <- XGB_Importance[order(XGB_Importance$Gain), ]
              
              XGB_Importance$Feature<- factor(XGB_Importance$Feature, 
                                              levels = XGB_Importance$Feature[
                                                order(XGB_Importance$Gain)])
              lipids <- as.character(XGB_Importance$Feature )
              
              #add colors to x axis labels
              lm_classes <- as.factor(XGB_Importance$Fatty_Acid)
              names(lm_classes) <- XGB_Importance$Feature
              
              dha_index<-which((lm_classes=="DHA")== TRUE)
              n_three_index<-which((lm_classes=="n3DPA")== TRUE)
              epa_index<-which((lm_classes=="EPA")== TRUE)
              aa_index<-which((lm_classes=="AA")== TRUE)
              
              lm_colors <- NULL
              lm_colors[dha_index] <- "blue"
              lm_colors[n_three_index] <- "brown"
              lm_colors[epa_index] <- "darkgoldenrod1"
              lm_colors[aa_index] <- "darkslategray"
              
              XGB_VIP_Plot <- ggplot(data = XGB_Importance, mapping = aes(x = Feature, y = Gain, color = Fatty_Acid)) + geom_point(size = 3) +
                scale_y_continuous(name = "Gain") +
                labs(x = "Features", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
                scale_x_discrete(labels = parse(text = lipids)) +
                scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                     "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
                coord_flip() +
                theme(axis.title = element_text(size = 20),
                      axis.text.x  = element_text(size = 15, hjust = 0.5),
                      axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                      legend.position = "top",
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
              
              
              #save VIP plot
              assign(paste("XGB_VIP_", names(dataframes_list)[j], sep = ""), XGB_VIP_Plot)
              
              #Save Xgboost model 
              assign(paste("xgb_Mod_", names(dataframes_list)[j], sep = ""), xgb_Mod1)
            }
            #---> MACHINE LEARNING (Classyfire R): 
            
            # Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
            # and creates a novelty detection for the creation of the model. 
            
            # The idea is to create several models and see which one fits the best. The models will be based on the whole
            # lipid profiles and the different groups based on substrates. 
            
            # "cfBuild" to create the SVM:
            # Clasyfire requieres matrix: 
            if (SVM_ML){
              #convert dataframe to matrix
              lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])
              
              #building SVM model
              support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses, 
                                                  bootNum = 70,ensNum = 70, cpus = 4) 
              #obtaining confusion matrix 
              conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale))
              
              # SVM table:
              Support_vector_table <- data.frame(machine_learning = "SVM",
                                                 groups = names(dataframes_list)[j],
                                                 percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                                 sensitivity = conf_matrix[4, 3]/100,
                                                 specificity = conf_matrix[1, 3]/100,
                                                 TP_per = conf_matrix[4, 3],
                                                 FP_per = conf_matrix[2, 3],
                                                 TN_per = conf_matrix[1, 3],
                                                 FN_per = conf_matrix[3, 3],
                                                 stringsAsFactors = FALSE)
              final_table <- rbind(final_table, Support_vector_table)
              
              # Ensemble Plot:
              ensAcc   <- getAcc(support_lmprofiles_scale)$Test
              meanVal  <- ensAcc[1]
              for (i in 2:length(ensAcc)) {
                meanVal <- c(meanVal, mean(ensAcc[1:i]))
              } 
              ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc), 
                                          AvgAcc = meanVal)
              ensemble_plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) + 
                geom_point(aes(colour = AvgAcc), size = 5) + 
                geom_line(linetype = "dotted", size = 1) +
                ggtitle(names(dataframes_list)[j]) +
                scale_x_continuous(name = "Ensemble interaction") +
                scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
                theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                      axis.title.y = element_text(size = 25, colour = "black"),
                      axis.title.x = element_text(size = 25, colour = "black"),
                      axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                      axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                      legend.position = ("none"),
                      panel.background = element_rect(fill = "white"),
                      axis.ticks.length = unit(0.4, "cm"))
              assign(paste("ensemble_", names(dataframes_list)[j], sep = ""), ensemble_plot)
              
              #Save SVM model 
              assign(paste("svm_Mod_", names(dataframes_list)[j], sep = ""), support_lmprofiles_scale)
            }
            
            #---> ELASTIC NET REGRESSION (caret R):
            if (LA_ML){
              # Get the explanatory variables as a matrix:
              explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
              
              # LASSO Analysis:
              model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet", 
                                 trControl = trainControl("boot", number = 70))
              # Get the confusion matrix of the model:
              conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
              
              # Final Elastic net model table:
              net_table <- data.frame(machine_learning = "GLMNET",
                                      groups = names(dataframes_list)[j],
                                      percentage_accuracy = max(model_net$results$Accuracy)*100,
                                      sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                                      specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                                      TP_per = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                      FP_per = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                      TN_per = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                      FN_per = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                      stringsAsFactors = FALSE)
              final_table <- rbind(final_table, net_table)
              
              # Parameter tuning figure: 
              scaleFUN <- function(x) sprintf("%.2f", x)
              
              boot_net_plot <- ggplot(model_net, highlight = TRUE) +
                scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
                scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
                ggtitle(names(dataframes_list)[j]) +
                scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
                scale_shape_manual(values=c(16, 16, 16)) +
                labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
                guides(shape = FALSE) +
                theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                      axis.title.y = element_text(size = 25, colour = "black"),
                      axis.title.x = element_text(size = 25, colour = "black"),
                      axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
                      axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep.
                      legend.title = element_text(size = 25),
                      legend.text  = element_text(size = 20),
                      legend.key = element_rect(fill = NA),
                      legend.key.size = unit(1.3, "cm"), 
                      panel.background = element_rect(fill = "white"),
                      axis.ticks.length = unit(0.4, "cm"))
              assign(paste("net_", names(dataframes_list)[j], sep = ""), boot_net_plot)
              
              #rename and save elastic net model 
              assign(paste("la_Mod_", names(dataframes_list)[j], sep = ""), model_net)
            }
            
            #---> BAYESIAN MODEL (Caret R): 
            
            if (BC_ML){
              explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
              
              bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm", 
                                trControl = trainControl("boot", number = 70))
              conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)
              
              # Final Bayesian model table:
              bay_table <- data.frame(machine_learning = "BAYES",
                                      groups = names(dataframes_list)[j],
                                      percentage_accuracy = max(bayesian$results$Accuracy)*100,
                                      sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                                      specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                                      TP_per = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                      FP_per = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                      TN_per = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                      FN_per = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                      stringsAsFactors = FALSE)
              final_table <- rbind(final_table, bay_table)
              
              #rename and save bayesian model 
              assign(paste("bc_Mod_", names(dataframes_list)[j], sep = ""), bayesian)
            }
            
            incProgress(1/length(dataframes_list), detail = paste(names(dataframes_list)[j]))
          }
          
          final_table <- final_table[-1, ]
          table <- data.frame(Model = factor(final_table$groups, levels = unique(final_table$groups)),
                              methodology = final_table$machine_learning,
                              accuracy = round(final_table$percentage_accuracy, 0))
          accuracy <- ggplot(data = table, aes(x = Model, y = accuracy, fill = methodology)) +
            geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black", size = 2) +
            scale_fill_manual(values=c("dodgerblue2","firebrick2",'goldenrod1', 'lightslategray', "purple")) + 
            geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
                      aes(label = paste(accuracy, "%", sep = "")), size = 8) +
            scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20),
                               expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
            coord_cartesian(ylim = c(1, 120)) +
            theme(axis.title = element_text(size = 40),
                  axis.title.x = element_blank(),
                  axis.text.x  =  element_text(size = 40, hjust = 0.5, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.0), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.0), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_blank(),
                  legend.position = "top",
                  legend.key.size = unit(1.3, "cm"), 
                  legend.text  = element_text(size = 25),
                  legend.spacing.x = unit(1, "cm"),
                  axis.ticks.length = unit(0.4, "cm"))
          
          #Put figures and models together together: 
          if (RF_ML){
            rf_plot_list <- list("ALL LM" = `ntree_ALL LM`)
            rf_VIPplot_list <- list("ALL LM" = `rf_VIP_ALL LM`)
            rf_Mod_list<-list("ALL LM" = `rf_Mod_ALL LM`)
          }
          if (XGB_ML){
            XGB_VIPplot_list <- list("ALL LM" = `XGB_VIP_ALL LM`)
            xgb_Mod_list<-list("ALL LM" = `xgb_Mod_ALL LM`)}
          if (SVM_ML){
            svm_plot_list <- list("ALL LM" = `ensemble_ALL LM`)
            svm_Mod_list<-list("ALL LM" = `svm_Mod_ALL LM`)
          }
          if (LA_ML){
            net_plot_list <- list("ALL LM" = `net_ALL LM`)
            la_Mod_list<-list("ALL LM" = `la_Mod_ALL LM`)}
          if (BC_ML){bc_Mod_list<-list("ALL LM" = `bc_Mod_ALL LM`)}
          
          
          column <- 2
          if (exists('DHA') == TRUE) {
            if (RF_ML){
              rf_plot_list[["DHA"]] <- ntree_DHA
              rf_VIPplot_list[["DHA"]] <- rf_VIP_DHA
              rf_Mod_list[["DHA"]] <- rf_Mod_DHA
            }
            
            if (XGB_ML){
              XGB_VIPplot_list[["DHA"]] <- XGB_VIP_DHA
              xgb_Mod_list[["DHA"]] <- xgb_Mod_DHA
            }
            if (SVM_ML){
              svm_plot_list[["DHA"]] <- ensemble_DHA
              svm_Mod_list[["DHA"]] <- svm_Mod_DHA}
            if (LA_ML){
              net_plot_list[["DHA"]] <- net_DHA
              la_Mod_list[["DHA"]] <- la_Mod_DHA}
            if (BC_ML){bc_Mod_list[["DHA"]] <- bc_Mod_DHA}
            column <- column + 1}
          
          if (exists('n3DPA') == TRUE) {
            if (RF_ML){
              rf_plot_list[["n3DPA"]] <- ntree_n3DPA
              rf_VIPplot_list[["n3DPA"]] <- rf_VIP_n3DPA
              rf_Mod_list[["n3DPA"]] <- rf_Mod_n3DPA
            }
            if (XGB_ML){
              XGB_VIPplot_list[["n3DPA"]] <- XGB_VIP_n3DPA
              xgb_Mod_list[["n3DPA"]] <- xgb_Mod_n3DPA}
            if (SVM_ML){
              svm_plot_list[["n3DPA"]] <- ensemble_n3DPA
              svm_Mod_list[["n3DPA"]] <- svm_Mod_n3DPA}
            if (LA_ML){
              net_plot_list[["n3DPA"]] <- net_n3DPA
              la_Mod_list[["n3DPA"]] <- la_Mod_n3DPA}
            if (BC_ML){bc_Mod_list[["n3DPA"]] <- bc_Mod_n3DPA}
            
            column <- column + 1}
          
          if (exists('EPA') == TRUE) {
            if (RF_ML){
              rf_plot_list[["EPA"]] <- ntree_EPA
              rf_VIPplot_list[["EPA"]] <- rf_VIP_EPA
              rf_Mod_list[["EPA"]] <- rf_Mod_EPA
            }
            if (XGB_ML){
              XGB_VIPplot_list[["EPA"]] <- XGB_VIP_EPA
              xgb_Mod_list[["EPA"]] <- xgb_Mod_EPA}
            if (SVM_ML){
              svm_plot_list[["EPA"]] <- ensemble_EPA
              svm_Mod_list[["EPA"]] <- svm_Mod_EPA}
            if (LA_ML){
              net_plot_list[["EPA"]] <- net_EPA
              la_Mod_list[["EPA"]] <- la_Mod_EPA}
            if (BC_ML){bc_Mod_list[["EPA"]] <- bc_Mod_EPA}
            
            column <- column + 1}
          
          if (exists('AA') == TRUE) {
            if (RF_ML){
              rf_plot_list[["AA"]] <- ntree_AA
              rf_VIPplot_list[["AA"]] <- rf_VIP_AA
              rf_Mod_list[["AA"]] <- rf_Mod_AA}
            if (XGB_ML){
              XGB_VIPplot_list[["AA"]] <- XGB_VIP_AA
              xgb_Mod_list[["AA"]] <- xgb_Mod_AA}
            if (SVM_ML){
              svm_plot_list[["AA"]] <- ensemble_AA
              svm_Mod_list[["AA"]] <- svm_Mod_AA}
            if (LA_ML){
              net_plot_list[["AA"]] <- net_AA
              la_Mod_list[["AA"]] <- la_Mod_AA}
            if (BC_ML){bc_Mod_list[["AA"]] <- bc_Mod_AA}
            
            column <- column + 1}
          if (column >= 6) {column <- 4}
          
          #if (RF_ML){grid.arrange(grobs = rf_plot_list, ncol = column/2)}
          else if (SVM_ML){grid.arrange(grobs = svm_plot_list, ncol = column/2)}
          else if (LA_ML){grid.arrange(grobs = net_plot_list, ncol = column/2)}
          
          # Final table otputs :
          final_table$`% Accuracy Score` <- round(final_table$percentage_accuracy, 0)
          final_table$Sensitivity <- round(final_table$sensitivity, 2)
          final_table$Specificity <- round(final_table$specificity, 2)
          final_table$TP <- round(final_table$TP_per, 0)
          final_table$FP <- round(final_table$FP_per, 0)
          final_table$TN <- round(final_table$TN_per, 0)
          final_table$FN <- round(final_table$FN_per, 0)
          
          final_table <- final_table[, c(1, 2, 10:16)]
          colnames(final_table)[1] <- "Machine Learning Methodology"
          colnames(final_table)[2] <- "Model"
          
          
          ##outputs for accuracy tab 
          output$AccPlot_Title <- renderUI({req(input$Go_ML); h2("% Accuracy Score Figure for different ML models", align = "center") })
          output$AccTable_Title <- renderUI({req(input$Go_ML); h2("Model Table Summary", align = "center") })
          
          output$Accuracy_ML <- renderPlot({return(accuracy)})
          output$downloadAcc_ML_Plot <- downloadHandler(filename = function(){paste("Accuracy_Plot",'.png',sep='')},
                                                        content = function(file){
                                                          ggsave(file,plot= accuracy, width = 15, height = 10)})
          
          output$ML_Table <- renderDataTable({return(final_table)})
          output$downloadAcc_ML_Table <- downloadHandler(filename = function(){"Accuracy_Table.csv"}, 
                                                         content = function(fname){
                                                           write.csv(final_table, fname)})
          
          
          #outputs if Random forest model is run
          if (RF_ML){
            #ui outputs to get the page names 
            output$RF_Plot_Title <- renderUI({req(input$Go_ML); h2("Optimal Parameters RandomForest", align = "center") })
            output$RF_VIPPlot_Title <- renderUI({req(input$Go_ML); h2("Importance of Variance Plot", align = "center") })
            
            #random forest optimal parameters plots output and download button 
            output$RF_Plot <- renderPlot({return(grid.arrange(grobs = rf_plot_list, ncol = 1))})
            
            output$RF_Plot1 <- renderPlot({ rf_plot_list[[1]] })
            output$RF_Plot2 <- renderPlot({ rf_plot_list[[2]] })
            output$RF_Plot3 <- renderPlot({ rf_plot_list[[3]] })
            output$RF_Plot4 <- renderPlot({ rf_plot_list[[4]] })
            output$RF_Plot5 <- renderPlot({ rf_plot_list[[5]] })
            
            rf_plot_grid<-grid.arrange(grobs = rf_plot_list, ncol = 1)
            output$download_RF_Plot <- downloadHandler(filename = function(){paste("Random_Forest_Plot",input$RF_Mod_Num1,'.png',sep='')},
                                                       content = function(file){
                                                         ggsave(file, plot = rf_plot_list[[match(input$RF_Mod_Num1, name_List)]], width = 12, height = 12)})
            
            #random forest importance plot output and download button
            output$RF_VIPPlot <- renderPlot({return(grid.arrange(grobs = rf_VIPplot_list, ncol = 2))})
            
            output$RF_VIPPlot1 <- renderPlot({ rf_VIPplot_list[[1]] })
            output$RF_VIPPlot2 <- renderPlot({ rf_VIPplot_list[[2]] })
            output$RF_VIPPlot3 <- renderPlot({ rf_VIPplot_list[[3]] })
            output$RF_VIPPlot4 <- renderPlot({ rf_VIPplot_list[[4]] })
            output$RF_VIPPlot5 <- renderPlot({ rf_VIPplot_list[[5]] })
            
            
            output$download_RF_VIPPlot <- downloadHandler(filename = function(){paste("Random_Forest_VIP_Plot",input$RF_Mod_Num,'.png',sep='')},
                                                          content = function(file){
                                                            ggsave(file,plot= rf_VIPplot_list[[match(input$RF_Mod_Num, name_List)]], width = 12, height = 12)})
            
            #downloading model information for paramters tab and importance tab
            output$download_RF_Mod <- downloadHandler(filename = function(){paste("Random_Forest_Model_",input$RF_Mod_Num,'.rds',sep='')},
                                                      content = function(fname){
                                                        saveRDS(rf_Mod_list[[match(input$RF_Mod_Num, name_List)]], fname)})
            
            output$download_RF_Mod1 <- downloadHandler(filename = function(){paste("Random_Forest_Model_",input$RF_Mod_Num1,'.rds',sep='')},
                                                       content = function(fname){
                                                         saveRDS(rf_Mod_list[[2]], compress = TRUE, fname)})
            
          }
          
          #outputs for xgboost tabs 
          if (XGB_ML){
            #ui output page titles - only appears once action button clicked
            output$XGB_Plot_Title <- renderUI({req(input$Go_ML); h2("Top 5 models for XGBoost Test Error Rate", align = "center") })
            output$XGB_VIPPlot_Title <- renderUI({req(input$Go_ML); h2("Importance of Variance Plot", align = "center") })
            
            
            #Parameters plots and tablefor xgboost parameters page 
            output$XGB_Plot1 <- renderPlotly({return(`nround_ALL LM`)})
            if (exists('DHA') == TRUE){output$XGB_Plot2 <- renderPlotly({return(nround_DHA)})}
            if (exists('n3DPA') == TRUE) {output$XGB_Plot3 <- renderPlotly({return(nround_n3DPA)})}
            if (exists('EPA') == TRUE) {output$XGB_Plot4 <- renderPlotly({return(nround_EPA)})}
            if (exists('AA') == TRUE) {output$XGB_Plot5 <- renderPlotly({return(nround_AA)})}
            XGB_Model_Table<-XGB_Model_Table[-1,]
            XGB_Model_Table[,-1] <-round(XGB_Model_Table[,-1],2)
            output$XGB_Models<-renderDataTable(XGB_Model_Table, rownames = FALSE,
                                               options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
            output$download_XGB_Table <- downloadHandler(filename = function(){"XGB_Model_Table.csv"}, 
                                                         content = function(fname){
                                                           write.csv(XGB_Model_Table, fname)})
            
            #xgboost importance plots output 
            output$XGB_VIPPlot <- renderPlot({return(grid.arrange(grobs = XGB_VIPplot_list, ncol = 2))})
            
            output$XGB_VIPPlot1 <- renderPlot({ XGB_VIPplot_list[[1]] })
            output$XGB_VIPPlot2 <- renderPlot({ XGB_VIPplot_list[[2]] })
            output$XGB_VIPPlot3 <- renderPlot({ XGB_VIPplot_list[[3]] })
            output$XGB_VIPPlot4 <- renderPlot({ XGB_VIPplot_list[[4]] })
            output$XGB_VIPPlot5 <- renderPlot({ XGB_VIPplot_list[[5]] })
            
            # Downloading the VIP Plots as PNG files:
            output$download_XGB_VIPPlot <- downloadHandler(filename = function(){paste("XGBoost_VIP_Plot",input$XGB_Mod_Num1,'.png',sep='')},
                                                           content = function(file){
                                                             ggsave(file,plot= XGB_VIPplot_list[[match(input$XGB_Mod_Num1, name_List)]], width = 12, height = 12)})
            # Downloading the XGB Plots as PNG files:
            output$download_XGB_Plot <- downloadHandler(filename = function(){paste("XGBoost_Plot",input$XGB_Mod_Num,'.png',sep='')},
                                                        content = function(file){
                                                          ggsave(file,plot= XGB_VIPplot_list[[match(input$XGB_Mod_Num, name_List)]], width = 12, height = 12)})
            
            
            # Downloading model information for parameters tab:
            output$download_XGB_Mod <- downloadHandler(filename = function(){paste("XGBoost_Model_",input$XGB_Mod_Num,'.rds',sep='')},
                                                       content = function(fname){
                                                         saveRDS(xgb_Mod_list[[match(input$XGB_Mod_Num, name_List)]], fname)})
            # Downloading model information from the importance tab:
            output$download_XGB_Mod1 <- downloadHandler(filename = function(){paste("XGBoost_Model_",input$XGB_Mod_Num1,'.rds',sep='')},
                                                        content = function(fname){
                                                          saveRDS(xgb_Mod_list[[match(input$XGB_Mod_Num1, name_List)]], fname)})
            
          }
          
          
          #all the outputs for svm model 
          if (SVM_ML){
            
            #ui output for SVM title page 
            output$SVM_Plot_Title <- renderUI({req(input$Go_ML); h2("Optimal Parameters SVM", align = "center") })
            
            #plot output and download for svm optima parameters 
            output$SVM_Plot <- renderPlot({return(grid.arrange(grobs = svm_plot_list, ncol = 1))})
            
            output$SVM_Plot1 <- renderPlot({ svm_plot_list[[1]] })
            output$SVM_Plot2 <- renderPlot({ svm_plot_list[[2]] })
            output$SVM_Plot3 <- renderPlot({ svm_plot_list[[3]] })
            output$SVM_Plot4 <- renderPlot({ svm_plot_list[[4]] })
            output$SVM_Plot5 <- renderPlot({ svm_plot_list[[5]] })
            
            
            svm_plot_grid<-grid.arrange(grobs = svm_plot_list, ncol = 1)
            output$download_SVM_Plot <- downloadHandler(filename = function(){paste("SVM_Ensemble_Plot",input$SVM_Mod_Num,'.png',sep='')},
                                                        content = function(file){
                                                          ggsave(file, plot = svm_plot_list[[match(input$SVM_Mod_Num, name_List)]], width = 12, height = 12)})
            
            #download button for svm model
            output$download_SVM_Mod <- downloadHandler(filename = function(){paste("SVM_Model_",input$SVM_Mod_Num,'.rds',sep='')},
                                                       content = function(fname){
                                                         saveRDS(svm_Mod_list[[match(input$SVM_Mod_Num, name_List)]], fname)})
            
          }
          
          #outputs for elastic net regression model 
          if (LA_ML){
            #Elastic Net title 
            output$LA_Plot_Title <- renderUI({req(input$Go_ML); h2("Optimal Parameters Elastic Net Regression", align = "center") })
            
            #plot output for elastic net regression 
            output$LA_Plot <- renderPlot({return(grid.arrange(grobs = net_plot_list, ncol = 1))})
            
            output$LA_Plot1 <- renderPlot({ net_plot_list[[1]] })
            output$LA_Plot2 <- renderPlot({ net_plot_list[[2]] })
            output$LA_Plot3 <- renderPlot({ net_plot_list[[3]] })
            output$LA_Plot4 <- renderPlot({ net_plot_list[[4]] })
            output$LA_Plot5 <- renderPlot({ net_plot_list[[5]] })
            
            net_plot_grid<-grid.arrange(grobs = net_plot_list, ncol = 1)
            output$download_LA_Plot <- downloadHandler(filename = function(){paste("Elasti_Net_Plot",input$LA_Mod_Num,'.png',sep='')},
                                                       content = function(file){
                                                         ggsave(file, plot = net_plot_list[[match(input$LA_Mod_Num, name_List)]], width = 12, height = 12)})
            
            #download button for LA model
            output$download_LA_Mod <- downloadHandler(filename = function(){paste("ElasticNet_Model_",input$LA_Mod_Num,'.rds',sep='')},
                                                      content = function(fname){
                                                        saveRDS(la_Mod_list[[match(input$LA_Mod_Num, name_List)]], fname)})
          }
          
          #output for bayesian model 
          if (BC_ML){
            
            #download button for BC model
            output$download_BC_Mod <- downloadHandler(filename = function(){paste("BayesianGLM_Model_",input$BC_Mod_Num,'.rds',sep='')},
                                                      content = function(fname){
                                                        saveRDS(bc_Mod_list[[match(input$BC_Mod_Num, name_List)]], fname)})
          }
          end <- Sys.time()
          
          time_taken <- difftime(end, beginning, units = "secs")
          formatted_time <- format_time(time_taken)
          
          # Get CPU usage (number of cores used)
          cpu_cores <- detectCores(logical = FALSE)
          
          # Get memory usage
          memory_usage <- mem_used()
          
        })
        
        output$system_info <- renderUI({
          HTML(paste0(
            "<b>Time Taken:</b> ", formatted_time, " seconds<br>",
            "<b>CPU Cores Used:</b> ", cpu_cores, "<br>",
            "<b>Memory Used:</b> ", round(memory_usage / (1024 ^ 2), 2), " MB<br>"
          ))
        })
        
      }
      output$error_ML_coln <- renderUI(NULL)
    },
    error = function(e) {
      # Global error handling
      output$error_ML_coln <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred, 
          please ensure that your file follows the appropriate format.
          For more information please check the available file from the download example file button.", e$message
        )
      })
    }
    )
  })
  
  observeEvent(input$Build_ML, {
    
    tryCatch({
      # Call the lipid mediator profiling file variant:
      inFile <- input$LM_ML_File2
      
      print(inFile)
      
      if(is.null(inFile))
        
        return(NULL)
      
      
      
      set.seed(415) # To get same results even with the random part.
      options(digits = 3) # To get only part of the decimals. 
      
      # Open the file:
      
      
      lm_profile = read.table(inFile$datapath, 
                              sep=input$sep_LM_ML2,
                              header = TRUE,
                              row.names = 1,
                              stringsAsFactors = FALSE)
      
      Met_Rename<-c()
      #remove spaces from column names 
      colnames(lm_profile)<-gsub("\\.","", colnames(lm_profile)) 
      
      for (x in 1:ncol(lm_profile)){
        if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
        }else{
          Met_Rename[x] <-colnames(lm_profile)[x]
        }
      }
      colnames(lm_profile)<-Met_Rename
      
      # Separates the profiles data by lipid mediators types:
      
      substrates <- unique(unname(unlist(lm_profile[1, -1])))
      
      # Creates data frames for each subset:
      
      dataframes_list <- list("ALL LM" = lm_profile[, -1])
      
      for (i in 1:length(substrates)) {
        substrate <- lm_profile[ ,lm_profile[1, ]== substrates[[i]]]
        assign(substrates[[i]], substrate)
        
      }
      
      if (exists('DHA') == TRUE) {dataframes_list[["DHA"]] <- DHA}
      if (exists('n3DPA') == TRUE) {dataframes_list[["n3DPA"]] <- n3DPA}
      if (exists('EPA') == TRUE) {dataframes_list[["EPA"]] <- EPA}
      if (exists('AA') == TRUE) {dataframes_list[["AA"]] <- AA}
      
      
      # Creates the data frame that is going to collect all the info regarding the models:
      
      final_table <- data.frame(machine_learning = "del",
                                groups = "del",
                                percentage_accuracy = 1,
                                sensitivity = 1,
                                specificity = 1,
                                TP_per = 1,
                                FP_per = 1,
                                TN_per = 1,
                                FN_per = 1,
                                stringsAsFactors = FALSE)
      #dataframe with top 5 models for XGBoost
      XGB_Model_Table <- data.frame(model = "del",
                                    Accuracy = 1,
                                    specificity = 1,
                                    sensitivity = 1,
                                    TP_per = 1, 
                                    FP_per = 1,
                                    TN_per = 1,
                                    FN_per = 1, 
                                    nrounds = 1,
                                    eta= 1,
                                    max_depth= 1,
                                    gamma= 1,
                                    min_child_weight= 1,
                                    subsample= 1, 
                                    colsample_bytree= 1)
      
      
      #progress bar
      #withProgress(message = 'Building Models', style = style, value = 0.1, {
      # Sys.sleep(0.25)})
      
      
      withProgress(message = 'Building Models:', detail = "This may take a while.", value = 0, {
        
        for (j in 1:length(dataframes_list)) {
          
          lm_profiles <- dataframes_list[[j]]
          
          # Save all the values as numeric:
          lm_profile_number <- sapply(lm_profiles[-1, ], function(x) as.numeric(x))
          row.names(lm_profile_number) <- row.names(lm_profiles[-1, ])
          
          # Scale the data because is better when you are running Machine Learning models:
          lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))
          
          # If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors
          # replace the NA for zeros. 
          lm_profiles_scale[is.na(lm_profiles_scale)] <- 0
          
          #--------------------
          ### change 0 to 1/5 of lowest value 
          
          
          lm_profile_zero<-lm_profiles_scale[-1,]
          
          #get vector the same length as number of columns 
          cols<- 1:(ncol(lm_profile_zero))
          #replace zeros for each column to 1/5 the smallest value for each column 
          
          lm_profile_zero[cols] <- lapply(lm_profile_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
          
          #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
          lm_profile_zero<-as.matrix(lm_profile_zero)
          lm_profile_zero[is.infinite(lm_profile_zero)] <- NA
          min_index <- which.min(lm_profile_zero)
          zero_replace <- (lm_profile_zero[min_index]/5)
          lm_profile_zero <- as.data.frame(lm_profile_zero)
          lm_profile_zero[is.na(lm_profile_zero)] <- zero_replace  
          
          #add first row to dataframe 
          lm_profiles_scale<-rbind(lm_profiles_scale[1,] ,lm_profile_zero)
          
          #exclude columns with the same value for all rows
          #lm_profile_scale<-Filter(var, lm_profile_scale[-1,])
          
          #---------------------
          # Add the classification variable to the data frame (Responder and non responder):
          lm_profiles_scale$responses <- factor(lm_profile[-1, ]$groups)
          
          # Make sure that column names do not represent a problem to randomForest making them a valid name to R.
          names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))
          
          #---> MACHINE LEARNING (randomForest R):
          if (input$RF_Build == TRUE){
            # In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model. 
            oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.
            for (mtry in 1:(ncol(lm_profiles_scale) - 1)) { 
              rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                                    importance = TRUE, ntree = input$RF_ntrees)
              oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[input$RF_ntrees])*100)
            }
            # Define the best mtry according to the best prediction value. 
            final_mtry <- which.max(oob_error)
            
            # Run the model again with the right mtry value. 
            rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                                 importance = TRUE, ntree = input$RF_ntrees)
            # Get the confusion matrix of the model, sensitivity and specificity: 
            confusion_lm_profiles <- as.data.frame(rf_lm_profiles_final$confusion)
            confusion_lm_profiles[is.na(confusion_lm_profiles)] <- 0
            
            # Calculates sensitivity, specificity and AUC.
            sensitivity_lm_profiles <- confusion_lm_profiles[2, 2]/(confusion_lm_profiles[2, 2] + confusion_lm_profiles[2, 1])
            specificity_lm_profiles <- confusion_lm_profiles[1, 1]/(confusion_lm_profiles[1, 1] + confusion_lm_profiles[1, 2])
            
            # Final table for random forest:
            no_oob_error_table <- data.frame(machine_learning = "RF",
                                             groups = names(dataframes_list)[j],
                                             percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[input$RF_ntrees])*100),
                                             sensitivity = sensitivity_lm_profiles,
                                             specificity = specificity_lm_profiles,
                                             TP_per = (1 - confusion_lm_profiles[2, 3])*100,
                                             FP_per = confusion_lm_profiles[1, 3]*100,
                                             TN_per = (1 -confusion_lm_profiles[1, 3])*100,
                                             FN_per = confusion_lm_profiles[2, 3]*100,
                                             stringsAsFactors = FALSE)
            final_table <- rbind(final_table, no_oob_error_table)
            # Number of trees plot:
            tree_table <- data.frame(Ensemble = c(1:input$RF_ntrees),
                                     OBB = rf_lm_profiles_final$err.rate[, 1],
                                     AvgAcc = 100 - ((rf_lm_profiles_final$err.rate[, 1])*100))
            ntree_plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) +
              geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
              ggtitle(names(dataframes_list)[j]) +
              scale_x_continuous(name = "Trees") +
              scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
              theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 25, colour = "black"),
                    axis.title.x = element_text(size = 25, colour = "black"),
                    axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                    legend.position = ("none"),
                    panel.background = element_rect(fill = "white"),
                    axis.ticks.length = unit(0.4, "cm"))
            assign(paste("ntree_", names(dataframes_list)[j], sep = ""), ntree_plot)
            
            
            #get importance dataframe 
            rf_Importance<-as.data.frame(importance(rf_lm_profiles_final, type = 1))
            rf_Importance$lipid_mediators <- rownames(rf_Importance)
            
            Fatty_Acid <- 1:length(rf_Importance$lipid_mediators)
            cbind(rf_Importance, Fatty_Acid)
            
            
            # #fatty acid points 
            for (i in 1:length(rf_Importance$lipid_mediators)){
              rf_Importance$Fatty_Acid[i]<-SPMs_FA[rf_Importance$lipid_mediators[i]]
            }
            
            #change the names to be in the proper format 
            rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)] <- paste(gsub("^X", "'", rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)]), "'", sep = "")
            
            rf_Importance$lipid_mediators <- gsub("A4'", "A'[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
            rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
            
            rf_Importance$lipid_mediators <- gsub("TXB2", "TXB[2]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("PGD2", "PGD[2]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("PGE2", "PGE[2]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("LXB4", "LXB[4]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("LXA4", "LXA[4]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("LTB4", "LTB[4]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("LTC4", "LTC[4]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("LTE4", "LTE[4]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("LTD4", "LTD[4]", rf_Importance$lipid_mediators)
            
            
            rf_Importance$lipid_mediators <- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required
            rf_Importance$lipid_mediators <- gsub("PGF2a", "PGF[2~a]", rf_Importance$lipid_mediators)
            rf_Importance$lipid_mediators <- gsub("n-3-dpa", "[n-3]~DPA", rf_Importance$lipid_mediators)
            
            #rf_Importance$lipid_mediators <- gsub("4$", "[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
            #rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
            #rf_Importance$lipid_mediators<- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required 
            #rf_Importance$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", rf_Importance$lipid_mediators) # Transform n-3 DPA as a subscript: 
            
            # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
            # (this with the purpose of get the LM in decreasing order in the figure):
            
            rf_Importance <- rf_Importance[order(rf_Importance$MeanDecreaseAccuracy), ]
            
            rf_Importance$lipid_mediators <- factor( rf_Importance$lipid_mediators, 
                                                     levels = rf_Importance$lipid_mediators[
                                                       order(rf_Importance$MeanDecreaseAccuracy)])
            
            lipids <- as.character(rf_Importance$lipid_mediators)
            #add colors to SPM names
            lm_classes <- as.factor(rf_Importance$Fatty_Acid)
            names(lm_classes) <- rf_Importance$lipid_mediators
            
            dha_index<-which((lm_classes=="DHA")== TRUE)
            n_three_index<-which((lm_classes=="n3DPA")== TRUE)
            epa_index<-which((lm_classes=="EPA")== TRUE)
            aa_index<-which((lm_classes=="AA")== TRUE)
            
            lm_colors <- NULL
            lm_colors[dha_index] <- "blue"
            lm_colors[n_three_index] <- "brown"
            lm_colors[epa_index] <- "darkgoldenrod1"
            lm_colors[aa_index] <- "darkslategray"
            
            #rf_Importance1<- subset(rf_Importance, MeanDecreaseAccuracy >1)
            rf_VIP_plot <- ggplot(data = rf_Importance, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = Fatty_Acid)) + geom_point(size = 3) +
              scale_y_continuous(name = "Mean Decrease Accuracy") +
              labs(x = "Lipid Mediators", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
              scale_x_discrete(labels = parse(text = lipids)) + 
              scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                   "EPA" ="darkgoldenrod1", "AA" = "darkslategray"
              )) +
              coord_flip() + 
              theme(axis.title = element_text(size = 20),
                    axis.text.x  = element_text(size = 15, hjust = 0.5),
                    axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                    legend.position = "top",
                    aspect.ratio = 2/1,
                    legend.title = element_text(size = 20),
                    legend.text  = element_text(size = 15),
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
            
            
            assign(paste("rf_VIP_", names(dataframes_list)[j], sep = ""), rf_VIP_plot)
            
            
          }
          
          #---> MACHINE LEARNING (EXTREME GRADIENT BOOSTING)
          if (input$XGB_Build == TRUE){
            df<-lm_profiles_scale
            #need to convert responses to integer of 0 and 1
            
            df$responses<-as.factor(df$responses)
            responses<-df$responses
            label<-as.integer(df$responses) -1
            
            #remove labels from dataframe
            df$responses<-NULL
            df<-sapply(df, function(x) as.numeric(x))
            
            #
            #split it into training and testing dataset with 75/20
            n = nrow(df) #get number of rows for dataframe
            train.index = sample(n,floor(0.7*n)) #randomize rows to get dataset
            train.data = as.matrix(df[train.index,])
            train.label = label[train.index]
            test.data = as.matrix(df[-train.index,])
            test.label = label[-train.index]
            
            #transform dataset to xgb matrix
            xgb.train = xgb.DMatrix(data=train.data,label=train.label)
            xgb.test = xgb.DMatrix(data=test.data, label=test.label)
            
            #create grid with all the options for xgboost
            searchGridSubCol <- expand.grid(subsample = c(0.5, 0.75, 1),
                                            colsample_bytree = c(0.5, 0.75, 1),
                                            max_depth = c(2, 4, 6, 8, 10),
                                            gamma_val = c(0, 1, 5, 10),
                                            eta = c(0.01, 0.1, 0.2, 0.3),
                                            min_child = seq(1)
            )
            
            #run through each combination of the grid
            system.time(
              ErrorsHyperparameters <- apply(searchGridSubCol, 1, function(parameterList){
                
                #Extract Parameters to test
                currentSubsampleRate <- parameterList[["subsample"]]
                currentColsampleRate <- parameterList[["colsample_bytree"]]
                currentDepth <- parameterList[["max_depth"]]
                currentEta <- parameterList[["eta"]]
                currentGamma <- parameterList[["gamma_val"]]
                currentMinChild <- parameterList[["min_child"]]
                
                #run selected parameter through the model
                xgboostModelCV <- xgb.cv(data =  xgb.train, nrounds = input$XGB_nrounds, nfold = 5, showsd = TRUE,
                                         metrics = "error", verbose = TRUE, "eval_metric" = "error",
                                         "objective" = "binary:logistic", "max.depth" = currentDepth, "eta" = currentEta,
                                         "subsample" = currentSubsampleRate, "colsample_bytree" = currentColsampleRate
                                         , print_every_n = 10, booster = "gbtree",
                                         early_stopping_rounds = 10, "gamma" = currentGamma, "min_child_weight" = currentMinChild)
                
                #have error evaluation score as dataframe
                xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
                test_error <- tail(xvalidationScores$test_error_mean, 1)
                train_error <- tail(xvalidationScores$train_error_mean,1)
                output <- return(c(test_error, train_error, currentSubsampleRate, currentColsampleRate, currentDepth, currentEta, currentGamma, currentMinChild))})) #add to final table called output
            
            
            output_Error <- as.data.frame(t(ErrorsHyperparameters))
            
            varnames <- c("TestError", "TrainError", "SubSampRate", "ColSampRate", "Depth", "eta", "gamma", "currentMinChild")
            names(output_Error) <- varnames #change the variable names
            top_output<-output_Error[order(output_Error$TestError, output_Error$TrainError), ] #order it by the lowest error rate
            
            #empty table to put the top 10 model information into
            Final_Acc <- data.frame(model = "test",
                                    Accuracy = 0,
                                    specificity = 0,
                                    sensitivity = 0,
                                    TP_per = 0,
                                    FP_per = 0,
                                    TN_per = 0,
                                    FN_per = 0,
                                    nrounds = 0,
                                    eta= 0,
                                    max_depth= 0,
                                    gamma= 0,
                                    min_child_weight=0,
                                    subsample=0,
                                    colsample_bytree= 0
            )
            
            top<-1:50
            
            #loop to build models for the parameters with the lowest error rate
            for (i in top){
              #parameter list with the lowest output
              params <- list(booster = "gbtree", objective = "binary:logistic", eta= top_output[i,6], gamma= top_output[i,7],
                             max_depth= top_output[i,5], min_child_weight=top_output[i,8], subsample=top_output[i,3],
                             colsample_bytree= top_output[i,4])
              #Build model based on the training data
              xgb<-xgb.train(data = xgb.train, params = params, nrounds = input$XGB_nrounds, eval.metric = "error", early.stop.rounds=10,  silent = 0)
              
              #run model on test data = see how good it is at predicting
              pred<-predict(xgb, test.data)
              
              #obtain confusion matrix
              pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1
              pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
              pred<-as.numeric(pred) #makes it numeric
              test.label_num <-as.numeric(test.label) #makes it numeric
              pred_y <- cut(pred, 2, label = c("Non-responder", "Responder")) #converts to factor of responder or non responder
              test_y <-cut(test.label_num, 2, label = c("Non-responder", "Responder")) #convers to factor of responder or non responder
              ConMatrix <- confusionMatrix(reference = test_y, data = pred_y) #build confusion matrix
              TP = as.numeric(ConMatrix$table[2,2])
              FP = as.numeric(ConMatrix$table[2,1])
              TN = as.numeric(ConMatrix$table[1,1])
              FN = as.numeric(ConMatrix$table[1,2])
              
              #convert build table with accuracy for each model
              Table <- data.frame(model = names(dataframes_list)[j],
                                  Accuracy = as.numeric(ConMatrix$overall[1])*100,
                                  specificity = as.numeric(ConMatrix$byClass[2]),
                                  sensitivity = as.numeric(ConMatrix$byClass[1]),
                                  TP_per = (TP/(TP+FN)*100),
                                  FP_per = (FP/(FP+TN)*100),
                                  TN_per = (TN/(FP+TN)*100),
                                  FN_per = (FN/(TP+FN)*100),
                                  nrounds = input$XGB_nrounds,
                                  eta= top_output[i,6],
                                  max_depth= top_output[i,5],
                                  gamma= top_output[i,7],
                                  min_child_weight=top_output[i,8],
                                  subsample=top_output[i,3],
                                  colsample_bytree= top_output[i,4]
              )
              Final_Acc<-rbind(Final_Acc, Table)
              
            }
            
            #remove the first row of random values
            Final_Acc<-Final_Acc[-1,]
            
            #order the top 10 models based on highest accuracy, sensitivity, and specificity
            Top_Acc<-Final_Acc[with(Final_Acc, order(-Accuracy, -sensitivity, -specificity)), ]
            
            Top_Acc<-Top_Acc[1:5,]
            
            Top_Acc1 <- Top_Acc
            top5<-1:5
            for (i in top5){
              Top_Acc1[i,1] <-paste(names(dataframes_list)[j],"_Model", i,  sep = "")
            }
            
            #selecting the model number based on the users input 
            if ( names(dataframes_list)[j] == "ALL LM"){
              Model_Num<-input$XGB_ALL_LM_Model
            }
            if ( names(dataframes_list)[j] == "DHA"){
              Model_Num<-input$XGB_DHA_Model
            }
            if ( names(dataframes_list)[j] == "n3DPA"){
              Model_Num<-input$XGB_n3DPA_Model
            }
            if ( names(dataframes_list)[j] == "EPA"){
              Model_Num<-input$XGB_EPA_Model
            }
            if ( names(dataframes_list)[j] == "AA"){
              Model_Num<-input$XGB_AA_Model
            }
            
            # XGBoost Table:
            XGBoost_table <- data.frame(machine_learning = "XGB",
                                        groups = names(dataframes_list)[j],
                                        percentage_accuracy =Top_Acc[Model_Num,2],
                                        sensitivity = Top_Acc[Model_Num,3],
                                        specificity = Top_Acc[Model_Num,4],
                                        TP_per = Top_Acc[Model_Num,5],
                                        FP_per = Top_Acc[Model_Num,6],
                                        TN_per = Top_Acc[Model_Num,7],
                                        FN_per = Top_Acc[Model_Num,8],
                                        stringsAsFactors = FALSE)
            
            
            
            final_table <- rbind(final_table, XGBoost_table)
            XGB_Model_Table <- rbind(XGB_Model_Table, Top_Acc1)
            
            #preparing nrounds plots for top 5 models
            iterations<-1:input$XGB_nrounds
            top_table<- data.frame(iterations)
            
            for (i in top5){
              
              #paramters for the i model
              params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[i,10], gamma= Top_Acc[i,12],
                             max_depth= Top_Acc[i,11], min_child_weight= Top_Acc[i,13], subsample= Top_Acc[i,14])
              
              #run the cv to get the error rate for each nround
              xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = input$XGB_nrounds, nfold = 5, showsd = T,
                               stratified = T, early_stop_round = 20, maximize = F, metrics = "error", verbose = FALSE)
              # xgb.name1<-paste0("Train_model_", i, sep = "")
              # xgb.name2<-paste0("Test_model_", i, sep = "")
              # assign(paste(names(dataframes_list)[j], "Train_model_", i,sep = ""), xgb.name1)
              # assign(paste(names(dataframes_list)[j], "Test_model_", i ,sep = ""), xgb.name2)
              #Table1<-data.frame(train = as.numeric(xgbcv$evaluation_log$train_error_mean), test = as.numeric(xgbcv$evaluation_log$test_error_mean))
              Table1<-data.frame( test = as.numeric(xgbcv$evaluation_log$test_error_mean))
              colnames(Table1)[1] <-paste(names(dataframes_list)[j],"_Model_", i, sep = "")
              #colnames(Table1)[2] <-paste(names(dataframes_list)[j],"Test_model_", i, sep = "")
              top_table<-cbind(top_table, Table1) #add to table
              
            }
            
            #plot using plotly
            nround_plot<-plot_ly(data = top_table, x = ~iterations)
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[2])), name = colnames(top_table)[2], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[3])), name = colnames(top_table)[3], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[4])), name = colnames(top_table)[4], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[5])), name = colnames(top_table)[5], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[6])), name = colnames(top_table)[6], mode = 'lines')
            nround_plot <- nround_plot %>% layout(title = paste(names(dataframes_list)[j], "Plot", sep =" "), xaxis = list(title = 'Number of Rounds'),
                                                  yaxis = list(title = 'Test Error Rate'))
            
            #makes download quality of plot better
            nround_plot <- nround_plot%>% config(toImageButtonOptions = list(format = "jpeg", width = 1500, height = 750))
            
            assign(paste("nround_", names(dataframes_list)[j], sep = ""), nround_plot)
            
            #importance plot
            #build model for the top model 
            params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[Model_Num,10], gamma= Top_Acc[Model_Num,12], 
                           max_depth= Top_Acc[Model_Num,11], min_child_weight= Top_Acc[Model_Num, 13], subsample= Top_Acc[Model_Num, 14], colsample_bytree= Top_Acc[Model_Num, 15])
            xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = input$XGB_nrounds, eval.metric = "error", early.stop.rounds=10,  silent = 0)
            
            
            XGB_Importance = xgb.importance(model = xgb_Mod1) #important matrix
            
            Fatty_Acid <- 1:length(XGB_Importance$Feature)
            cbind(XGB_Importance, Fatty_Acid)
            
            #fatty acid points 
            for (i in 1:length(XGB_Importance$Feature)){
              XGB_Importance$Fatty_Acid[i]<-SPMs_FA[XGB_Importance$Feature[i]]
            }
            
            # #change the names to be in the proper format
            XGB_Importance$Feature[grep("^X", XGB_Importance$Feature)] <- paste(gsub("^X", "'", XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)]), "'", sep = "")
            
            XGB_Importance$Feature <- gsub("A4'", "A'[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 5
            XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature) #B4 converted to subscript 4
            
            XGB_Importance$Feature <- gsub("TXB2", "TXB[2]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("PGD2", "PGD[2]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("PGE2", "PGE[2]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("LXB4", "LXB[4]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("LXA4", "LXA[4]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("LTB4", "LTB[4]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("LTC4", "LTC[4]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("LTE4", "LTE[4]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("LTD4", "LTD[4]", XGB_Importance$Feature)
            
            
            XGB_Importance$Feature <- gsub("\\.", "-", XGB_Importance$Feature) # Replace "." with "-" when required
            XGB_Importance$Feature <- gsub("PGF2a", "PGF[2~a]", XGB_Importance$Feature)
            XGB_Importance$Feature <- gsub("n-3-dpa", "[n-3]~DPA", XGB_Importance$Feature)
            
            #XGB_Importance$Feature <- gsub("4$", "[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 4
            #XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature ) #B4 converted to subscript 4
            #XGB_Importance$Feature  <- gsub("\\.", "-", XGB_Importance$Feature ) # Replace "." with "-" when required
            #XGB_Importance$Feature <- gsub("n.3.DPA", "[n-3~DPA]", XGB_Importance$Feature ) # Transform n-3 DPA as a subscript:
            
            
            # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
            # (this with the purpose of get the LM in decreasing order in the figure):
            
            XGB_Importance <- XGB_Importance[order(XGB_Importance$Gain), ]
            
            XGB_Importance$Feature<- factor(XGB_Importance$Feature, 
                                            levels = XGB_Importance$Feature[
                                              order(XGB_Importance$Gain)])
            lipids <- as.character(XGB_Importance$Feature )
            
            #add colors to x axis labels
            lm_classes <- as.factor(XGB_Importance$Fatty_Acid)
            names(lm_classes) <- XGB_Importance$Feature
            
            dha_index<-which((lm_classes=="DHA")== TRUE)
            n_three_index<-which((lm_classes=="n3DPA")== TRUE)
            epa_index<-which((lm_classes=="EPA")== TRUE)
            aa_index<-which((lm_classes=="AA")== TRUE)
            
            lm_colors <- NULL
            lm_colors[dha_index] <- "blue"
            lm_colors[n_three_index] <- "brown"
            lm_colors[epa_index] <- "darkgoldenrod1"
            lm_colors[aa_index] <- "darkslategray"
            
            XGB_VIP_Plot <- ggplot(data = XGB_Importance, mapping = aes(x = Feature, y = Gain, color = Fatty_Acid)) + geom_point(size = 3) +
              scale_y_continuous(name = "Gain") +
              labs(x = "Features", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
              scale_x_discrete(labels = parse(text = lipids)) +
              scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                   "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
              coord_flip() +
              theme(axis.title = element_text(size = 20),
                    axis.text.x  = element_text(size = 15, hjust = 0.5),
                    axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                    legend.position = "top",
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
            
            # #build model for the top model
            # params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[Model_Num, 10], gamma= Top_Acc[Model_Num, 12],
            #                max_depth= Top_Acc[Model_Num, 11], min_child_weight= Top_Acc[Model_Num, 13], subsample= Top_Acc[Model_Num, 14])
            # xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = input$XGB_nrounds, eval.metric = "error", early.stop.rounds=10,  silent = 0)
            # 
            # #get the importance informamation
            # XGB_Importance = xgb.importance(model = xgb_Mod1)
            # #
            # XGB_Importance$Features <- XGB_Importance$Features
            # #
            # 
            # XGB_VIP_Plot<-xgb.ggplot.importance(XGB_Importance, rel_to_first = TRUE, n_clusters = 1:2, xlab = "Relative importance")
            # XGB_VIP_Plot<-XGB_VIP_Plot + ggplot2::labs(y = "Frequency",  title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
            #   scale_fill_discrete(name = "Class", labels = c("Non_Responder", "Responsder"))
            assign(paste("XGB_VIP_", names(dataframes_list)[j], sep = ""), XGB_VIP_Plot)
          }
          #---> MACHINE LEARNING (Classyfire R):
          
          # Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
          # and creates a novelty detection for the creation of the model.
          
          # The idea is to create several models and see which one fits the best. The models will be based on the whole
          # lipid profiles and the different groups based on substrates.
          
          # "cfBuild" to create the SVM:
          # Clasyfire requieres matrix:
          if (input$SVM_Build == TRUE){
            #convert dataframe to matrix
            lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])
            
            #building SVM model
            support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses,
                                                bootNum = input$SVM_BootNum, ensNum = input$SVM_Ensemble, cpus = 4)
            #obtaining confusion matrix
            conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale))
            
            # SVM table:
            Support_vector_table <- data.frame(machine_learning = "SVM",
                                               groups = names(dataframes_list)[j],
                                               percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                               sensitivity = conf_matrix[4, 3]/100,
                                               specificity = conf_matrix[1, 3]/100,
                                               TP_per = conf_matrix[4, 3],
                                               FP_per = conf_matrix[2, 3],
                                               TN_per = conf_matrix[1, 3],
                                               FN_per = conf_matrix[3, 3],
                                               stringsAsFactors = FALSE)
            final_table <- rbind(final_table, Support_vector_table)
            
            # Ensemble Plot:
            ensAcc   <- getAcc(support_lmprofiles_scale)$Test
            meanVal  <- ensAcc[1]
            for (i in 2:length(ensAcc)) {
              meanVal <- c(meanVal, mean(ensAcc[1:i]))
            }
            ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc),
                                        AvgAcc = meanVal)
            ensemble_plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) +
              geom_point(aes(colour = AvgAcc), size = 5) +
              geom_line(linetype = "dotted", size = 1) +
              ggtitle(names(dataframes_list)[j]) +
              scale_x_continuous(name = "Ensemble interaction") +
              scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
              theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 25, colour = "black"),
                    axis.title.x = element_text(size = 25, colour = "black"),
                    axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep.
                    legend.position = ("none"),
                    panel.background = element_rect(fill = "white"),
                    axis.ticks.length = unit(0.4, "cm"))
            assign(paste("ensemble_", names(dataframes_list)[j], sep = ""), ensemble_plot)
          }
          
          #---> ELASTIC NET REGRESSION (caret R):
          if (input$LA_Build == TRUE){
            # Get the explanatory variables as a matrix:
            explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
            
            # LASSO Analysis:
            model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet",
                               trControl = trainControl("boot", number = input$LA_BootNum))
            # Get the confusion matrix of the model:
            conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
            
            # Final Elastic net model table:
            net_table <- data.frame(machine_learning = "GLMNET",
                                    groups = names(dataframes_list)[j],
                                    percentage_accuracy = max(model_net$results$Accuracy)*100,
                                    sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                                    specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                                    TP_per = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                    FP_per = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                    TN_per = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                    FN_per = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                    stringsAsFactors = FALSE)
            final_table <- rbind(final_table, net_table)
            
            # Parameter tuning figure:
            scaleFUN <- function(x) sprintf("%.2f", x)
            
            boot_net_plot <- ggplot(model_net, highlight = TRUE) +
              scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
              scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
              ggtitle(names(dataframes_list)[j]) +
              scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
              scale_shape_manual(values=c(16, 16, 16)) +
              labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
              guides(shape = FALSE) +
              theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 25, colour = "black"),
                    axis.title.x = element_text(size = 25, colour = "black"),
                    axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep.
                    legend.title = element_text(size = 25),
                    legend.text  = element_text(size = 20),
                    legend.key = element_rect(fill = NA),
                    legend.key.size = unit(1.3, "cm"),
                    panel.background = element_rect(fill = "white"),
                    axis.ticks.length = unit(0.4, "cm"))
            assign(paste("net_", names(dataframes_list)[j], sep = ""), boot_net_plot)
          }
          
          #---> BAYESIAN MODEL (Caret R):
          
          if (input$BC_Build == TRUE){
            explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
            
            bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm",
                              trControl = trainControl("boot", number = input$BC_BootNum))
            conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)
            
            # Final Bayesian model table:
            bay_table <- data.frame(machine_learning = "BAYES",
                                    groups = names(dataframes_list)[j],
                                    percentage_accuracy = max(bayesian$results$Accuracy)*100,
                                    sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                                    specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                                    TP_per = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                    FP_per = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                    TN_per = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                    FN_per = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                    stringsAsFactors = FALSE)
            final_table <- rbind(final_table, bay_table)
          }
          
          incProgress(1/length(dataframes_list), detail = paste(names(dataframes_list)[j]))
        }
        
        final_table <- final_table[-1, ]
        table <- data.frame(Model = factor(final_table$groups, levels = unique(final_table$groups)),
                            methodology = final_table$machine_learning,
                            accuracy = round(final_table$percentage_accuracy, 0))
        accuracy <- ggplot(data = table, aes(x = Model, y = accuracy, fill = methodology)) +
          geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black", size = 2) +
          scale_fill_manual(values=c("dodgerblue2","firebrick2",'goldenrod1', 'lightslategray', "purple")) + 
          geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
                    aes(label = paste(accuracy, "%", sep = "")), size = 8) +
          scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20),
                             expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
          coord_cartesian(ylim = c(1, 120)) +
          theme(axis.title = element_text(size = 40),
                axis.title.x = element_blank(),
                axis.text.x  =  element_text(size = 40, hjust = 0.5, colour = "black"), # Put color to the labels
                axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
                axis.line = element_line(colour = 'black', size = 1.0), # Color and thickness of axis
                axis.ticks = element_line(colour = "black", size = 1.0), # Color and thickness of every axis sep. 
                panel.background = element_rect(fill = "white"),
                legend.title = element_blank(),
                legend.position = "top",
                legend.key.size = unit(1.3, "cm"), 
                legend.text  = element_text(size = 25),
                legend.spacing.x = unit(1, "cm"),
                axis.ticks.length = unit(0.4, "cm"))
        
        #Put figures together: 
        if (input$RF_Build == TRUE){
          rf_plot_list <- list("ALL LM" = `ntree_ALL LM`)
          rf_VIPplot_list <- list("ALL LM" = `rf_VIP_ALL LM`)
        }
        if (input$SVM_Build == TRUE){svm_plot_list <- list("ALL LM" = `ensemble_ALL LM`)}
        if (input$LA_Build == TRUE){net_plot_list <- list("ALL LM" = `net_ALL LM`)}
        if (input$XGB_Build == TRUE){XGB_VIPplot_list <- list("ALL LM" = `XGB_VIP_ALL LM`)}
        
        column <- 2
        if (exists('DHA') == TRUE) {
          if (input$RF_Build == TRUE){
            rf_plot_list[["DHA"]] <- ntree_DHA
            rf_VIPplot_list[["DHA"]] <- rf_VIP_DHA
          }
          if (input$SVM_Build == TRUE){svm_plot_list[["DHA"]] <- ensemble_DHA}
          if (input$LA_Build == TRUE){net_plot_list[["DHA"]] <- net_DHA}
          if (input$XGB_Build == TRUE){XGB_VIPplot_list[["DHA"]] <- XGB_VIP_DHA}
          column <- column + 1}
        if (exists('n3DPA') == TRUE) {
          if (input$RF_Build == TRUE){
            rf_plot_list[["n3DPA"]] <- ntree_n3DPA
            rf_VIPplot_list[["n3DPA"]] <- rf_VIP_n3DPA
          }
          if (input$SVM_Build == TRUE){svm_plot_list[["n3DPA"]] <- ensemble_n3DPA}
          if (input$LA_Build == TRUE){net_plot_list[["n3DPA"]] <- net_n3DPA}
          if (input$XGB_Build == TRUE){XGB_VIPplot_list[["n3DPA"]] <- XGB_VIP_n3DPA}
          column <- column + 1}
        if (exists('EPA') == TRUE) {
          if (input$RF_Build == TRUE){
            rf_plot_list[["EPA"]] <- ntree_EPA
            rf_VIPplot_list[["EPA"]] <- rf_VIP_EPA
          }
          if (input$SVM_Build == TRUE){svm_plot_list[["EPA"]] <- ensemble_EPA}
          if (input$LA_Build == TRUE){net_plot_list[["EPA"]] <- net_EPA}
          if (input$XGB_Build == TRUE){XGB_VIPplot_list[["EPA"]] <- XGB_VIP_EPA}
          column <- column + 1}
        if (exists('AA') == TRUE) {
          if (input$RF_Build == TRUE){
            rf_plot_list[["AA"]] <- ntree_AA
            rf_VIPplot_list[["AA"]] <- rf_VIP_AA}
          if (input$SVM_Build == TRUE){svm_plot_list[["AA"]] <- ensemble_AA}
          if (input$LA_Build == TRUE){net_plot_list[["AA"]] <- net_AA}
          if (input$XGB_Build == TRUE){XGB_VIPplot_list[["AA"]] <- XGB_VIP_AA}
          column <- column + 1}
        if (column >= 6) {column <- 4}
        
        #if (RF_ML){grid.arrange(grobs = rf_plot_list, ncol = column/2)}
        else if (input$SVM_Build == TRUE){grid.arrange(grobs = svm_plot_list, ncol = column/2)}
        else if (input$LA_Build == TRUE){grid.arrange(grobs = net_plot_list, ncol = column/2)}
        
        # Final table cute:
        final_table$`% Accuracy Score` <- round(final_table$percentage_accuracy, 0)
        final_table$Sensitivity <- round(final_table$sensitivity, 2)
        final_table$Specificity <- round(final_table$specificity, 2)
        final_table$TP <- round(final_table$TP_per, 0)
        final_table$FP <- round(final_table$FP_per, 0)
        final_table$TN <- round(final_table$TN_per, 0)
        final_table$FN <- round(final_table$FN_per, 0)
        
        final_table <- final_table[, c(1, 2, 10:16)]
        colnames(final_table)[1] <- "Machine Learning Methodology"
        colnames(final_table)[2] <- "Model"
        
        if (input$RF_Build == TRUE){output$Build_AccPlot_Title <- renderUI({req(input$Build_ML); h2("% Accuracy Score Figure for different ML models", align = "center") })}
        if (input$RF_Build == TRUE){output$Build_RF_Plot_Title <- renderUI({req(input$Build_ML); h2("Optimal Parameters RandomForest", align = "center") })}
        if (input$RF_Build == TRUE){output$Build_RF_VIPPlot_Title <- renderUI({req(input$Build_ML); h2("Importance of Variance Plot", align = "center") })}
        if (input$SVM_Build == TRUE){output$Build_SVM_Plot_Title <- renderUI({req(input$Build_ML); h2("Optimal Parameters SVM", align = "center") })}
        if (input$XGB_Build == TRUE){
          output$Build_XGB_Plot_Title <- renderUI({req(input$Build_ML); h2("Top 5 models for XGBoost test error rate", align = "center") })
          output$Build_XGB_VIPPlot_Title <- renderUI({req(input$Build_ML); h2("Importance of Variance Plot", align = "center") })
        }
        output$Build_LA_Plot_Title <- renderUI({req(input$Build_ML); h2("Optimal Parameters Elastic Net", align = "center") })
        output$Build_AccTable_Title <- renderUI({req(input$Build_ML); h2("Model Table Summary", align = "center") })
        output$Build_Accuracy_ML <- renderPlot({return(accuracy)})
        output$downloadBuild_Acc_ML_Plot <- downloadHandler(filename = function(){paste("Accuracy_Plot",'.png',sep='')},
                                                            content = function(file){
                                                              ggsave(file,plot= accuracy, width = 15, height = 10)})
        
        if (input$RF_Build == TRUE){
          output$Build_RF_Plot <- renderPlot({return(grid.arrange(grobs = rf_plot_list, ncol = 1))})
          
          # Building the Random Forests plots individually for the Optimisation section:
          output$Build_RF_Plot1 <- renderPlot({ rf_plot_list[[1]] })
          output$Build_RF_Plot2 <- renderPlot({ rf_plot_list[[2]] })
          output$Build_RF_Plot3 <- renderPlot({ rf_plot_list[[3]] })
          output$Build_RF_Plot4 <- renderPlot({ rf_plot_list[[4]] })
          output$Build_RF_Plot5 <- renderPlot({ rf_plot_list[[5]] })
          
          rf_plot_grid<-grid.arrange(grobs = rf_plot_list, ncol = 1)
          output$downloadBuild_RF_Plot <- downloadHandler(filename = function(){paste("Random_Forest_Plot",'.png',sep='')},
                                                          content = function(file){
                                                            ggsave(file, plot= rf_plot_grid, width = 12, height = 20)})
        }
        if (input$RF_Build == TRUE){
          output$Build_RF_VIPPlot <- renderPlot({return(grid.arrange(grobs = rf_VIPplot_list, ncol = 2))})
          
          # Building the Random Forest VIP Plots individually for the Optimisation section:
          output$Build_RF_VIPPlot1 <- renderPlot({ rf_VIPplot_list[[1]] })
          output$Build_RF_VIPPlot2 <- renderPlot({ rf_VIPplot_list[[2]] })
          output$Build_RF_VIPPlot3 <- renderPlot({ rf_VIPplot_list[[3]] })
          output$Build_RF_VIPPlot4 <- renderPlot({ rf_VIPplot_list[[4]] })
          output$Build_RF_VIPPlot5 <- renderPlot({ rf_VIPplot_list[[5]] })
          
          output$downloadBuild_RF_VIPPlot <- downloadHandler(filename = function(){paste("Random_Forest_VIP_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= rf_VIPplot_list[[input$Build_RF_VIPPlot_Num]], width = 12, height = 12)})
        }
        if (input$SVM_Build == TRUE){
          output$Build_SVM_Plot <- renderPlot({return(grid.arrange(grobs = svm_plot_list, ncol = 1))})
          
          # Creating individual SVM Plots for the Optimisation section:
          output$Build_SVM_Plot1 <- renderPlot({ svm_plot_list[[1]] })
          output$Build_SVM_Plot2 <- renderPlot({ svm_plot_list[[2]] })
          output$Build_SVM_Plot3 <- renderPlot({ svm_plot_list[[3]] })
          output$Build_SVM_Plot4 <- renderPlot({ svm_plot_list[[4]] })
          output$Build_SVM_Plot5 <- renderPlot({ svm_plot_list[[5]] })
          
          svm_plot_grid<-grid.arrange(grobs = svm_plot_list, ncol = 1)
          output$downloadBuild_SVM_Plot <- downloadHandler(filename = function(){paste("SVM_Ensemble_Plot",'.png',sep='')},
                                                           content = function(file){
                                                             ggsave(file,plot= svm_plot_grid, width = 12, height = 20)})
        }
        if (input$LA_Build == TRUE){
          output$Build_LA_Plot <- renderPlot({return(grid.arrange(grobs = net_plot_list, ncol = 1))})
          
          # Building the individual net linear regression models for the Optimisation section:
          output$Build_LA_Plot1 <- renderPlot({ net_plot_list[[1]] })
          output$Build_LA_Plot2 <- renderPlot({ net_plot_list[[2]] })
          output$Build_LA_Plot3 <- renderPlot({ net_plot_list[[3]] })
          output$Build_LA_Plot4 <- renderPlot({ net_plot_list[[4]] })
          output$Build_LA_Plot5 <- renderPlot({ net_plot_list[[5]] })
          
          net_plot_grid<-grid.arrange(grobs = net_plot_list, ncol = 1)
          output$downloadBuild_LA_Plot <- downloadHandler(filename = function(){paste("LASSO_Analysis_Plot",'.png',sep='')},
                                                          content = function(file){
                                                            ggsave(file,plot= net_plot_grid, width = 12, height = 20)})
        }
        if (input$XGB_Build == TRUE){
          output$Build_XGB_Plot1 <- renderPlotly({return(`nround_ALL LM`)})
          output$Build_XGB_VIPPlot <- renderPlot({return(grid.arrange(grobs = XGB_VIPplot_list, ncol = 2))})
          
          
          # Building all Vip plots from the list:
          output$Buiild_XGB_VIPPlot1 <- renderPlotly({ XGB_VIPplot_list[[1]] })
          output$Buiild_XGB_VIPPlot2 <- renderPlotly({ XGB_VIPplot_list[[2]] })
          output$Buiild_XGB_VIPPlot3 <- renderPlotly({ XGB_VIPplot_list[[3]] })
          output$Buiild_XGB_VIPPlot4 <- renderPlotly({ XGB_VIPplot_list[[4]] })
          output$Buiild_XGB_VIPPlot5 <- renderPlotly({ XGB_VIPplot_list[[5]] })
          
          output$downloadBuild_XGB_VIPPlot <- downloadHandler(filename = function(){paste("XGBoost_VIP_Plot",'.png',sep='')},
                                                              content = function(file){
                                                                ggsave(file,plot= XGB_VIPplot_list[[input$Build_XGB_VIPPlot_Num]], width = 12, height = 12)})
          if (exists('DHA') == TRUE){
            output$Build_XGB_Plot2 <- renderPlotly({return(nround_DHA)})
          }
          if (exists('n3DPA') == TRUE) {
            output$Build_XGB_Plot3 <- renderPlotly({return(nround_n3DPA)})
          }
          if (exists('EPA') == TRUE) {
            output$Build_XGB_Plot4 <- renderPlotly({return(nround_EPA)})
          }
          if (exists('AA') == TRUE) {
            output$Build_XGB_Plot5 <- renderPlotly({return(nround_AA)})}
          XGB_Model_Table<-XGB_Model_Table[-1,]
          XGB_Model_Table[,-1] <-round(XGB_Model_Table[,-1],2)
          output$XGB_Build_Models<-renderDataTable(XGB_Model_Table, rownames = FALSE,
                                                   options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
          
          output$downloadBuild_XGB_Table <- downloadHandler(filename = function(){"XGB_Model_Table.csv"}, 
                                                            content = function(fname){
                                                              write.csv(XGB_Model_Table, fname)})
        }
        
        
        output$Build_ML_Table <- renderDataTable({return(final_table)})
        output$downloadBuild_Acc_ML_Table <- downloadHandler(filename = function(){"Accuracy_Table.csv"}, 
                                                             content = function(fname){
                                                               write.csv(final_table, fname)})
        
        
        
      })
      output$error_ML_Opt <- renderUI(NULL)
    },
    
    error = function(e) {
      # General error handling
      output$error_ML_Opt <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          This section uses the same file format as the first Machine Learning section.
          For more information, please look at the example file available through the 
          download example file button located in the ML Model page.", e$message
        )
      })
    }
    
    )
  })
  
  
  observe({
    tryCatch({
      req(input$LM_ML_File2_5)
      
      lm_profile = read.table(input$LM_ML_File2_5$datapath, 
                              sep=input$sep_LM_ML2_5,
                              header = TRUE,
                              row.names = 1,
                              stringsAsFactors = FALSE)
      
      Met_Rename<-c()
      #remove spaces from column names 
      colnames(lm_profile)<-gsub("\\.","", colnames(lm_profile)) 
      
      for (x in 1:ncol(lm_profile)){
        if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
        }else{
          Met_Rename[x] <-colnames(lm_profile)[x]
        }
      }
      colnames(lm_profile)<-Met_Rename
      
      #get the column names and update 
      Col_NamesML<-unique(as.character(unname(unlist(colnames(lm_profile)))))
      updateSelectInput(session, "Group_ML2_5", choices = Col_NamesML, selected = Col_NamesML[[1]])
      updateSelectInput(session, "MetName_ML", choices = Col_NamesML, selected = Col_NamesML[[1]])
      
      output$error_ML_Build <- renderUI(NULL)
      
},
    
    error = function(e) {
      # General error handling
      output$error_ML_Build <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          This section uses the same file format as the first Machine Learning section.
          For more information, please look at the example file available through the 
          download example file button located in the ML Model page.", e$message
        )
      })
    })
    
    
  })
  
  observeEvent(input$Build_ML2_5, {
    
    tryCatch({
      
      # Call the lipid mediator profiling file variant:
      inFile <- input$LM_ML_File2_5
      
      print(inFile)
      
      if(is.null(inFile))
        
        return(NULL)
      
      
      set.seed(415) # To get same results even with the random part.
      options(digits = 3) # To get only part of the decimals.
      
      # Open the file:
      
      
      lm_profile = read.table(inFile$datapath,
                              sep=input$sep_LM_ML2_5,
                              header = TRUE,
                              row.names = 1,
                              stringsAsFactors = FALSE)
      
      Met_Rename<-c()
      #remove spaces from column names 
      colnames(lm_profile)<-gsub("\\.","", colnames(lm_profile)) 
      
      for (x in 1:ncol(lm_profile)){
        if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
        }else{
          Met_Rename[x] <-colnames(lm_profile)[x]
        }
      }
      colnames(lm_profile)<-Met_Rename
      
      
      
      #remove fatty acid info row 
      lm_profile<-lm_profile[-1,]
      
      #obtain the responses info column
      ResponseName<-as.character(input$Group_ML2_5)
      Responses<-lm_profile[,ResponseName]
      
      #obtain columns of interest 
      lm_profileMet<-lm_profile[,(names(lm_profile) %in% input$MetName_ML)]
      
      # Creates the data frame that is going to collect all the info regarding the models:
      
      final_table <- data.frame(machine_learning = "del",
                                Accuracy = 1,
                                Sensitivity = 1,
                                Specificity = 1,
                                TP = 1,
                                FP = 1,
                                TN = 1,
                                FN = 1,
                                stringsAsFactors = FALSE)
      
      #progress bar
      #withProgress(message = 'Building Models', style = style, value = 0.1, {
      # Sys.sleep(0.25)})
      
      
      withProgress(message = 'Building Models:', detail = "This may take a while.", value = 0, {
        
        # Save all the values as numeric:
        lm_profile_number <- sapply(lm_profileMet, function(x) as.numeric(x))
        row.names(lm_profile_number) <- row.names(lm_profileMet)
        
        # Scale the data because is better when you are running Machine Learning models:
        lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))
        
        # If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors
        # replace the NA for zeros.
        lm_profiles_scale[is.na(lm_profiles_scale)] <- 0
        #--------------------
        ### change 0 to 1/5 of lowest value
        lm_profile_zero<-lm_profiles_scale[-1,]
        
        #get vector the same length as number of columns
        cols<- 1:(ncol(lm_profile_zero))
        
        #replace zeros for each column to 1/5 the smallest value for each column
        lm_profile_zero[cols] <- lapply(lm_profile_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
        
        #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
        lm_profile_zero<-as.matrix(lm_profile_zero)
        lm_profile_zero[is.infinite(lm_profile_zero)] <- NA
        min_index <- which.min(lm_profile_zero)
        zero_replace <- (lm_profile_zero[min_index]/5)
        lm_profile_zero <- as.data.frame(lm_profile_zero)
        lm_profile_zero[is.na(lm_profile_zero)] <- zero_replace
        
        #add first row to dataframe
        lm_profiles_scale<-rbind(lm_profiles_scale[1,] ,lm_profile_zero)
        
        #remove columns that have all the same value 
        #lm_profiles_scale<-Filter(var, lm_profiles_scale[-1,])
        
        #---------------------
        # Add the classification variable to the data frame (Responder and non responder):
        
        lm_profiles_scale$responses <- factor(Responses)
        
        # Make sure that column names do not represent a problem to randomForest making them a valid name to R.
        names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))
        
        #---> MACHINE LEARNING (randomForest R):
        if (input$RF_BuildMet == TRUE){
          
          # In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model.
          oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.
          for (mtry in 1:(ncol(lm_profiles_scale) - 1)) {
            rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry,
                                                  importance = TRUE, ntree = 10000)
            oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
          }
          # Define the best mtry according to the best prediction value.
          final_mtry <- which.max(oob_error)
          
          # Run the model again with the right mtry value.
          RF_Model <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry,
                                   importance = TRUE, ntree = 10000)
          # Get the confusion matrix of the model, sensitivity and specificity:
          RF_ConfMat <- as.data.frame(RF_Model$confusion)
          RF_ConfMat[is.na(RF_ConfMat)] <- 0
          
          # Calculates sensitivity, specificity and AUC.
          Sensi_RF <- RF_ConfMat[2, 2]/(RF_ConfMat[2, 2] + RF_ConfMat[2, 1])
          Speci_RF <- RF_ConfMat[1, 1]/(RF_ConfMat[1, 1] + RF_ConfMat[1, 2])
          
          # Final table for random forest:
          RF_Table <- data.frame(machine_learning = "RF",
                                 Accuracy = 100 - ((RF_Model$err.rate[10000])*100),
                                 Sensitivity = Sensi_RF,
                                 Specificity = Speci_RF,
                                 TP = (1 - RF_ConfMat[2, 3])*100,
                                 FP = RF_ConfMat[1, 3]*100,
                                 TN = (1 - RF_ConfMat[1, 3])*100,
                                 FN = RF_ConfMat[2, 3]*100,
                                 stringsAsFactors = FALSE)
          
          final_table <- rbind(final_table, RF_Table)
          # Number of trees plot:
          tree_table <- data.frame(Ensemble = c(1:10000),
                                   OBB = RF_Model$err.rate[, 1],
                                   AvgAcc = 100 - ((RF_Model$err.rate[, 1])*100))
          
          RF_Plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) +
            geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
            ggtitle("Random Forest nTrees Plot") +
            scale_x_continuous(name = "Trees") +
            scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
            theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 25, colour = "black"),
                  axis.title.x = element_text(size = 25, colour = "black"),
                  axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep.
                  legend.position = ("none"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(0.4, "cm"))
          #get importance dataframe
          rf_Importance<-as.data.frame(importance(RF_Model, type = 1))
          rf_Importance$lipid_mediators <- rownames(rf_Importance)
          
          Fatty_Acid <- 1:length(rf_Importance$lipid_mediators)
          cbind(rf_Importance, Fatty_Acid)
          
          #fatty acid points
          for (i in 1:length(rf_Importance$lipid_mediators)){
            rf_Importance$Fatty_Acid[i]<-SPMs_FA[rf_Importance$lipid_mediators[i]]
          }
          #change the names to be in the proper format
          rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)] <- paste(gsub("^X", "'", rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)]), "'", sep = "")
          
          rf_Importance$lipid_mediators <- gsub("A4'", "A'[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
          rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
          
          rf_Importance$lipid_mediators <- gsub("TXB2", "TXB[2]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("PGD2", "PGD[2]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("PGE2", "PGE[2]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("LXB4", "LXB[4]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("LXA4", "LXA[4]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("LTB4", "LTB[4]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("LTC4", "LTC[4]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("LTE4", "LTE[4]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("LTD4", "LTD[4]", rf_Importance$lipid_mediators)
          
          
          rf_Importance$lipid_mediators <- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required
          rf_Importance$lipid_mediators <- gsub("PGF2a", "PGF[2~a]", rf_Importance$lipid_mediators)
          rf_Importance$lipid_mediators <- gsub("n-3-dpa", "[n-3]~DPA", rf_Importance$lipid_mediators)
          
          #rf_Importance$lipid_mediators <- gsub("4$", "[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
          #rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
          #rf_Importance$lipid_mediators<- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required
          #rf_Importance$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", rf_Importance$lipid_mediators) # Transform n-3 DPA as a subscript:
          
          # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
          # (this with the purpose of get the LM in decreasing order in the figure):
          rf_Importance <- rf_Importance[order(rf_Importance$MeanDecreaseAccuracy), ]
          rf_Importance$lipid_mediators <- factor( rf_Importance$lipid_mediators,
                                                   levels = rf_Importance$lipid_mediators[order(rf_Importance$MeanDecreaseAccuracy)])
          
          lipids <- as.character(rf_Importance$lipid_mediators)
          
          #add colors to SPM names
          lm_classes <- as.factor(rf_Importance$Fatty_Acid)
          names(lm_classes) <- rf_Importance$lipid_mediators
          
          dha_index<-which((lm_classes=="DHA")== TRUE)
          n_three_index<-which((lm_classes=="n3DPA")== TRUE)
          epa_index<-which((lm_classes=="EPA")== TRUE)
          aa_index<-which((lm_classes=="AA")== TRUE)
          
          lm_colors <- NULL
          lm_colors[dha_index] <- "blue"
          lm_colors[n_three_index] <- "brown"
          lm_colors[epa_index] <- "darkgoldenrod1"
          lm_colors[aa_index] <- "darkslategray"
          
          #rf_Importance1<- subset(rf_Importance, MeanDecreaseAccuracy >1)
          RF_VIP_Plot <- ggplot(data = rf_Importance, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = Fatty_Acid)) + geom_point(size = 3) +
            scale_y_continuous(name = "Mean Decrease Accuracy") +
            labs(x = "Lipid Mediators", title = "Radnom Forest VIP Plot" )+
            scale_x_discrete(labels = parse(text = lipids)) +
            scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                 "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
            coord_flip() +
            theme(axis.title = element_text(size = 20),
                  axis.text.x  = element_text(size = 15, hjust = 0.5),
                  axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                  legend.position = "top",
                  aspect.ratio = 2/1,
                  legend.title = element_text(size = 20),
                  legend.text  = element_text(size = 15),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
          
          #outputs for random forest 
          output$BuildMet_RF_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("Optimal paramters for Random Forest", align = "center") })
          output$BuildMet_RF_VIPPlot_Title <- renderUI({req(input$Build_ML2_5); h2("Importance of Variance Plot", align = "center") })
          
          output$BuildMet_RF_Plot <- renderPlot(RF_Plot)
          output$BuildMet_RF_VIPPlot <- renderPlot(RF_VIP_Plot)
          
          
          
          output$downloadBuildMet_RF_Plot <- downloadHandler(filename = function(){paste("RF_ntrees_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= RF_Plot, width = 12, height = 8)})
          output$downloadBuildMet_RF_VIPPlot <- downloadHandler(filename = function(){paste("RF_VIP_Plot",'.png',sep='')},
                                                                content = function(file){
                                                                  ggsave(file,plot= RF_VIP_Plot, width = 12, height = 8)})
          output$downloadBuildMet_RF_Mod <- downloadHandler(filename = function(){"RandomForest_Model.rds"},
                                                            content = function(fname){
                                                              saveRDS(RF_Model, compress = TRUE, fname)})
          output$downloadBuildMet_RF_Mod1 <- downloadHandler(filename = function(){"RandomForest_Model.rds"},
                                                             content = function(fname){
                                                               saveRDS(RF_Model, fname)})
        }
        
        #---> MACHINE LEARNING (EXTREME GRADIENT BOOSTING)
        if (input$XGB_BuildMet == TRUE){
          df<-lm_profiles_scale
          
          #need to convert responses to integer of 0 and 1
          df$responses<-as.factor(df$responses)
          responses<-df$responses
          label<-as.integer(df$responses) -1
          
          #remove labels from dataframe
          df$responses<-NULL
          df<-sapply(df, function(x) as.numeric(x))
          
          #split it into training and testing dataset with 75/20
          n = nrow(df) #get number of rows for dataframe
          train.index = sample(n,floor(0.7*n)) #randomize rows to get dataset
          train.data = as.matrix(df[train.index,])
          train.label = label[train.index]
          test.data = as.matrix(df[-train.index,])
          test.label = label[-train.index]
          
          #transform dataset to xgb matrix
          xgb.train = xgb.DMatrix(data=train.data,label=train.label)
          xgb.test = xgb.DMatrix(data=test.data, label=test.label)
          
          # #create grid with all the options for xgboost
          # searchGridSubCol <- expand.grid(subsample = c(0.5, 0.75, 1),
          #                                 colsample_bytree = c(0.5, 0.75, 1),
          #                                 max_depth = c(2, 4, 6, 8, 10),
          #                                 gamma_val = c(0, 1, 5, 10),
          #                                 eta = c(0.01, 0.1, 0.2, 0.3),
          #                                 min_child = seq(1))
          
          searchGridSubCol <- expand.grid(subsample = c(0.6, 0.8, 1), 
                                          colsample_bytree = c(0.6, 0.8, 1),
                                          max_depth = c(3, 4, 5),
                                          gamma_val = c(0.5, 1, 1.5, 2, 5), 
                                          eta = c(0.1, 0.2, 0.3),
                                          min_child = c(1, 5, 10)
          )
          
          #run through each combination of the grid
          system.time(
            ErrorsHyperparameters <- apply(searchGridSubCol, 1, function(parameterList){
              #Extract Parameters to test
              currentSubsampleRate <- parameterList[["subsample"]]
              currentColsampleRate <- parameterList[["colsample_bytree"]]
              currentDepth <- parameterList[["max_depth"]]
              currentEta <- parameterList[["eta"]]
              currentGamma <- parameterList[["gamma_val"]]
              currentMinChild <- parameterList[["min_child"]]
              
              #run selected parameter through the model
              xgboostModelCV <- xgb.cv(data =  xgb.train, nrounds = 10000, nfold = 5, showsd = TRUE,
                                       metrics = "error", verbose = TRUE, "eval_metric" = "error",
                                       "objective" = "binary:logistic", "max.depth" = currentDepth, "eta" = currentEta,
                                       "subsample" = currentSubsampleRate, "colsample_bytree" = currentColsampleRate, 
                                       print_every_n = 10, booster = "gbtree",
                                       early_stopping_rounds = 10, "gamma" = currentGamma, "min_child_weight" = currentMinChild)
              
              #have error evaluation score as dataframe
              xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
              test_error <- tail(xvalidationScores$test_error_mean, 1)
              train_error <- tail(xvalidationScores$train_error_mean,1)
              output <- return(c(test_error, train_error, currentSubsampleRate, currentColsampleRate, currentDepth, currentEta, currentGamma, currentMinChild))})) #add to final table called output
          
          output_Error <- as.data.frame(t(ErrorsHyperparameters))
          varnames <- c("TestError", "TrainError", "SubSampRate", "ColSampRate", "Depth", "eta", "gamma", "currentMinChild")
          names(output_Error) <- varnames #change the variable names
          top_output<-output_Error[order(output_Error$TestError, output_Error$TrainError), ] #order it by the lowest error rate
          
          #empty table to put the top 10 model information into
          Final_Acc <- data.frame(machine_learning = "test",
                                  Accuracy = 0,
                                  Specificity = 0,
                                  Sensitivity = 0,
                                  TP = 0,
                                  FP = 0,
                                  TN = 0,
                                  FN = 0,
                                  nrounds = 0,
                                  eta= 0,
                                  max_depth= 0,
                                  gamma= 0,
                                  min_child_weight=0,
                                  subsample=0,
                                  colsample_bytree= 0)
          top<-1:50
          #loop to build models for the parameters with the lowest error rate
          for (i in top){
            
            #parameter list with the lowest output
            params <- list(booster = "gbtree", objective = "binary:logistic", eta= top_output[i,6], gamma= top_output[i,7],
                           max_depth= top_output[i,5], min_child_weight=top_output[i,8], subsample=top_output[i,3],
                           colsample_bytree= top_output[i,4])
            
            #Build model based on the training data
            xgb<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
            
            #run model on test data = see how good it is at predicting
            pred<-predict(xgb, test.data)
            
            #obtain confusion matrix
            pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1
            pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
            pred<-as.numeric(pred) #makes it numeric
            test.label_num <-as.numeric(test.label) #makes it numeric
            Responses_List<-unique(Responses)
            pred_y <- cut(pred, 2, label = c(Responses_List[1], Responses_List[2])) #converts to factor of responder or non responder
            test_y <-cut(test.label_num, 2, label = c(Responses_List[1], Responses_List[2])) #convers to factor of responder or non responder
            ConMatrix <- confusionMatrix(reference = test_y, data = pred_y) #build confusion matrix
            TP = as.numeric(ConMatrix$table[2,2])
            FP = as.numeric(ConMatrix$table[2,1])
            TN = as.numeric(ConMatrix$table[1,1])
            FN = as.numeric(ConMatrix$table[1,2])
            
            #convert build table with accuracy for each model
            Table <- data.frame(machine_learning = "XGBoost",
                                Accuracy = as.numeric(ConMatrix$overall[1])*100,
                                Specificity = as.numeric(ConMatrix$byClass[2]),
                                Sensitivity = as.numeric(ConMatrix$byClass[1]),
                                TP = (TP/(TP+FN)*100),
                                FP = (FP/(FP+TN)*100),
                                TN = (TN/(FP+TN)*100),
                                FN = (FN/(TP+FN)*100),
                                nrounds = 10000,
                                eta= top_output[i,6],
                                max_depth= top_output[i,5],
                                gamma= top_output[i,7],
                                min_child_weight=top_output[i,8],
                                subsample=top_output[i,3],
                                colsample_bytree= top_output[i,4])
            Final_Acc<-rbind(Final_Acc, Table)
          }
          
          #remove the first row of random values
          Final_Acc<-Final_Acc[-1,]
          
          #order the top 10 models based on highest accuracy, sensitivity, and specificity
          Top_Acc<-Final_Acc[with(Final_Acc, order(-Accuracy, -Sensitivity, -Specificity)), ]
          
          
          #preparing nrounds plots for top 5 models
          iterations<-1:10000
          top_table<- data.frame(iterations)
          
          params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[1,10], gamma= Top_Acc[1,12],
                         max_depth= Top_Acc[1,11], min_child_weight= Top_Acc[1,13], subsample= Top_Acc[1,14])
          
          #run the cv to get the error rate for each nround
          xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = 10000, nfold = 5, showsd = T,
                           stratified = T, early_stop_round = 20, maximize = F, metrics = "error", verbose = FALSE)
          # 
          Table1<-data.frame( test = as.numeric(xgbcv$evaluation_log$test_error_mean))
          colnames(Table1)[1] <-"XGBoost_Model"
          top_table<-cbind(top_table, Table1) #add to table
          
          #plot using plotly
          nround_plot<-plot_ly(data = top_table, x = ~iterations)
          nround_plot <-  nround_plot %>% add_trace(y = ~XGBoost_Model, name = colnames(top_table)[2], mode = 'lines')
          nround_plot <- nround_plot %>% layout(title =  "XGBoost Plot", xaxis = list(title = 'Number of Rounds'),
                                                yaxis = list(title = 'Test Error Rate'))
          #makes download quality of plot better
          nround_plot <- nround_plot%>% config(toImageButtonOptions = list(format = "jpeg", width = 1500, height = 750))
          
          
          
          #importance plot
          #build model for the top model
          params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[1,10], gamma= Top_Acc[1,12],
                         max_depth= Top_Acc[1,11], min_child_weight= Top_Acc[1,13], subsample= Top_Acc[1,14], colsample_bytree= Top_Acc[1,15])
          xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
          
          
          XGB_Importance = xgb.importance(model = xgb_Mod1) #important matrix
          
          Fatty_Acid <- 1:length(XGB_Importance$Feature)
          cbind(XGB_Importance, Fatty_Acid)
          
          #fatty acid points
          for (i in 1:length(XGB_Importance$Feature)){
            XGB_Importance$Fatty_Acid[i]<-SPMs_FA[XGB_Importance$Feature[i]]
          }
          
          # #change the names to be in the proper format
          XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)] <- paste(gsub("^X", "'", XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)]), "'", sep = "")
          
          XGB_Importance$Feature <- gsub("A4'", "A'[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 5
          XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature) #B4 converted to subscript 4
          
          XGB_Importance$Feature <- gsub("TXB2", "TXB[2]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("PGD2", "PGD[2]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("PGE2", "PGE[2]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("LXB4", "LXB[4]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("LXA4", "LXA[4]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("LTB4", "LTB[4]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("LTC4", "LTC[4]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("LTE4", "LTE[4]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("LTED", "LTD[4]", XGB_Importance$Feature)
          
          XGB_Importance$Feature <- gsub("\\.", "-", XGB_Importance$Feature) # Replace "." with "-" when required
          XGB_Importance$Feature <- gsub("PGF2a", "PGF[2~a]", XGB_Importance$Feature)
          XGB_Importance$Feature <- gsub("n-3-dpa", "[n-3]~DPA", XGB_Importance$Feature)
          
          #XGB_Importance$Feature <- gsub("4$", "[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 4
          #XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature ) #B4 converted to subscript 4
          #XGB_Importance$Feature  <- gsub("\\.", "-", XGB_Importance$Feature ) # Replace "." with "-" when required
          #XGB_Importance$Feature <- gsub("n.3.DPA", "[n-3~DPA]", XGB_Importance$Feature ) # Transform n-3 DPA as a subscript:
          
          
          # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
          # (this with the purpose of get the LM in decreasing order in the figure):
          
          XGB_Importance <- XGB_Importance[order(XGB_Importance$Gain), ]
          
          XGB_Importance$Feature<- factor(XGB_Importance$Feature,
                                          levels = XGB_Importance$Feature[
                                            order(XGB_Importance$Gain)])
          lipids <- as.character(XGB_Importance$Feature )
          
          #add colors to x axis labels
          lm_classes <- as.factor(XGB_Importance$Fatty_Acid)
          names(lm_classes) <- XGB_Importance$Feature
          
          dha_index<-which((lm_classes=="DHA")== TRUE)
          n_three_index<-which((lm_classes=="n3DPA")== TRUE)
          epa_index<-which((lm_classes=="EPA")== TRUE)
          aa_index<-which((lm_classes=="AA")== TRUE)
          
          lm_colors <- NULL
          lm_colors[dha_index] <- "blue"
          lm_colors[n_three_index] <- "brown"
          lm_colors[epa_index] <- "darkgoldenrod1"
          lm_colors[aa_index] <- "darkslategray"
          
          XGB_VIP_Plot <- ggplot(data = XGB_Importance, mapping = aes(x = Feature, y = Gain, color = Fatty_Acid)) + geom_point(size = 3) +
            scale_y_continuous(name = "Gain") +
            labs(x = "Features", title = "XGBoost VIP Plot") +
            scale_x_discrete(labels = parse(text = lipids)) +
            scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                 "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
            coord_flip() +
            theme(axis.title = element_text(size = 20),
                  axis.text.x  = element_text(size = 15, hjust = 0.5),
                  axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                  legend.position = "top",
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
          XGBoost_Table<-Top_Acc[1, 1:8]
          final_table<-rbind(final_table, XGBoost_Table)
          
          #outputs plots and downloads 
          output$BuildMet_XGB_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("XGBoost nrounds vs Test Error Rate Plot", align = "center") })
          output$BuildMet_XGB_VIPPlot_Title <- renderUI({req(input$Build_ML2_5); h2("Importance of Variance Plot", align = "center") })
          
          output$BuildMet_XGB_Plot <- renderPlotly(nround_plot)
          output$BuildMet_XGB_VIPPlot <- renderPlot(XGB_VIP_Plot)
          output$downloadBuildMet_XGB_VIPPlot <- downloadHandler(filename = function(){paste("XGBoost_VIP_Plot",'.png',sep='')},
                                                                 content = function(file){
                                                                   ggsave(file,plot= XGB_VIP_Plot, width = 12, height = 8)})
          output$downloadBuildMet_XGB_Mod <- downloadHandler(filename = function(){"XGBoost_Model.rds"},
                                                             content = function(fname){
                                                               saveRDS(xgb_Mod1, fname)})
          output$downloadBuildMet_XGB_Mod1 <- downloadHandler(filename = function(){"XGBoost_Model.rds"},
                                                              content = function(fname){
                                                                saveRDS(xgb_Mod1, fname)})
        }
        
        
        #---> MACHINE LEARNING (Classyfire R):
        
        # Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
        # and creates a novelty detection for the creation of the model.
        
        # The idea is to create several models and see which one fits the best. The models will be based on the whole
        # lipid profiles and the different groups based on substrates.
        
        # "cfBuild" to create the SVM:
        # Clasyfire requieres matrix:
        if (input$SVM_BuildMet == TRUE){
          #convert dataframe to matrix
          lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])
          
          #building SVM model
          support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses,
                                              bootNum = 70, ensNum = 70, cpus = 4)
          #obtaining confusion matrix
          conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale))
          
          # SVM table:
          Support_vector_table <- data.frame(machine_learning = "SVM",
                                             Accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                             Sensitivity = conf_matrix[4, 3]/100,
                                             Specificity = conf_matrix[1, 3]/100,
                                             TP = conf_matrix[4, 3],
                                             FP = conf_matrix[2, 3],
                                             TN = conf_matrix[1, 3],
                                             FN = conf_matrix[3, 3],
                                             stringsAsFactors = FALSE)
          final_table <- rbind(final_table, Support_vector_table)
          
          # Ensemble Plot:
          ensAcc   <- getAcc(support_lmprofiles_scale)$Test
          meanVal  <- ensAcc[1]
          for (i in 2:length(ensAcc)) {
            meanVal <- c(meanVal, mean(ensAcc[1:i]))
          }
          ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc),
                                      AvgAcc = meanVal)
          SVM_Plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) +
            geom_point(aes(colour = AvgAcc), size = 5) +
            geom_line(linetype = "dotted", size = 1) +
            ggtitle("SVM Ensemble Plot") +
            scale_x_continuous(name = "Ensemble interaction") +
            scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
            theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 25, colour = "black"),
                  axis.title.x = element_text(size = 25, colour = "black"),
                  axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep.
                  legend.position = ("none"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(0.4, "cm"))
          
          output$BuildMet_SVM_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("Optimal Parameters SVM", align = "center") })
          output$BuildMet_SVM_Plot <- renderPlot(SVM_Plot)
          output$downloadBuildMet_SVM_Plot <- downloadHandler(filename = function(){paste("SVM_Ensemble_Plot",'.png',sep='')},
                                                              content = function(file){
                                                                ggsave(file,plot= SVM_plot, width = 12, height = 8)})
          output$downloadBuildMet_SVM_Mod <- downloadHandler(filename = function(){"SVM_Model.rds"},
                                                             content = function(fname){
                                                               saveRDS(support_lmprofiles_scale, fname)})
          
        }
        
        #---> ELASTIC NET REGRESSION (caret R):
        if (input$LA_BuildMet == TRUE){
          # Get the explanatory variables as a matrix:
          explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
          
          # LASSO Analysis:
          model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet",
                             trControl = trainControl("boot", number = 70))
          # Get the confusion matrix of the model:
          conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
          
          # Final Elastic net model table:
          net_table <- data.frame(machine_learning = "ElasticNet",
                                  Accuracy = max(model_net$results$Accuracy)*100,
                                  Sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                                  Specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                                  TP = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                  FP = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                  TN = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                  FN = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                  stringsAsFactors = FALSE)
          final_table <- rbind(final_table, net_table)
          
          # Parameter tuning figure:
          scaleFUN <- function(x) sprintf("%.2f", x)
          
          Elastic_Net_Plot <- ggplot(model_net, highlight = TRUE) +
            scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
            scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
            ggtitle("Elastic Net Parameters Plot") +
            scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
            scale_shape_manual(values=c(16, 16, 16)) +
            labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
            guides(shape = FALSE) +
            theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 25, colour = "black"),
                  axis.title.x = element_text(size = 25, colour = "black"),
                  axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep.
                  legend.title = element_text(size = 25),
                  legend.text  = element_text(size = 20),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(0.4, "cm"))
          
          output$BuildMet_LA_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("Optimal Parameters Elastic Net Regression", align = "center") })
          output$BuildMet_LA_Plot <- renderPlot(Elastic_Net_Plot)
          output$downloadBuildMet_LA_Plot <- downloadHandler(filename = function(){paste("ElasticNet_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= Elastic_Net_Plot, width = 12, height = 8)})
          output$downloadBuildMet_LA_Mod <- downloadHandler(filename = function(){"ElasticNet_Model.rds"},
                                                            content = function(fname){
                                                              saveRDS(model_net, fname)})
          
        }
        
        #---> BAYESIAN MODEL (Caret R):
        
        if (input$BC_BuildMet == TRUE){
          explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
          
          bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm",
                            trControl = trainControl("boot", number = 70))
          conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)
          
          # Final Bayesian model table:
          bay_table <- data.frame(machine_learning = "Bayes GLM",
                                  Accuracy = max(bayesian$results$Accuracy)*100,
                                  Sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                                  Specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                                  TP = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                  FP = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                  TN = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                  FN = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                  stringsAsFactors = FALSE)
          final_table <- rbind(final_table, bay_table)
          
          output$downloadBuildMet_BC_Mod <- downloadHandler(filename = function(){"BayesGLM_Model.rds"},
                                                            content = function(fname){
                                                              saveRDS(bayesian, fname)})
        }
        
        final_table <- final_table[-1, ]
        final_table$Accuracy <- round(final_table$Accuracy, 0)
        final_table$Sensitivity <- round(final_table$Sensitivity, 2)
        final_table$Specificity <- round(final_table$Specificity, 2)
        final_table$TP <- round(final_table$TP, 0)
        final_table$FP <- round(final_table$FP, 0)
        final_table$TN <- round(final_table$TN, 0)
        final_table$FN <- round(final_table$FN, 0)
        
        accuracy <- ggplot(data = final_table, aes(x = machine_learning, y = Accuracy, fill = machine_learning)) +
          geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black", size = 2) +
          scale_fill_manual(values=c("dodgerblue2","firebrick2",'goldenrod1', 'lightslategray', "purple")) + 
          geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
                    aes(label = paste(Accuracy, "%", sep = "")), size = 8) +
          scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20),
                             expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
          coord_cartesian(ylim = c(1, 120)) +
          theme(axis.title = element_text(size = 20),
                axis.title.x = element_blank(),
                axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
                axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"), # Put color to the labels
                axis.line = element_line(colour = 'black', size = 1.0), # Color and thickness of axis
                axis.ticks = element_line(colour = "black", size = 1.0), # Color and thickness of every axis sep. 
                panel.background = element_rect(fill = "white"),
                legend.title = element_blank(),
                legend.position = "top",
                legend.key.size = unit(1.3, "cm"), 
                legend.text  = element_text(size = 25),
                legend.spacing.x = unit(1, "cm"),
                axis.ticks.length = unit(0.4, "cm"))
        
        
        # Final table cute:
        
        colnames(final_table)[1] <- "Machine Learning Methodology"
        
        output$BuildMet_AccPlot_Title <- renderUI({req(input$Build_ML2_5); h2("% Accuracy Score Figure for different ML models", align = "center") })
        output$BuildMet_AccPlot <- renderPlot(accuracy)
        output$downloadBuildMet_AccPlot <- downloadHandler(filename = function(){paste("Accuracy_Plot",'.png',sep='')},
                                                           content = function(file){
                                                             ggsave(file, plot= accuracy, width = 12, height = 8)})
        
        output$BuildMet_ML_Table <- renderDataTable({return(final_table)})
        output$downloadBuildMet_Acc_ML_Table <- downloadHandler(filename = function(){"Accuracy_Table.csv"},
                                                                content = function(fname){
                                                                  write.csv(final_table, fname)})
        
      })
      output$error_ML_Build <- renderUI(NULL)
    },
    error = function(e) {
      # General error handling
      output$error_ML_Build <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          This section uses the same file format as the first Machine Learning section.
          For more information, please look at the example file available through the 
          download example file button located in the ML Model page.", e$message
        )
      })
    })  
    
    
  })
  
  observe({
    
      req(input$LM_ML_File3)
    tryCatch({
      
      lm_profile = read.table(input$LM_ML_File3$datapath, 
                              sep=input$sep_LM_ML3,
                              header = TRUE,
                              row.names = 1,
                              stringsAsFactors = FALSE)
      
      Met_Rename<-c()
      #remove spaces from column names 
      colnames(lm_profile)<-gsub("\\.","", colnames(lm_profile)) 
      
      for (x in 1:ncol(lm_profile)){
        if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
        }else{
          Met_Rename[x] <-colnames(lm_profile)[x]
        }
      }
      colnames(lm_profile)<-Met_Rename
      
      #get the column names and update 
      Col_NamesML<-unique(as.character(unname(unlist(colnames(lm_profile)))))
      updateSelectInput(session, "Group_ML3", choices = Col_NamesML, selected = Col_NamesML[[1]])
      updateSelectInput(session, "MetName_ML3", choices = Col_NamesML, selected = Col_NamesML[[1]])
      
      output$error_ML_Run <- renderUI(NULL)
    },
    error = function(e) {
      # General error handling
      output$error_ML_Run <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          For more information, please look at the example file available through the 
          download example file button.", e$message
        )
      })
    }
    )
    
  })
  
  observeEvent(input$Run_ML,{
    
    tryCatch({
      
      
      # Call the lipid mediator profiling file variant:
      inFile <- input$LM_ML_File3
      set.seed(415) # To get same results even with the random part.
      options(digits = 3) # To get only part of the decimals.
      
      # Open the file:
      lm_profiles = read.csv(input$LM_ML_File3$datapath, 
                             sep=input$sep_LM_ML3,
                             header = TRUE,
                             row.names = 1,
                             stringsAsFactors = FALSE)
      Met_Rename<-c()
      #remove spaces from column names 
      colnames(lm_profiles)<-gsub("\\.","", colnames(lm_profiles)) 
      
      for (x in 1:ncol(lm_profiles)){
        if (colnames(lm_profiles)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(lm_profiles)[x]]]
        }else{
          Met_Rename[x] <-colnames(lm_profiles)[x]
        }
      }
      colnames(lm_profiles)<-Met_Rename
      
      
      #obtain the responses info column
      ResponseName<-as.character(input$Group_ML3)
      
      #Zero handing and scaling data 
      label <- lm_profiles[,ResponseName]
      if (input$Choose_MetName3 == TRUE){
        Val_Data_Raw<- lm_profiles[,(names(lm_profiles) %in% input$MetName_ML3)]
      } else{
        Val_Data_Raw<- lm_profiles[,!(names(lm_profiles) %in% ResponseName)]
      }
      
      Val_Data_Raw<-sapply(Val_Data_Raw, function(x) as.numeric(x)) #convert to numeric
      row.names(Val_Data_Raw) <- row.names(label)
      Val_Data_Scale<-as.data.frame(scale(Val_Data_Raw, center = FALSE, scale = TRUE)) #scale data
      Val_Data_Scale[is.na(Val_Data_Scale)] <- 0 #convert NA to 0
      
      #zero handling to make 0 to 1/5 lowest value
      
      Val_Data_zero<-Val_Data_Scale
      
      
      #get vector the same length as number of columns
      cols<- 1:(ncol(Val_Data_zero))
      
      
      #replace zeros for each column to 1/5 the smallest value for each column
      
      Val_Data_zero[cols] <- lapply(Val_Data_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
      
      #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
      Val_Data_zero<-as.matrix(Val_Data_zero)
      Val_Data_zero[is.infinite(Val_Data_zero)] <- NA
      min_index <- which.min(Val_Data_zero)
      zero_replace <- (Val_Data_zero[min_index]/5)
      Val_Data_zero <- as.data.frame(Val_Data_zero)
      Val_Data_zero[is.na(Val_Data_zero)] <- zero_replace
      
      Val_Data<-Val_Data_zero
      
      #remove columns with all the same values 
      #Val_Data<-Filter(var, Val_Data_zero)
      
      #make sure column names do not cause issue in random forest
      names(label) <- make.names(names(label))
      names(Val_Data) <- make.names(names(Val_Data))
      
      # #No data preprocessing 
      # label <- lm_profiles[,1]
      # Val_Data<- lm_profiles[,-1]
      # Val_Data<-sapply(Val_Data, function(x) as.numeric(x)) #convert to numeric
      
      #output$Run_ValData<-renderTable(Val_Data)
      ROC_List<- list()
      ROC_Names<-list()
      pos = 1
      
      #final ROC Table
      ROC_Table <- data.frame(Model = "del",
                              AUC = 1,
                              Accuracy = 1,
                              Sensitivity = 1,
                              Specificity = 1,
                              TP_per = 1,
                              FP_per = 1,
                              TN_per = 1,
                              FN_per = 1,
                              stringsAsFactors = FALSE)
      
      if (input$RF_Run == TRUE){
        RF_Model <- readRDS(input$RF_Mod_File$datapath)
        pred_RF<-as.data.frame(predict(RF_Model, Val_Data, type = 'prob'))
        RF_ROC<-pROC::roc(label,pred_RF$Non_Responder)
        
        RF_ROCPlot<- ggroc(RF_ROC, legacy.axes = TRUE) +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
          theme_bw()
        
        output$Run_RF_ROCPlot<-renderPlot({return(RF_ROCPlot)})
        
        output$downloadRun_RF_ROCPlot <- downloadHandler(filename = function(){paste("RandomForest_ROC_Plot",'.png',sep='')},
                                                         content = function(file){
                                                           ggsave(file,plot= RF_ROCPlot, width = 15, height = 10)})
        
        # ROC_List<-append(ROC_List, RF_ROC)
        # ROC_Names<-append(ROC_Names, "Random Forest")
        
        ROC_List[[pos]] <- RF_ROC
        ROC_Names[[pos]] <- "RF"
        pos = pos + 1
        
        #get table information---confusion matrix 
        # Get the confusion matrix of the model, sensitivity and specificity: 
        RF_ConfMat<- as.data.frame(RF_Model$confusion)
        RF_ConfMat[is.na(RF_ConfMat)] <- 0
        
        # Calculates sensitivity, specificity and AUC.
        RF_Sensi <- RF_ConfMat[2, 2]/(RF_ConfMat[2, 2] + RF_ConfMat[2, 1])
        RF_Speci <- RF_ConfMat[1, 1]/(RF_ConfMat[1, 1] + RF_ConfMat[1, 2])
        
        # Final table for random forest:
        RF_Table <- data.frame(Model = "Random Forest",
                               AUC = auc(RF_ROC),
                               Accuracy = 100 - ((RF_Model$err.rate[10000])*100),
                               Sensitivity = RF_Sensi,
                               Specificity = RF_Speci,
                               TP_per = (1 - RF_ConfMat[2, 3])*100,
                               FP_per = RF_ConfMat[1, 3]*100,
                               TN_per = (1 - RF_ConfMat[1, 3])*100,
                               FN_per = RF_ConfMat[2, 3]*100,
                               stringsAsFactors = FALSE)
        ROC_Table <- rbind(ROC_Table, RF_Table)
      }
      
      if (input$XGB_Run == TRUE){
        XGB_Model <- readRDS(input$XGB_Mod_File$datapath)
        pred_XGB<-as.data.frame(predict(XGB_Model, as.matrix(Val_Data),  type = 'prob'))
        XGB_ROC<-pROC::roc(label, pred_XGB[,1])
        
        XGB_ROCPlot<- ggroc(XGB_ROC, legacy.axes = TRUE) +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
          theme_bw()
        
        #XGB_ROCPlot<-plot(XGB_ROC)
        
        output$Run_XGB_ROCPlot<-renderPlot({return(XGB_ROCPlot)})
        
        output$downloadRun_XGB_ROCPlot <- downloadHandler(filename = function(){paste("XGBoost_ROC_Plot",'.png',sep='')},
                                                          content = function(file){
                                                            ggsave(file,plot= XGB_ROCPlot, width = 15, height = 10)})
        
        # ROC_List<-append(ROC_List, XGB_ROC)
        # ROC_Names<-append(ROC_Names, "XGBoost")
        ROC_List[[pos]] <- XGB_ROC
        ROC_Names[[pos]] <- "XGB"
        pos = pos + 1
        
        #get table information for Confusion matrix
        pred<-pred_XGB[,1]
        pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1 
        pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
        pred<-as.numeric(pred) #makes it numeric
        label_xgb<-as.factor(label)
        label_xgb<-as.integer(label_xgb) -1
        label_xgb_num<-as.numeric(label_xgb)
        label_list<-unique(label) #list of unique values 
        pred_y <- cut(pred, 2, label = c(label_list[1], label_list[2])) #converts to factor of responder or non responder 
        label_y <-cut(label_xgb_num, 2, label = c(label_list[1], label_list[2])) #convers to factor of responder or non responder 
        XGB_ConfMat <- confusionMatrix(reference = label_y, data = pred_y) #build confusion matrix
        
        TP = as.numeric(XGB_ConfMat$table[2,2])
        FP = as.numeric(XGB_ConfMat$table[2,1])
        TN = as.numeric(XGB_ConfMat$table[1,1])
        FN = as.numeric(XGB_ConfMat$table[1,2])
        
        #convert build table with accuracy for each model
        XGB_Table <- data.frame(Model = "Extreme Gradient Boosting",
                                AUC = auc(XGB_ROC),
                                Accuracy = as.numeric(XGB_ConfMat$overall[1])*100,
                                Specificity = as.numeric(XGB_ConfMat$byClass[2]),
                                Sensitivity = as.numeric(XGB_ConfMat$byClass[1]),
                                TP_per = (TP/(TP+FN)*100),
                                FP_per = (FP/(FP+TN)*100),
                                TN_per = (TN/(FP+TN)*100),
                                FN_per = (FN/(TP+FN)*100),
                                stringsAsFactors = FALSE)
        
        ROC_Table <- rbind(ROC_Table, XGB_Table)
      }
      
      if (input$SVM_Run == TRUE){
        SVM_Model <- readRDS(input$SVM_Mod_File$datapath)
        pred_SVM<-cfPredict(SVM_Model, Val_Data)
        names(pred_SVM)<-c("prediction", "Coef Score")
        label_list<-unique(label) #list of unique values (the group names )
        
        pred_SVM$Responder[pred_SVM$prediction == label_list[1]] <- pred_SVM$`Coef Score`[pred_SVM$prediction == label_list[1]]/100 
        pred_SVM$Responder[pred_SVM$prediction == label_list[2]] <- 1 - (pred_SVM$`Coef Score`[pred_SVM$prediction == label_list[2]]/100)
        pred_SVM$Non_Responder <- 1 - pred_SVM$Responder
        
        SVM_ROC<-pROC::roc(label,pred_SVM$Non_Responder)
        
        SVM_ROCPlot<- ggroc(SVM_ROC, legacy.axes = TRUE) +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
          theme_bw()
        
        output$Run_SVM_ROCPlot<-renderPlot({return(SVM_ROCPlot)})
        
        output$downloadRun_SVM_ROCPlot <- downloadHandler(filename = function(){paste("SVM_ROC_Plot",'.png',sep='')},
                                                          content = function(file){
                                                            ggsave(file,plot= SVM_ROCPlot, width = 15, height = 10)})
        
        
        # ROC_List<-append(ROC_List, SVM_ROC)
        # ROC_Names<-append(ROC_Names, "SVM")
        
        ROC_List[[pos]] <- SVM_ROC
        ROC_Names[[pos]] <- "SVM"
        pos = pos + 1
        
        #get  confusion matrix informaiton
        #obtaining confusion matrix 
        SVM_ConfMat <- as.data.frame(getConfMatr(SVM_Model))
        
        # SVM table:
        SVM_Table <- data.frame(Model = "Support Vector Machine",
                                AUC = auc(SVM_ROC),
                                Accuracy = getAvgAcc(SVM_Model)$Test,
                                Sensitivity = SVM_ConfMat[4, 3]/100,
                                Specificity = SVM_ConfMat[1, 3]/100,
                                TP_per = SVM_ConfMat[4, 3],
                                FP_per = SVM_ConfMat[2, 3],
                                TN_per = SVM_ConfMat[1, 3],
                                FN_per = SVM_ConfMat[3, 3],
                                stringsAsFactors = FALSE)
        ROC_Table <- rbind(ROC_Table, SVM_Table)
      }
      
      
      if (input$LA_Run == TRUE){
        LA_Model <- readRDS(input$LA_Mod_File$datapath)
        pred_LA<-as.data.frame(predict(LA_Model, Val_Data, type = 'prob'))
        LA_ROC<-pROC::roc(label,pred_LA$Non_Responder)
        
        LA_ROCPlot<- ggroc(LA_ROC, legacy.axes = TRUE) +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
          theme_bw()
        
        output$Run_LA_ROCPlot<-renderPlot({return(LA_ROCPlot)})
        
        output$downloadRun_LA_ROCPlot <- downloadHandler(filename = function(){paste("ElasticNet_ROC_Plot",'.png',sep='')},
                                                         content = function(file){
                                                           ggsave(file,plot= LA_ROCPlot, width = 15, height = 10)})
        
        
        # ROC_List<-append(ROC_List, LA_ROC)
        # ROC_Names<-append(ROC_Names, "Elastic Net")
        
        ROC_List[[pos]] <- LA_ROC
        ROC_Names[[pos]] <- "ELR"
        pos = pos + 1
        
        #confusion matrxi
        LA_ConfMat <- as.data.frame(confusionMatrix(LA_Model, "none")$table)
        
        # Final Elastic net model table:
        LA_Table <- data.frame(Model = "Elastic Net Liner Regression Model",
                               AUC = auc(LA_ROC),
                               Accuracy = max(LA_Model$results$Accuracy)*100,
                               Sensitivity = (LA_ConfMat[4, 3]/(LA_ConfMat[4, 3] + LA_ConfMat[3, 3])),
                               Specificity = (LA_ConfMat[1, 3]/(LA_ConfMat[1, 3] + LA_ConfMat[2, 3])),
                               TP_per = (LA_ConfMat[4, 3]/(LA_ConfMat[4, 3] + LA_ConfMat[3, 3]))*100,
                               FP_per = (LA_ConfMat[2, 3]/(LA_ConfMat[1, 3] + LA_ConfMat[2, 3]))*100,
                               TN_per = (LA_ConfMat[1, 3]/(LA_ConfMat[1, 3] + LA_ConfMat[2, 3]))*100,
                               FN_per = (LA_ConfMat[3, 3]/(LA_ConfMat[4, 3] + LA_ConfMat[3, 3]))*100,
                               stringsAsFactors = FALSE)
        
        ROC_Table <- rbind(ROC_Table, LA_Table)
        
        
      }
      
      if (input$BC_Run == TRUE){
        BC_Model <- readRDS(input$BC_Mod_File$datapath)
        pred_BC<-as.data.frame(predict(BC_Model, Val_Data, type = 'prob'))
        BC_ROC<-pROC::roc(label,pred_BC$Non_Responder)
        
        BC_ROCPlot<- ggroc(BC_ROC, legacy.axes = TRUE) +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
          theme_bw()
        
        output$Run_BC_ROCPlot<-renderPlot({return(BC_ROCPlot)})
        
        output$downloadRun_BC_ROCPlot <- downloadHandler(filename = function(){paste("BayesGLM_ROC_Plot",'.png',sep='')},
                                                         content = function(file){
                                                           ggsave(file,plot= BC_ROCPlot, width = 15, height = 10)})
        
        
        # ROC_List<-append(ROC_List, BC_ROC)
        # ROC_Names<-append(ROC_Names, "Bayes")
        # 
        ROC_List[[pos]] <- BC_ROC
        ROC_Names[[pos]] <- "Bayes"
        pos = pos + 1
        
        #confusion matrix
        BC_ConfMat <- as.data.frame(confusionMatrix(BC_Model, "none")$table)
        
        # Final Elastic net model table:
        BC_Table <- data.frame(Model = "Bayes Linear Regression Model",
                               AUC = auc(BC_ROC),
                               Accuracy = max(BC_Model$results$Accuracy)*100,
                               Sensitivity = (BC_ConfMat[4, 3]/(BC_ConfMat[4, 3] + BC_ConfMat[3, 3])),
                               Specificity = (BC_ConfMat[1, 3]/(BC_ConfMat[1, 3] + BC_ConfMat[2, 3])),
                               TP_per = (BC_ConfMat[4, 3]/(BC_ConfMat[4, 3] + BC_ConfMat[3, 3]))*100,
                               FP_per = (BC_ConfMat[2, 3]/(BC_ConfMat[1, 3] + BC_ConfMat[2, 3]))*100,
                               TN_per = (BC_ConfMat[1, 3]/(BC_ConfMat[1, 3] + BC_ConfMat[2, 3]))*100,
                               FN_per = (BC_ConfMat[3, 3]/(BC_ConfMat[4, 3] + BC_ConfMat[3, 3]))*100,
                               stringsAsFactors = FALSE)
        
        ROC_Table <- rbind(ROC_Table, BC_Table)
      }
      
      
      Final_ROCPlot<-ggroc(ROC_List, legacy.axes = TRUE) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
        scale_color_discrete(name = "ML_Models", labels = ROC_Names) +
        theme_bw()
      
      output$ROC_Plot<-renderPlot({return(Final_ROCPlot)})
      
      output$downloadRun_ROCPlot <- downloadHandler(filename = function(){paste("ROC_Plot",'.png',sep='')},
                                                    content = function(file){
                                                      ggsave(file,plot= Final_ROCPlot, width = 15, height = 10)})
      
      
      #round table and remove first column
      ROC_Table<-ROC_Table[-1,]
      ROC_Table$Accuracy <- round(as.numeric(ROC_Table$Accuracy), 0)
      ROC_Table$AUC <- round(as.numeric(ROC_Table$AUC), 2)
      ROC_Table$Sensitivity <- round(ROC_Table$Sensitivity, 2)
      ROC_Table$Specificity <- round(ROC_Table$Specificity, 2)
      ROC_Table$TP_per <- round(as.numeric(ROC_Table$TP_per), 0)
      ROC_Table$FP_per <- round(as.numeric(ROC_Table$FP_per), 0)
      ROC_Table$TN_per <- round(as.numeric(ROC_Table$TN_per), 0)
      ROC_Table$FN_per <- round(as.numeric(ROC_Table$FN_per), 0)
      
      output$ROC_Table<-renderDataTable(ROC_Table, rownames = FALSE)
      
      output$downloadRun_ROCTable <- downloadHandler(filename = function(){"ROC_Accuracy_Table.csv"}, 
                                                     content = function(fname){
                                                       write.csv(ROC_Table, fname)})
      output$error_ML_Run <- renderUI(NULL)
    },
    error = function(e) {
      output$error_ML_Run <- renderUI({
        div(
          class = "alert alert-danger", role = "alert",
          strong("Error: "), "An unexpected error occurred,
          Please check the format of your file and try again.
          For more information, please look at the example file available through the 
          download example file button.", e$message
        )
      })
    }
    )
    
    
  })
  
  reactive({
    inFile <- input$SPM_val
    
    print(inFile)
    
    if(is.null(inFile))
      
      return(NULL)
    
    # Open the file:
    
    df_prime <- read.table(inFile$datapath, 
                           header=input$header,
                           sep=input$sep)
    
    
    required_columns <- c("adj.P.Val", "vip.scores")
    
    # Iterate over the required columns
    for (col in required_columns) {
      # Check if the column is missing
      if (!(col %in% colnames(df_prime))) {
        # Add the missing column with NA values
        df_prime[[col]] <- NA
      }
    }
    
    # Assign the modified dataframe to df
    df <- df_prime
    
    df$adj.P.Val[is.na(df$adj.P.Val)] <- 1
    df$vip.scores[is.na(df$vip.scores)] <- 0
    
    # Using dictionary to add subscripts to SPM names in the file uploaded by the user
    df <- df %>% mutate(name = recode(name, !!!SPM_sub))
    
    
  })
  
  # Button to begin network generation process ----
  observeEvent(input$Go, {
    
    inFile <- input$SPM_val
    
    print(inFile)
    
    if(is.null(inFile))
      
      return(NULL)
    
    df_prime <- read.csv(inFile$datapath,
                         header = input$header,
                         sep = input$sep)
    
    required_columns <- c("adj.P.Val", "vip.scores")
    
    # Iterate over the required columns
    for (col in required_columns) {
      # Check if the column is missing
      if (!(col %in% colnames(df_prime))) {
        # Add the missing column with NA values
        df_prime[[col]] <- NA
      }
    }
    
    # Assign the modified dataframe to df
    df <- df_prime
    
    df$adj.P.Val[is.na(df$adj.P.Val)] <- 1
    df$vip.scores[is.na(df$vip.scores)] <- 0
    
    
    
    
    # Building the data for the Node color changes in Human model ----
    EPA_nodes <- reactiveVal(EPA_data)
    AA_nodes <- reactiveVal(AA_data)
    DHA_nodes <- reactiveVal(DHA_data)
    DPA_nodes <- reactiveVal(DPA_data)
    
    # Building the data for the Node color changes in Mouse model ----
    EPA_data_mouse <- EPA_data
    EPA_nodes_mouse <- reactiveVal(EPA_data)
    AA_nodes_mouse <- reactiveVal(AA_data)
    DHA_nodes_mouse <- reactiveVal(DHA_data)
    DPA_nodes_mouse <- reactiveVal(DPA_data)
    
    # Reactive value to store clicked node ID
    clicked_node_id <- reactiveVal(NULL)
    
    observeEvent(input$clickNode, {
      clicked_node_id(input$clickNode)
    })
    
    observeEvent(input$update_node_color, {
      if (!is.null(clicked_node_id())) {
        selected_node_id <- clicked_node_id()
        selected_color <- input$col1
        
        updated_EPA_nodes <- EPA_nodes()
        updated_EPA_nodes[updated_EPA_nodes$id == selected_node_id, "color.background"] <- selected_color
        
        updated_AA_nodes <- AA_nodes()
        updated_AA_nodes[updated_AA_nodes$id == selected_node_id, "color.background"] <- selected_color
        
        updated_DHA_nodes <- DHA_nodes()
        updated_DHA_nodes[updated_DHA_nodes$id == selected_node_id, "color.background"] <- selected_color
        
        updated_DPA_nodes <- DPA_nodes()
        updated_DPA_nodes[updated_DPA_nodes$id == selected_node_id, "color.background"] <- selected_color
        
        EPA_nodes(updated_EPA_nodes)
        AA_nodes(updated_AA_nodes)
        DHA_nodes(updated_DHA_nodes)
        DPA_nodes(updated_DPA_nodes)
        
        updated_EPA_nodes_mouse <- EPA_nodes_mouse()
        updated_EPA_nodes_mouse[updated_EPA_nodes_mouse$id == selected_node_id, "color.background"] <- selected_color
        
        updated_AA_nodes_mouse <- AA_nodes_mouse()
        updated_AA_nodes_mouse[updated_AA_nodes_mouse$id == selected_node_id, "color.background"] <- selected_color
        
        updated_DHA_nodes_mouse <- DHA_nodes_mouse()
        updated_DHA_nodes_mouse[updated_DHA_nodes_mouse$id == selected_node_id, "color.background"] <- selected_color
        
        updated_DPA_nodes_mouse <- DPA_nodes_mouse()
        updated_DPA_nodes_mouse[updated_DPA_nodes_mouse$id == selected_node_id, "color.background"] <- selected_color
        
        EPA_nodes_mouse(updated_EPA_nodes_mouse)
        AA_nodes_mouse(updated_AA_nodes_mouse)
        DHA_nodes_mouse(updated_DHA_nodes_mouse)
        DPA_nodes_mouse(updated_DPA_nodes_mouse)
      }
    })
    
    # Building the data for edge colour changes in Human model: ----
    EPA_edges <- reactiveVal(EPA_direction)
    AA_edges <- reactiveVal(AA_direction)
    DHA_edges <- reactiveVal(DHA_direction)
    DPA_edges <- reactiveVal(DPA_direction)
    
    # Building the data for edge colour changes in Mouse model: ----
    EPA_edges <- reactiveVal(EPA_direction)
    AA_edges <- reactiveVal(AA_direction)
    DHA_edges <- reactiveVal(DHA_direction)
    DPA_edges <- reactiveVal(DPA_direction)
    
    # Reactive value to store clicked node ID
    clicked_edge_id <- reactiveVal(NULL)
    
    observeEvent(input$clickEdge, {
      clicked_edge_id(input$clickEdge)
    })
    
    observeEvent(input$update_edge_color, {
      if (!is.null(clicked_edge_id())) {
        selected_edge_id <- clicked_edge_id()
        selected_edge_color <- input$col2
        
        updated_EPA_edges <- EPA_edges()
        updated_EPA_edges[updated_EPA_edges$id == selected_edge_id, "color"] <- selected_edge_color
        
        updated_AA_edges <- AA_edges()
        updated_AA_edges[updated_AA_edges$id == selected_edge_id, "color"] <- selected_edge_color
        
        updated_DHA_edges <- DHA_edges()
        updated_DHA_edges[updated_DHA_edges$id == selected_edge_id, "color"] <- selected_edge_color
        
        updated_DPA_edges <- DPA_edges()
        updated_DPA_edges[updated_DPA_edges$id == selected_edge_id, "color"] <- selected_edge_color
        
        EPA_edges(updated_EPA_edges)
        AA_edges(updated_AA_edges)
        DHA_edges(updated_DHA_edges)
        DPA_edges(updated_DPA_edges)
      }
    })
    
    # Reactive expression to calculate unmatched_rows
    unmatched_rows <- reactive({
      
      
      EPA_dict <- c("RvE1", "RvE2", "RvE3", "RvE4")
      AA_dict <- c("TXB2", "PGD2", "PGE2", "PGF2a", 
                   "LXA4", "15-oxo-LXA4", "13,14-dehydro-15-oxo LXA4", "15-epi-LXA4", 
                   "15-epi-LXB4", "LXB4", "5S,15S-diHETE", "6-trans-LTB4", 
                   "6-trans-12-epi-LTB4", "LTB4", "20-OH-LTB4", "20-COOH-LTB4", 
                   "5S,12-diHETE", "LTC4", "LTD4", "LTE4")
      DHA_dict <- c("PCTR1", "PCTR2", "PCTR3", "PD1", "17R-PD1", "PDX", "22-OH-PD1", "RvD1", "RvD2", 
                    "RvD3", "RvD4", "RvD5", "RvD6", "17R- RvD3", "17R- RvD1", "4S,14S-diHDHA", 
                    "7S,14S-diHDHA", "MaR2", "MaR1", "14-oxo-MaR1", "22-OH-MaR1", "MCTR1", 
                    "MCTR2", "MCTR3")
      DPA_dict <- c("PD1 n-3 DPA", "PD2 n-3 DPA", "10S,17S-diHPA", "RvT1", "RvT2", "RvT3", "RvT4", 
                    "RvD5 n-3 DPA", "RvD2 n-3 DPA", "RvD1 n-3 DPA", "MaR1 n-3 DPA", 
                    "MaR2 n-3 DPA", "7S,14S-diHDPA")
      
      
      SPM_dict <- c(EPA_dict, AA_dict, DHA_dict, DPA_dict)
      
      # Ensure that input_names is reactive or obtained from user input
      input_names <- df[[1]]  # Replace this with your actual input handling if needed

      unused_SPMS <- SPM_dict[!SPM_dict %in% input_names]
      unmatched_SPMs <- input_names[!input_names %in% SPM_dict]
      unmatched_SPMs <- data.frame(unmatched_SPMs)
      
      # Identifying EPA and DHA SPMs
      Identified_EPA_SPMs <- EPA_direction %>%
        filter(to %in% EPA_dict) %>%
        select("Unmatched_SPMs" = to, "Derived From" = from) %>%
        mutate("Fatty acid precursor" = "EPA")

      Identified_DHA_SPMs <- DHA_direction %>%
        filter(to %in% DHA_dict) %>%
        select("Unmatched_SPMs" = to, "Derived From" = from) %>%
        mutate("Fatty acid precursor" = "DHA")
      
      Identified_DPA_SPMs <- DPA_direction %>%
        filter(to %in% DPA_dict) %>%
        select("Unmatched_SPMs" = to, "Derived From" = from) %>%
        mutate("Fatty acid precursor" = "DPA")
      
      Identified_AA_SPMs <- AA_direction %>%
        filter(to %in% AA_dict) %>%
        select("Unmatched_SPMs" = to, "Derived From" = from) %>%
        mutate("Fatty acid precursor" = "AA")
      
      # Combine EPA and DHA SPMs into a single data frame
      combined_df <- rbind(Identified_EPA_SPMs, Identified_DHA_SPMs, Identified_AA_SPMs, Identified_DPA_SPMs) %>%
        select("Fatty acid precursor", everything())
      
      # Filter combined_df to keep only unmatched rows
      unmatched <- combined_df %>%
        filter(Unmatched_SPMs %in% unused_SPMS) 
      
      
      return(unmatched)
      
    })
    
    unmatched_SPMs <- reactive({
      
      
      EPA_dict <- c("RvE1", "RvE2", "RvE3", "RvE4")
      AA_dict <- c("TXB2", "PGD2", "PGE2", "PGF2a", 
                   "LXA4", "15-oxo-LXA4", "13,14-dehydro-15-oxo LXA4", "15-epi-LXA4", 
                   "15-epi-LXB4", "LXB4", "5S,15S-diHETE", "6-trans-LTB4", 
                   "6-trans-12-epi-LTB4", "LTB4", "20-OH-LTB4", "20-COOH-LTB4", 
                   "5S,12S-diHETE", "LTC4", "LTD4", "LTE4")
      DHA_dict <- c("PCTR1", "PCTR2", "PCTR3", "PD1", "17R-PD1", "PDX", "22-OH-PD1", "RvD1", "RvD2", 
                    "RvD3", "RvD4", "RvD5", "RvD6", "17R-RvD3", "17R-RvD1", "4S,14S-diHDHA", 
                    "7S,14S-diHDHA", "MaR2", "MaR1", "14-oxo-MaR1", "22-OH-MaR1", "MCTR1", 
                    "MCTR2", "MCTR3")
      DPA_dict <- c("PD1n-3DPA", "PD2n-3DPA", "10S,17S-diHDPA", "RvT1", "RvT2", "RvT3", "RvT4", 
                    "RvD5n-3DPA", "RvD2n-3DPA", "RvD1n-3DPA", "MaR1n-3DPA", 
                    "MaR2n-3DPA", "7S,14S-diHDPA")
      
      
      SPM_dict <- c(EPA_dict, AA_dict, DHA_dict, DPA_dict)
      
      input_names <- df[[1]]  
      
      unmatched_SPMs <- input_names[!input_names %in% SPM_dict]
      
      unmatched_df <- df %>%
        filter(input_names %in% unmatched_SPMs)
      
      return(unmatched_df)
      
    })
    
    # Render the unmatched rows in the DataTable output
    output$unused_database <- renderDataTable({unmatched_rows()})
    output$unused_file <- renderDataTable({unmatched_SPMs()})
    output$unused_database2 <- renderDataTable({unmatched_rows()})
    output$unused_file2 <- renderDataTable({unmatched_SPMs()})
    
    
    
    # Building the for loops: ----
    observe({
      

      # Output the results
     
      # EPA network for loops for automatic changes: ----
      
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  EPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  EPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  EPA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                EPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                EPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                EPA_data$shape[row_index] <- "triangleDown"
              } else {
                EPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(EPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  EPA_data$size[row_index] <- y
                } else {
                  EPA_data$size[row_index] <- 10
                }
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                EPA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05){ 
                EPA_direction$color[row_index] <- "grey"
              } else {
                EPA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05 && log2FC != 0) {
                EPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  EPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  EPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  EPA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                EPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                EPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                EPA_data$shape[row_index] <- "triangleDown"
              } else {
                EPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(EPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  EPA_data$size[row_index] <- y
                } else {
                  EPA_data$size[row_index] <- 10
                } 
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC != 0) {
                  EPA_direction$color[row_index] <- "black"
                } else if (log2FC == 0){
                  EPA_direction$color[row_index] <- "grey"
                } else {
                  EPA_direction$color[row_index] <- "lightgrey"
                }
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                EPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      EPA_nodes(EPA_data)
      EPA_edges(EPA_direction)
      
      # AA network for loops for automatic changes: ----
      
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  AA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  AA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  AA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                AA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                AA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                AA_data$shape[row_index] <- "triangleDown"
              } else {
                AA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(AA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC != 0) {
                x <- df$logFC[i]
                y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                if (y > 40) {
                  y = 40
                }
                AA_data$size[row_index] <- y
              } else {
                AA_data$size[row_index] <- 10
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                AA_direction$color[row_index] <- "black"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05 && log2FC != 0) {
                AA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  AA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  AA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  AA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                AA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                AA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                AA_data$shape[row_index] <- "triangleDown"
              } else {
                AA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(AA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC != 0) {
                x <- df$logFC[i]
                y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                if (y > 40) {
                  y = 40
                }
                AA_data$size[row_index] <- y
              } else {
                AA_data$size[row_index] <- 10
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                AA_direction$color[row_index] <- "black"
              } else if (vip_scores == 0) {
                AA_direction$color[row_index] <- "grey"
              } else {
                AA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                AA_direction$width[row_index] <- 2.5
              } else {
                
              }
            }
          }
        }
      }
      AA_nodes(AA_data)
      AA_edges(AA_direction)
      
      # DHA network for loops for automatic changes: ----
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  DHA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DHA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DHA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                DHA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DHA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DHA_data$shape[row_index] <- "triangleDown"
              } else {
                DHA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DHA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DHA_data$size[row_index] <- y
                } else {
                  DHA_data$size[row_index] <- 10
                }
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                DHA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05){
                DHA_direction$color[row_index] <- "grey"
              } else {
                DHA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                DHA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  DHA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DHA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DHA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                DHA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DHA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DHA_data$shape[row_index] <- "triangleDown"
              } else {
                DHA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DHA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DHA_data$size[row_index] <- y
                } else {
                  DHA_data$size[row_index] <- 10
                } 
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                DHA_direction$color[row_index] <- "black"
              } else if (vip_scores < 1) {
                DHA_direction$color[row_index] <- "grey"
              } else {
                DHA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                DHA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      DHA_nodes(DHA_data)
      DHA_edges(DHA_direction)
      
      # n-3 DPA network for loops for automatic changes: ----
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  DPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DPA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                DPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DPA_data$shape[row_index] <- "triangleDown"
              } else {
                DPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DPA_data$size[row_index] <- y
                } else {
                  DPA_data$size[row_index] <- 10
                }
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                DPA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05) {
                DPA_direction$color[row_index] <- "grey"
              } else {
                DPA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05 && log2FC != 0) {
                DPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  DPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DPA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                DPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DPA_data$shape[row_index] <- "triangleDown"
              } else {
                DPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC > 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DPA_data$size[row_index] <- y
                } else {
                  DPA_data$size[row_index] <- 10
                } 
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                DPA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05) {
                DHA_direction$color[row_index] <- "grey"
              } else {
                DHA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                DPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      DPA_nodes(DPA_data)
      DPA_edges(DPA_direction)
    })
    
    # Output for networks: ---- 
    output$network <- renderVisNetwork({
      
      # Creating the EPA network template ----
      
      # Building coordinate matrix for the node positions of the Network:
      EPA_coordinates = EPA_data[, c("x", "y")]
      EPA_coord_matrix = data.matrix(EPA_coordinates)
      
      EPA_visnet <- visNetwork(nodes = EPA_nodes(), EPA_edges(), height = "1000px", 
                               width = "100%", 
                               style = "font-family:Calibri body") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = EPA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20)) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "EPA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "")%>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      # Creating the AA network template: ----
      
      # Building coordinate matrix for the node positions of the Network:
      AA_coordinates = AA_data[, c("x", "y")]
      AA_coord_matrix = data.matrix(AA_coordinates)
      
      AA_visnet <- visNetwork(nodes = AA_nodes(), edges = AA_edges(), height = "1000px", 
                              width = "100%") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = AA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20)) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "AA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "") %>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      
      # Creating the DHA network template: ----
      
      # Building coordinate matrix for the node positions of the Network:
      DHA_coordinates = DHA_data[, c("x", "y")]
      DHA_coord_matrix = data.matrix(DHA_coordinates)
      
      DHA_visnet <- visNetwork(nodes = DHA_nodes(), edges = DHA_edges(), height = "1000px", 
                               width = "100%", 
                               style = "font-family:Calibri body") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = DHA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20), physics = FALSE) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "DHA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "")%>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      # Creating the n-3 DPA network template: ----
      
      # Building coordinate matrix for the node positions of the Network:
      DPA_coordinates = DPA_data[, c("x", "y")]
      DPA_coord_matrix = data.matrix(DPA_coordinates)
      
      DPA_visnet <- visNetwork(nodes = DPA_nodes(), edges = DPA_edges(), height = "1000px", 
                               width = "100%", 
                               style = "font-family:Calibri body") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = DPA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20)) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "DPA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "")%>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      
      
      
      
      
      # Network selection connection with button in front end: ----
      
      networks <- switch(input$selection,
                         "AA Network" = AA_visnet,
                         "EPA Network" = EPA_visnet,
                         "DHA Network" = DHA_visnet,
                         "n-3 DPA Network" = DPA_visnet)
      
      
    })
    
    observe({
      
      # Changing the pathway for the mouse model from ALOX12 to ALOX12/15:
      AA_data$label[AA_data$label == "ALOX12"] <- "ALOX12/15"
      DHA_data$label[DHA_data$label == "ALOX12"] <- "ALOX12/15"
      DPA_data$label[DPA_data$label == "ALOX12"] <- "ALOX12/15"
      
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  EPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  EPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  EPA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                EPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                EPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                EPA_data$shape[row_index] <- "triangleDown"
              } else {
                EPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(EPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  EPA_data$size[row_index] <- y
                } else {
                  EPA_data$size[row_index] <- 10
                }
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                EPA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05){ 
                EPA_direction$color[row_index] <- "grey"
              } else {
                EPA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05 && log2FC != 0) {
                EPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  EPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  EPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  EPA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                EPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(EPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                EPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                EPA_data$shape[row_index] <- "triangleDown"
              } else {
                EPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(EPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  EPA_data$size[row_index] <- y
                } else {
                  EPA_data$size[row_index] <- 10
                } 
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC != 0) {
                  EPA_direction$color[row_index] <- "black"
                } else if (log2FC == 0){
                  EPA_direction$color[row_index] <- "grey"
                } else {
                  EPA_direction$color[row_index] <- "lightgrey"
                }
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(EPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                EPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      EPA_nodes_mouse(EPA_data)
      EPA_edges(EPA_direction)
      
      # AA network for loops for automatic changes: ----
      
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  AA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  AA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  AA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                AA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                AA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                AA_data$shape[row_index] <- "triangleDown"
              } else {
                AA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(AA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC != 0) {
                x <- df$logFC[i]
                y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                if (y > 40) {
                  y = 40
                }
                AA_data$size[row_index] <- y
              } else {
                AA_data$size[row_index] <- 10
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                AA_direction$color[row_index] <- "black"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05 && log2FC != 0) {
                AA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  AA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  AA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  AA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                AA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(AA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                AA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                AA_data$shape[row_index] <- "triangleDown"
              } else {
                AA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(AA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC != 0) {
                x <- df$logFC[i]
                y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                if (y > 40) {
                  y = 40
                }
                AA_data$size[row_index] <- y
              } else {
                AA_data$size[row_index] <- 10
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                AA_direction$color[row_index] <- "black"
              } else if (vip_scores == 0) {
                AA_direction$color[row_index] <- "grey"
              } else {
                AA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(AA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                AA_direction$width[row_index] <- 2.5
              } else {
                
              }
            }
          }
        }
      }
      AA_nodes_mouse(AA_data)
      AA_edges(AA_direction)
      
      # DHA network for loops for automatic changes: ----
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  DHA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DHA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DHA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                DHA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DHA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DHA_data$shape[row_index] <- "triangleDown"
              } else {
                DHA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DHA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DHA_data$size[row_index] <- y
                } else {
                  DHA_data$size[row_index] <- 10
                }
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                DHA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05){
                DHA_direction$color[row_index] <- "grey"
              } else {
                DHA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                DHA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  DHA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DHA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DHA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                DHA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DHA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DHA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DHA_data$shape[row_index] <- "triangleDown"
              } else {
                DHA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DHA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DHA_data$size[row_index] <- y
                } else {
                  DHA_data$size[row_index] <- 10
                } 
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                DHA_direction$color[row_index] <- "black"
              } else if (vip_scores < 1) {
                DHA_direction$color[row_index] <- "grey"
              } else {
                DHA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DHA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                DHA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      DHA_nodes(DHA_data)
      DHA_edges(DHA_direction)
      
      # n-3 DPA network for loops for automatic changes: ----
      if (input$disp=="pval") {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (adj_p_val < 0.05) {
                if (log2FC > 0 ) {
                  DPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DPA_data$color.background[row_index] <- "grey"
                }
              } else if (adj_p_val > 0.05) {
                DPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DPA_data$shape[row_index] <- "triangleDown"
              } else {
                DPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                if (log2FC != 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DPA_data$size[row_index] <- y
                } else {
                  DPA_data$size[row_index] <- 10
                }
              } 
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05) {
                DPA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05) {
                DPA_direction$color[row_index] <- "grey"
              } else {
                DPA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (adj_p_val < 0.05 && log2FC != 0) {
                DPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      else {
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            # If matching rows are found, update colors for all matching rows
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              # Updating color based on adjusted p value and log2FC
              if (vip_scores > 1) {
                if (log2FC > 0 ) {
                  DPA_data$color.background[row_index] <- "red"
                } else if (log2FC < 0) {
                  DPA_data$color.background[row_index] <- "blue"
                } else if (log2FC == 0) {
                  DPA_data$color.background[row_index] <- "grey"
                }
              } else if (vip_scores < 1) {
                DPA_data$color.background[row_index] <- "grey"
              } else {
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          
          matching_rows <- which(DPA_data$id == id_value)
          
          if (length(matching_rows) > 0) {
            for (row_index in matching_rows) {
              
              vip_scores <- df$vip.scores[i]
              log2FC <- df$logFC[i]
              
              if (log2FC > 0) {
                DPA_data$shape[row_index] <- "triangle"
              } else if (log2FC < 0) {
                DPA_data$shape[row_index] <- "triangleDown"
              } else {
                DPA_data$shape[row_index] <- "dot"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          
          id_value <- df$name[i]
          matching_row <- which(DPA_data$id == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                if (log2FC > 0) {
                  x <- df$logFC[i]
                  y <- 2.5573 * x^2 - 0.1017 * x + 8.8962
                  if (y > 40) {
                    y = 40
                  }
                  DPA_data$size[row_index] <- y
                } else {
                  DPA_data$size[row_index] <- 10
                } 
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1) {
                DPA_direction$color[row_index] <- "black"
              } else if (adj_p_val > 0.05) {
                DHA_direction$color[row_index] <- "grey"
              } else {
                DHA_direction$color[row_index] <- "lightgrey"
              }
            }
          }
        }
        for (i in 1:nrow(df)) {
          id_value <- df$name[i]
          
          # Find the corresponding row in log2FC_data
          matching_row <- which(DPA_direction$to == id_value)
          
          # Check if there are any matching rows
          if (length(matching_row) > 0) {
            for (row_index in matching_row) {
              
              vip_scores <- df$vip.scores[i]
              adj_p_val <- df$adj.P.Val[i]
              log2FC <- df$logFC[i]
              
              if (vip_scores > 1 && log2FC != 0) {
                DPA_direction$width[row_index] <- 2.5
              } else {
              }
            }
          }
        }
      }
      DPA_nodes(DPA_data)
      DPA_edges(DPA_direction)
      
    })
    
    output$network_mouse <- renderVisNetwork({
      
      
      # Creating the EPA network template ----
      
      # Building coordinate matrix for the node positions of the Network:
      EPA_coordinates = EPA_data[, c("x", "y")]
      EPA_coord_matrix = data.matrix(EPA_coordinates)
      
      EPA_visnet <- visNetwork(nodes = EPA_nodes(), EPA_edges(), height = "1000px", 
                               width = "100%", 
                               style = "font-family:Calibri body") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = EPA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20)) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "EPA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "")%>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      # Creating the AA network template: ----
      
      # Building coordinate matrix for the node positions of the Network:
      AA_coordinates = AA_data[, c("x", "y")]
      AA_coord_matrix = data.matrix(AA_coordinates)
      
      AA_visnet <- visNetwork(nodes = AA_nodes(), edges = AA_edges(), height = "1000px", 
                              width = "100%") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = AA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20)) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "AA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "") %>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      
      # Creating the DHA network template: ----
      
      # Building coordinate matrix for the node positions of the Network:
      DHA_coordinates = DHA_data[, c("x", "y")]
      DHA_coord_matrix = data.matrix(DHA_coordinates)
      
      DHA_visnet <- visNetwork(nodes = DHA_nodes(), edges = DHA_edges(), height = "1000px", 
                               width = "100%", 
                               style = "font-family:Calibri body") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = DHA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20), physics = FALSE) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "DHA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "")%>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      # Creating the n-3 DPA network template: ----
      
      # Building coordinate matrix for the node positions of the Network:
      DPA_coordinates = DPA_data[, c("x", "y")]
      DPA_coord_matrix = data.matrix(DPA_coordinates)
      
      DPA_visnet <- visNetwork(nodes = DPA_nodes(), edges = DPA_edges(), height = "1000px", 
                               width = "100%", 
                               style = "font-family:Calibri body") %>% 
        visIgraphLayout(type = "full", layout = "layout.norm", 
                        layoutMatrix = DPA_coord_matrix, physics = FALSE) %>%  
        visNodes(font = list(size = 20)) %>%
        visOptions(manipulation = TRUE) %>%
        visPhysics(enabled = FALSE, repulsion = FALSE, 
                   stabilization = TRUE) %>%
        visExport(type = "png", name = "DPA-network", 
                  float = "left", label = "Export Network", 
                  background = "white", style = "")%>%
        visInteraction(multiselect = TRUE) %>%
        visEvents(click = "function(properties) {
                  Shiny.onInputChange('clickNode', properties.nodes[0]);
                  Shiny.onInputChange('clickEdge', properties.edges[0]);
                }")
      
      
      
      
      
      
      # Network selection connection with button in front end: ----
      
      networks <- switch(input$selection,
                         "AA Network" = AA_visnet,
                         "EPA Network" = EPA_visnet,
                         "DHA Network" = DHA_visnet,
                         "n-3 DPA Network" = DPA_visnet)
      
      
    })
    
    
  })
  
  
})

shinyApp(ui = ui, server = server)