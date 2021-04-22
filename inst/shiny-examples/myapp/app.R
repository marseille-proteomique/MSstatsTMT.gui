library(shiny)
library(shinydashboard)
library(tidyverse)
library(MSstatsTMT)
library(MSstats)
library(stringr)
library(matlab)
library(rio)
library(DT)
library(shinyMatrix)
library(shinycssloaders)
library(shinybusy)
library(gghighlight)
library(ggpubr)

#increase the max request size for uploading files
options(shiny.maxRequestSize = 1000*1024^2)
#set options for the spinner when things are loading
options(spinner.color = "#518CE2", spinner.color.background = "000000", spinner.size = 2)

ui <- dashboardPage(skin = "blue",

                    dashboardHeader(title = "MSstat TMT analysis", titleWidth = 300),

                    dashboardSidebar(width = 300,
                      sidebarMenu(
                        menuItem("Readme", tabName = "read", icon = icon("readme")),
                        menuItem("Your annotation file", tabName = "annot", icon = icon("edit")),
                        menuItem("MSstats format", tabName = "MS_format", icon = icon("file-upload")),
                        menuItem("Perform quantification", tabName = "MS_quant", icon = icon("calculator")),
                        menuItem("Group comparison", tabName = "Group_comp", icon = icon("not-equal")),
                        menuItem("Profile and QC plot", tabName = "Plot", icon = icon("chart-line"))
                      )
                    ),

                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "read",
                                h1("Readme"),
                                tags$br(),

                                HTML("<h4>This app allow any users to use the MSstatsTMT package in a 'user-friendly' way. <br>
                                      MSstatsTMT is a R package which provides statistical tools for detecting differentially
                                      abundant proteins in shotgun mass spectrometry-based proteomic experiments
                                      with tandem mass tag (TMT) labeling. The quantification is based on the peptide-level. <br>
                                      For more information, please check <a href='http://msstats.org/msstatstmt/'>MSstatsTMT vignette</a></h4>"),

                                tags$hr(),
                                h2("MSstats format"),
                                tags$br(),

                                HTML("<h5>The MSstatsTMT package needs file input; it can accept output from  four different software : <br>
                                     - Proteome Discoverer <br>
                                     - MaxQuant <br>
                                     - SpectroMine <br>
                                     - OpenMS <br>
                                     On the tab 'MSstats format', you can upload your data from this four software. You will see that for
                                     all the software, except OpenMS, you need to add an annotation file. This file had to be done by the user
                                     before the analysis. It has to contains seven columns named Run, Fraction, TechRepMixture, Channel,
                                     Condition, BioReplicate and Mixture (see the example below or evntually check
                                     <a href='http://www.bioconductor.org/packages/release/bioc/vignettes/MSstatsTMT/inst/doc/MSstatsTMT.html'>MSstatsTMT detailed vignette</a>.
                                     <br>
                                     You can create this file in the tab 'Your annotation file'. Note that if you use OpenMS data, you don't need an annotation file.
                                     In this tab you have to import your output file from the software you use, put the number of the channels, batch, fraction and mixture you have.
                                     You can tell the condition of each channel with the input matrix. 'Norm' is your control, it has to be name exactly like 'Norm'
                                     if you want to use it as a control. You can choose the name of the channel, but it has be coherant with your output file !
                                     See the examples donx below for more precision. Once you created it, you can download it in xlsx format and use it in the app.
                                     <br>
                                     Once your uploaded your files, you can select you can select your filtering criterias.
                                     Finally, you can click on the submit button to start the calculation; message should appear, informing you of the progress.<br>
                                     When the conversion is done, the table output will appear. You can continue your analysis or download this tab and use it in the app
                                     to continue your analysis later.</h5>"),
                                tags$br(),

                                fluidRow(column(3, checkboxInput("annot_exmq", "See a basic annotation file example for MaxQuant")),
                                         column(3, checkboxInput("annot_exmine", "See a basic annotation file example for Spectromine")),
                                         column(3, checkboxInput("annot_expd", "See a basic annotation file example for Proteome Discoverer")),
                                         column(3, checkboxInput("annot_ex2", "See a more complex annotation file example for MaxQuant"))
                                         ),

                                conditionalPanel(condition = "input.annot_exmq", DT::dataTableOutput("annot_ex_mq")),
                                conditionalPanel(condition = "input.annot_exmine", DT::dataTableOutput("annot_ex_mine")),
                                conditionalPanel(condition = "input.annot_expd", DT::dataTableOutput("annot_ex_pd")),
                                conditionalPanel(condition = "input.annot_ex2", DT::dataTableOutput("annot_ex_comp")),
                                tags$hr(),

                                h2("Quantification"),
                                tags$br(),

                                HTML("<h5>Once you converted your file in MSstats format, you can start the quantification.
                                      You can continue the analysis with the output from the app or import your own file; if so it needs to be in MSstats format. <br>
                                      Then you can choose your quantification criterias. You can choose to perform a global median
                                      normalization on peptide level data; perform a reference channel based normalization
                                      between MS runs, the reference channel need to be annotated by 'Norm'. It will be performed after
                                      protein-level summarization.  If data only has one run, then the function will not perform this step.
                                      You can also choose to remove the 'Norm channel' and the empty channels. <br>
                                      <br>
                                      Step in more details : <br>
                                      First, take log2 intensity and replace value less than 1 with 0 (will replace the -Inf).
                                      Then do the peptide-level global median normalization, calculate the median of each channel and
                                      run. Then calculate the difference between the channel medians and their median, finally
                                      add this difference to all intensities (still the log2 intensity).
                                      Then negative log2 intensity are replaced with NA.
                                      You have to choose the summarization method for protein-level :<br>
                                      - msstats : for each run, use msstats dataProcess
                                      (function from MSstats, see <a href='http://msstats.org'>MSstats vignette</a>).
                                      With this method you can impute missing values by Accelarated failure model; if not
                                      it uses minimum value to impute the missing value for each peptide precursor ion. <br>

                                      - MedianPolish : add NAs to make every protein appear in all channels,
                                      create annotation to make sure there are no missing channels for each protein in data,
                                      same thing for missing PSMs for each run and protein ; then perform medianPolish function (on log2 intensity) <br>

                                      - LogSum : remove all NAs ; calculate the logsum for each protein and channel <br>

                                      - Median : same thing as described above but with median calculation <br>

                                      Next, do the protein normalization : calculate the mean over multiple normalization channels,
                                      calculate difference between median for each abundance mean. Then add this difference to abundance.
                                      Finally get the results, and make sure they are in good format. <br>
                                      <br>
                                      When you choosed your filtering criterias, you can start the quantification with the button.
                                      Message should appear, informing you of the progress. <br>
                                      At last, a table with the results will appear. As for MSstats format, you can continue your analysis or download the results.
                                     </h5>"),

                                tags$hr(),
                                h2("Group comparison"),
                                tags$br(),

                                HTML("<h5>Once the quantification is done, you can go to the tab 'Group comparison'. You can also import your own data. <br>
                                      In this tab you can compare different condition, like 'sick' anc 'control', and obtain a table which contains for each
                                      protein a p-value and a log2 FC. Before proceeding, you can choose your criterias; remove or not the 'Norm' and the
                                      empty channels, choose to do a moderate t-statistic or not and choose the method for adjusting the p-values. Default
                                      is 'BH' for Benjamini-Hochberg. Moreover, you can perform your own comparison with a comparison matrix that you can
                                      make in an interactive way. 0 and the condition will not be taken, 1 and -1 is for comparing. <br>
                                      Once your data are uploaded and your criterias filled, you can start the calculation with the button. Finally, you can download
                                      the results. Note that you can order the proteins according to their p-value or log2 FC by clicking on the arrow next to the
                                      column names of the data table.
                                      <br>
                                      With this output data, you can plot a volcano according to the conditions you want and download it.
                                     </h5>"),

                                tags$hr(),
                                h2("Profile and QC plot"),
                                tags$br(),

                                HTML("<h5>For this tab, you need your quantification data and the MSstatsTMT data format from the first tab.
                                      As for the other tab, you can import your own files. <br>
                                      You can choose to see the profile plot and/or the QC plot of the protein you selected. For the profile plot,
                                      you will see the profile of each peptide from the protein you selected. You can choose to visualize the
                                      summary plot in order to see the protein profile (you can visualize both or one separately).
                                      Once your data are uploaded, you can click on the buttons to see the plots. You can also save the plots with
                                      the download button at the bottom of the tab.
                                     </h5>"),

                                tags$hr(),

                                HTML("<h4>If you have any questions or dealing with some bugs, feel free to contact me at <u>marco.gerault@gmail.com</u>
                                     </h4>")


                                ),

                        tabItem(tabName = "annot",

                                radioButtons("type_outA", "File output from : ",
                                             choices =  c("Proteome Discoverer" = "PD",
                                                          "MaxQuant" = "MQ",
                                                          "SpectroMine" = "SpM"), inline = TRUE, selected = "MQ"),

                                tags$hr(),
                                h2("Select your files in order to create your annotation file"),
                                tags$hr(),

                                conditionalPanel(condition = "input.type_outA != 'MQ'",
                                                          fluidRow(column(5, uiOutput("file1A"))
                                                                   )
                                                          ),
                                conditionalPanel(condition = "input.type_outA == 'MQ'",
                                                 fluidRow(column(5, fileInput("eviA",
                                                                              label= h3("Select the evidence file from MaxQuant output"),
                                                                              accept = c(".xlsx", ".txt", ".csv"))
                                                                 )
                                                          )
                                                 ),

                                fluidRow(column(3, numericInput("nchan", "How many channel do you have ?", value = 10, min = 1, step = 1)),
                                         column(3, numericInput("nfrac", "How many fraction do you have ?", value = 14, min = 1, step = 1)),
                                         column(3, numericInput("nbatch", "How many batch do you have ?", value = 2, min = 1, step = 1)),
                                         column(3, numericInput("ntech", "How many technical replicate mixture
                                                                do you have ?", value = 1, min = 1, step = 1))
                                         ),
                                fluidRow(column(5, uiOutput("condi_matui")),
                                         column(5, uiOutput("batch_matui"))
                                         ),

                                #check if files are well uploaded
                                conditionalPanel(condition = "output.eviA_fileup",
                                                 actionButton("but_anno", "Create your annotation file")),

                                tags$hr(),

                                DT::dataTableOutput("ANNOTATION_out"),

                                conditionalPanel(condition = "input.but_anno",
                                                 downloadButton("save_ANNOTATION", "Download your annotation file"))
                                ),

                        tabItem(tabName = "MS_format",
                                #shinyjs allows to display message from function on the app for better follow-up of calculation
                                shinyjs::useShinyjs(),
                                radioButtons("type_out", "File output from : ",
                                             choices =  c("Proteome Discoverer" = "PD",
                                                          "MaxQuant" = "MQ",
                                                          "SpectroMine" = "SpM",
                                                          "OpenMS" = "OpMS"), inline = TRUE, selected = "MQ"),

                                tags$hr(),
                                h2("Select your files"),
                                tags$hr(),

                                fluidRow(column(4, uiOutput("file1")),
                                         column(4, conditionalPanel(condition = "input.type_out != 'OpMS'",
                                                                    fileInput("anno",
                                                                              label= h3("Select your annotations file"),
                                                                              accept = c(".xlsx", ".txt", ".csv"))
                                                                    )
                                                ),
                                         column(4, conditionalPanel(condition = "input.type_out == 'MQ'",
                                                                    fileInput("evi",
                                                                              label= h3("Select the evidence file from MaxQuant output"),
                                                                              accept = c(".xlsx", ".txt", ".csv"))
                                                                    )
                                                )
                                ),

                                tags$hr(),
                                h2("Filtering criteria"),
                                tags$hr(),

                                fluidRow(column(2, checkboxInput("Unpep", "Remove peptides that are assigned for more than one proteins",
                                                                 TRUE)),
                                         column(2, checkboxInput("missPSM", "Remove PSM with any missing value within each Run",
                                                                 FALSE)),
                                         column(2, checkboxInput("miss_somePSM", "Remove features that have 1 or 2 measurments within each Run",
                                                                 TRUE)),
                                         column(2, checkboxInput("pr_onepep", "Remove proteins which have only 1 peptide and charge",
                                                                 FALSE)),
                                         column(4, selectInput("summa", "If multiple measurements for certain feature in certain run,
                                                     select the feature with the largest summation or maximal value",
                                                               choices = c("sum", "max"), selected = "sum"))
                                ),
                                conditionalPanel(condition = "input.type_out == 'MQ'",
                                                 checkboxInput("Onbysi", "Remove proteins only identified by site",
                                                               FALSE)),
                                conditionalPanel(condition = "input.type_out == 'PD'",
                                                 checkboxInput("NumPr", "Remove shared peptides by information of Proteins column in PSM sheet",
                                                               TRUE)),
                                conditionalPanel(condition = "input.type_out == 'SpM'",
                                                 fluidRow(column(6, checkboxInput("filt_Qv", "Filter out the intensities that have grater q-value cutoff in EG.Qvalue column.
                                                                                    Those intensities will be replaced by NA and will be considered as censored missing values for imputation purpose",
                                                                                  TRUE)),
                                                          column(4, numericInput("qv_cut", "Choose a q-value cutoff", 0.01, min = 0, max = 1, step = 0.01)))
                                ),
                                uiOutput("prID"),

                                tags$hr(),

                                #check if files are well uploaded
                                conditionalPanel("output.main_fileup & input.type_out == 'OpMS'",
                                                 actionButton("go1", "Submit")),
                                conditionalPanel("output.main_fileup & output.evi_fileup & output.anno_fileup",
                                                 actionButton("go2", "Submit")),
                                conditionalPanel("output.main_fileup & output.anno_fileup & input.type_out != 'MQ'",
                                                 actionButton("go3", "Submit")),

                                tags$hr(),

                                #print the message from the MSstatsTMT's functions
                                textOutput("diagP"), textOutput("diagQ"), textOutput("diagS"), textOutput("diagO"),

                                tags$hr(),
                                h2("Results"),
                                tags$hr(),

                                DT::dataTableOutput("MSform_out"),

                                conditionalPanel(condition = "output.MSform_up",
                                                 downloadButton("save_MSform", "Download table in csv format"))

                                ),

                        tabItem(tabName = "MS_quant",
                                checkboxInput("imp", "Import your own file (need to be in MSstats format)", value = FALSE),

                                conditionalPanel(condition = "input.imp == true",
                                                 fileInput("forqu", "Import your MStats format file")
                                                 ),

                                tags$hr(),
                                h2("Quantification criteria"),
                                tags$hr(),

                                fluidRow(column(4, checkboxInput("glo_norm", "Perform global median normalization on peptide level
                                                                              (equalizing the medians accross all the channels and MS runs)",
                                                                 TRUE)),
                                         column(4, checkboxInput("ref_norm", "Perform reference based normalization between MS runs on protein level data.
                                                                              Need at least one reference channel in each MS run, annotated by 'Norm' in condition column
                                                                              of annotation file.",
                                                                 TRUE)),
                                         column(4, numericInput("maxQu", "Maximum quantile for deciding censored missing value, 0.999 for instance.
                                                                          If NULL, assume missing values are censored",
                                                                 min = 0, max = 1, value = NULL))
                                         ),

                                fluidRow(column(3, checkboxInput("rm_normch", "Remove 'Norm' channels from protein level data", FALSE)),
                                         column(3, checkboxInput("rm_empch", "Remove 'Empty' channels from protein level data", TRUE)),
                                         column(3, selectInput("norm_meth", "Choose a summarization methods to protein-level",
                                                               choices = c("msstats", "MedianPolish", "Median", "LogSum"),
                                                               selected = "msstats")),
                                         conditionalPanel(condition = "input.norm_meth == 'msstats'",
                                                          column(3, checkboxInput("MBimp", "Imputes missing values by Accelated failure model;
                                                                                            (if false, uses minimum value to impute the missing
                                                                                             value for each peptide precursor ion)", TRUE)),)
                                         ),

                                tags$hr(),

                                #check if data are ready (files or data from MSstats format tab)
                                conditionalPanel(condition = "output.MSform_up | output.MSform_fileup",
                                                 actionButton("goqu", "Start quantification"),
                                                 checkboxInput("show_diag", "Hide message from quantification", FALSE)),

                                tags$hr(),
                                conditionalPanel(condition = "input.show_diag == false",
                                                 textOutput("diagQuant")),

                                tags$hr(),
                                h2("Results"),
                                tags$hr(),

                                DT::dataTableOutput("Quanti_output"),

                                conditionalPanel(condition = "output.MSQuanti_up",
                                                 downloadButton("save_MSQuanti", "Download table in csv format"))
                                ),

                        tabItem(tabName = "Group_comp",
                                tabsetPanel(type = "tabs",
                                            tabPanel("Calculation",
                                                     checkboxInput("imp_g", "Import your own quantification file (need to be in MSstats format)", value = FALSE),

                                                     conditionalPanel(condition = "input.imp_g == true",
                                                                      fileInput("forgrp", "Import your MStats format file")
                                                                      ),

                                                     tags$hr(),
                                                     h2("Comparison criteria"),
                                                     tags$hr(),

                                                     fluidRow(column(3, checkboxInput("rm_normch_grp", "Remove 'Norm' channels from protein level data", FALSE)),
                                                              column(3, checkboxInput("rm_empch_grp", "Remove 'Empty' channels from protein level data", TRUE)),
                                                              column(3, checkboxInput("moder", "Moderate t-statistic; if not, uses ordinary t-statistic", TRUE)),
                                                              column(3, selectInput("adj_meth", "Select an adjusted method for multiple comparison",
                                                                                    choices = c("holm", "hochberg", "hommel", "bonferroni",
                                                                                                "BH", "BY", "fdr", "none"),
                                                                                    selected = "BH"))
                                                              ),

                                                     #check id data are ready (files or quantificattion tab output)
                                                     conditionalPanel(condition = "output.MSQuanti_up | output.MSQUANT_fileup",
                                                                      checkboxInput("y_contrast", "Choose your own comparison; if not,
                                                                              will compare all possible pairs between two conditions", FALSE),

                                                                      conditionalPanel(condition = "input.y_contrast == true",
                                                                                       uiOutput("con_matui"),
                                                                                       uiOutput("one_batch"),
                                                                                       uiOutput("sel_batch")),

                                                                      actionButton("gogrco", "Start group comparison"),
                                                                      checkboxInput("show_diagGrC", "Hide message from group comparison", FALSE)
                                                                      ),

                                                     tags$hr(),
                                                     h2("Results"),
                                                     tags$hr(),

                                                     conditionalPanel(condition = "input.show_diagGrC == false",
                                                                      textOutput("diagGrC")),

                                                     tags$hr(),

                                                     DT::dataTableOutput("GrC_output"),

                                                     conditionalPanel(condition = "output.GrC_up",
                                                                      downloadButton("save_GrC", "Download table in csv format"))
                                                     ),

                                            tabPanel("Volcano plot",
                                                     checkboxInput("imp_volc", "Import your own file (need to be in MSstats format)", value = FALSE),

                                                     conditionalPanel(condition = "input.imp_volc == true",
                                                                      fileInput("forVOLC", "Import your file")
                                                     ),

                                                     tags$hr(),
                                                     h2("Criterias"),
                                                     tags$hr(),

                                                     fluidRow(column(3, numericInput("pv_thr", "-log10(p-value) threshold", 1.3,
                                                                                     min = 0, step = 0.1)),
                                                              column(5, checkboxInput("pv_corr", "Use the corrected p-values", value = FALSE),
                                                                     checkboxInput("pv_corrY", "Adjust the p-values", value = FALSE),
                                                                     conditionalPanel(condition = "input.pv_corr | input.pv_corrY",
                                                                                      selectInput("pv_adjmeth", "What was your correction method ?",
                                                                                                  choices = c("holm", "hochberg", "hommel", "bonferroni",
                                                                                                              "BH", "BY", "fdr", "none"),
                                                                                                  selected = "BH")
                                                                                      )
                                                                     ),
                                                              column(3, sliderInput("FC_thr", "log2FC threshold ", value = c(-1,1),
                                                                                     min = -10, max = 10, step = 0.1)),


                                                              ),
                                                     fluidRow(column(3, textInput("volc_title", "Choose a title for your Volcano plot", "Volcano plot")),
                                                              column(3, uiOutput("compar_ui")),
                                                              column(5, checkboxInput("curv", "Draw curvature", FALSE),
                                                                     conditionalPanel(condition = "input.curv",
                                                                                      numericInput("curv_val",
                                                                                                   "Choose a value for the curvature", 0.1,
                                                                                                   min = 0, step = 0.05)))
                                                              ),

                                                     conditionalPanel(condition = "output.volc_fileup",
                                                                      actionButton("VOLCANO_button", "See volcano plot")
                                                                      ),

                                                     tags$hr(),
                                                     h2("Volcano plot"),
                                                     tags$hr(),

                                                     conditionalPanel(condition = "output.volc_fileup",
                                                                      plotOutput("p_vol", height = "800px"),
                                                                      downloadButton("down_volc", "Download volcano plot"))

                                                     )
                                            )
                                ),

                        tabItem(tabName = "Plot",

                                checkboxInput("imp_P", "Import your own quantification file at peptide and protein level (need to be in MSstats format)", value = FALSE),

                                conditionalPanel(condition = "input.imp_P == true",
                                                 fluidRow(column(4, fileInput("forP_pep", "Import your quantitative data at peptide level")),
                                                          column(4, fileInput("forP_pro", "Import your quantitative data at protein level"))
                                                          )
                                                 ),

                                fluidRow(column(3, checkboxInput("Prof_Plot", "Visualize profile plot", TRUE)),
                                         column(3, checkboxInput("QC_Plot", "Visualize QC plot", FALSE)),
                                         column(3, selectizeInput("PROT", "Which protein ?", choices = NULL))
                                         ),

                                conditionalPanel(condition = "input.Prof_Plot == true",
                                                 fluidRow(column(3, checkboxInput("Orig_Plot", "Visualize original profile plot", TRUE)),
                                                          column(3, checkboxInput("summa_Plot", "Visualize summary profile plot", FALSE))
                                                          )
                                                 ),

                                #same thing as before
                                conditionalPanel(condition = "output.MSform_prof_fileup & output.MSQUANT_prof_fileup | output.MSform_up & output.MSQuanti_up",
                                                 fluidRow(conditionalPanel(condition = "input.Prof_Plot == true",
                                                                           column(3, actionButton("prof_button", "See profile plot"))
                                                                           ),
                                                          conditionalPanel(condition = "input.QC_Plot == true",
                                                                           column(3, actionButton("qc_button", "See QC plot")))
                                                          ),
                                                 tags$hr(),
                                                 h2("Plots"),
                                                 tags$hr(),


                                                 conditionalPanel(condition = "input.Prof_Plot == true",
                                                                  withSpinner(plotOutput("p_pr", height = "800px"), type = 6),
                                                                  downloadButton("down_prof", "Download profile plot as png")
                                                                  ),
                                                 conditionalPanel(condition = "input.QC_Plot == true",
                                                                  withSpinner(plotOutput("p_qc"), type = 6),
                                                                  downloadButton("down_qc", "Download QC plot as png")
                                                                  )
                                                 )

                                )
                        )
                      )
                    )

server <- function(input, output, session){
  ### README

  #print the annotation file example
  output$annot_ex_mq <- DT::renderDataTable({
    DT::datatable(MSstatsTMT::annotation.mq,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('annotation file MQ')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })

  output$annot_ex_mine <- DT::renderDataTable({
    DT::datatable(MSstatsTMT::annotation.mine,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('annotation file SpecMine')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })

  output$annot_ex_pd <- DT::renderDataTable({
    DT::datatable(MSstatsTMT::annotation.pd,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('annotation file PD')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })

  output$annot_ex_comp <- DT::renderDataTable({
    DT::datatable(mqpar25_annotation,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('Complex annotation file MQ')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })


  ### ANNOTATION file

  #change the label of the file.input according to the type of file selected
  file1A_name <- reactive({
    if (input$type_outA == "PD"){
      file1_name <- paste("Select the PSM file from Proteome Discoverer output")
    }
    if (input$type_outA == "MQ"){
      file1_name <- paste("Select the proteinGroup file from MaxQuant output")
    }
    if (input$type_outA == "SpM"){
      file1_name <- paste("Select the PSM file from SpectroMine output")
    }
    file1_name
  })
  output$file1A <- renderUI({
    fileInput("mainA",
              label= h3(file1A_name()), accept = c(".xlsx", ".txt", ".csv"))
  })


  #import the file (import_list support almost every type of data file)
  mainA_data <- reactive({
    File <- input$mainA
    if (is.null(File))
      return(NULL)
    import_list(File$datapath, fread = FALSE)[[1]]
  })

  #check if a file is upload
  output$mainA_fileup <- reactive({
    return(!is.null(mainA_data()))
  })
  outputOptions(output, "mainA_fileup", suspendWhenHidden = FALSE)

  evideA <- reactive({
    File <- input$eviA
    if (is.null(File))
      return(NULL)
    import_list(File$datapath, fread = FALSE)[[1]]
  })

  #check if a file is upload
  output$eviA_fileup <- reactive({
    return(!is.null(evideA()))
  })
  outputOptions(output, "eviA_fileup", suspendWhenHidden = FALSE)


  #display the matrix input according your criterias and file
  Condition_matrix <- reactive({
     MAT <- matrix(rep("Norm", input$nchan), ncol = 1)
     colnames(MAT) <-"Condition"
     rownames(MAT) <- paste("channel.", c(1:input$nchan), sep = "")

    MAT
  })
  output$condi_matui <- renderUI({
    matrixInput("condi_mat", value = Condition_matrix(), class = "character", cols = list(names = TRUE),
                rows = list(names = TRUE, editableNames = TRUE))
  })

  Batch_matrix <- reactive({
    MAT <- matrix(rep("batch.1", input$nbatch), ncol = 1)
    colnames(MAT) <-"Batch"

    MAT
  })
  output$batch_matui <- renderUI({
    matrixInput("batch_mat", value = Batch_matrix(), class = "character", cols = list(names = TRUE),
                rows = list(names = FALSE))
  })

  #make the annotation file in the good format
  ANNOTATION <- eventReactive(input$but_anno, {

    if (input$type_outA == "PD"){
      Run <- unique(mainA_data()$Spectrum.File)
    }
    if (input$type_outA == "MQ"){
      Run <- unique(evideA()$Raw.file)
    }
    if (input$type_outA == "SpM"){
      Run <- unique(mainA_data()$R.FileName)
    }

    ANNO <- data.frame(matrix(ncol = 7, nrow = length(Run)*input$nchan, 1)
                       )

    colnames(ANNO) <- c("Run", "Fraction", "TechRepMixture", "Channel",
                        "Condition", "Mixture", "BioReplicate")

    Run <- Run[order(Run)]
    Run <- Run[order(nchar(Run), Run)]
    Run <- rep(Run, each = input$nchan)
    ANNO$Run <- Run

    Fraction <- c(1:input$nfrac)
    Fraction <- rep(Fraction, each = input$nchan)
    Fraction <- rep(Fraction, input$nbatch)
    Fraction <- as.character(Fraction)
    ANNO$Fraction <- Fraction

    channel <- rownames(input$condi_mat)
    channel <- rep(channel, nrow(ANNO)/input$nchan)
    channel <- as.factor(channel)
    ANNO$Channel <- channel

    techrep <- rep(c(1:input$ntech),each = input$nchan)
    techrep <- rep(techrep, nrow(ANNO)/length(techrep))
    techrep <- as.character(techrep)
    ANNO$TechRepMixture <- techrep

    cond <- input$condi_mat[,1]
    cond <-  rep(cond, nrow(ANNO)/input$nchan)
    cond <- as.factor(cond)
    ANNO$Condition <- cond

    batch <- input$batch_mat[,1]
    batch <- rep(batch, each = nrow(ANNO)/length(batch))
    batch <- as.factor(batch)
    ANNO$Mixture <- batch

    ANNO$BioReplicate <- paste(ANNO$Mixture, "_", ANNO$Condition, sep = "")


    ANNO
  })

  #print the annotation tab
  output$ANNOTATION_out <- DT::renderDataTable({
    DT::datatable(ANNOTATION(),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('Your annotation file')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })

  #download the tab
  output$save_ANNOTATION <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_ANNOTATION", ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(ANNOTATION(), file)
    }
  )




  ### MSstats format

  #change the label of the file.input according to the type of file selected
  file1_name <- reactive({
    if (input$type_out == "PD"){
      file1_name <- paste("Select the PSM file from Proteome Discoverer output")
    }
    if (input$type_out == "MQ"){
      file1_name <- paste("Select the proteinGroup file from MaxQuant output")
    }
    if (input$type_out == "SpM"){
      file1_name <- paste("Select the PSM file from SpectroMine output")
    }
    if (input$type_out == "OpMS"){
      file1_name <- paste("Select the MSstatsTMT report from OpenMS")
    }
    file1_name
  })
  output$file1 <- renderUI({
    fileInput("main",
              label= h3(file1_name()), accept = c(".xlsx", ".txt", ".csv"))
  })

  #change the selectinput for the protein ID according to the type of file selected
  output$prID <- renderUI({
    if(input$type_out == "PD"){
      selectInput("prID_PD", "Select the protein ID you want for your protein name",
                  choices = c("Protein.Accessions", "Master.Protein.Accessions"),
                  selected = "Protein.Accessions")
    }
    else if (input$type_out == "MQ"){
      selectInput("prID_MQ", "Select the protein ID you want for your protein name",
                  choices = c("Proteins", "Leading.proteins", "Leading.razor.protein", "Gene.names"),
                  selected = "Proteins")
    }
    else{
      NULL
    }
  })

  ### MSstats format

  #when importing personnal data
  main_data <- reactive({
    File <- input$main
    if (is.null(File))
      return(NULL)
    import_list(File$datapath, fread = FALSE)[[1]]
  })

  #check if a file is upload
  output$main_fileup <- reactive({
    return(!is.null(main_data()))
  })
  outputOptions(output, "main_fileup", suspendWhenHidden = FALSE)

  evide <- reactive({
    File <- input$evi
    if (is.null(File))
      return(NULL)
    import_list(File$datapath, fread = FALSE)[[1]]
  })

  #check if a file is upload
  output$evi_fileup <- reactive({
    return(!is.null(evide()))
  })
  outputOptions(output, "evi_fileup", suspendWhenHidden = FALSE)

  annot <- reactive({
    File <- input$anno
    if (is.null(File))
      return(NULL)
    import_list(File$datapath, fread = FALSE)[[1]]
  })

  #check if a file is upload
  output$anno_fileup <- reactive({
    return(!is.null(annot()))
  })
  outputOptions(output, "anno_fileup", suspendWhenHidden = FALSE)

  #this both criterias cannot be TRUE as the same time
  observe({
    if(input$missPSM){
      updateCheckboxInput(session, "miss_somePSM", value = FALSE)
    }
  })
  observe({
    if(input$miss_somePSM){
      updateCheckboxInput(session, "missPSM", value = FALSE)
    }
  })


  #start the calculation for MSstatsTMT format
  MSstat_data <- eventReactive(input$go1 | input$go2 | input$go3, {
    if (input$type_out == "PD"){
      withCallingHandlers({
        shinyjs::html("diagP", "")
        if (input$summa == "sum"){
          MS <- PDtoMSstatsTMTFormat(main_data(), annot(), which.proteinid = input$prID_PD,
                                     useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                     rmPSM_withfewMea_withinRun =  input$miss_somePSM, useNumProteinsColumn = input$NumPr,
                                     rmProtein_with1Feature = input$pr_onepep)
        }
        else if (input$summa == "max"){
          MS <- PDtoMSstatsTMTFormat(main_data(), annot(), which.proteinid = input$prID_PD,
                                     useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                     rmPSM_withfewMea_withinRun =  input$miss_somePSM, useNumProteinsColumn = input$NumPr,
                                     rmProtein_with1Feature = input$pr_onepep,
                                     summaryforMultipleRows = max)
        }

      },
      #print the message from the function
      message = function(m) {shinyjs::html(id = "diagP", html = paste(m$message, "<br>", sep = ""), add = TRUE)
      }  #use java script to display message when running, without it, only display the last message
      )
    }
    if (input$type_out == "MQ"){
      withCallingHandlers({
        shinyjs::html("diagQ", "")
        if (input$summa == "sum"){
          MS <- MaxQtoMSstatsTMTFormat(evide(), main_data(), annot(), which.proteinid = input$prID_MQ,
                                       rmProt_Only.identified.by.site = input$Onbysi,
                                       useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                       rmPSM_withfewMea_withinRun =  input$miss_somePSM,
                                       rmProtein_with1Feature = input$pr_onepep)
        }
        else if (input$summa == "max"){
          MS <- MaxQtoMSstatsTMTFormat(evide(), main_data(), annot(), which.proteinid = input$prID_MQ,
                                       rmProt_Only.identified.by.site = input$Onbysi,
                                       useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                       rmPSM_withfewMea_withinRun =  input$miss_somePSM,
                                       rmProtein_with1Feature = input$pr_onepep,
                                       summaryforMultipleRows = max)
        }

        },
        #print the message from the function
        message = function(m) {shinyjs::html(id = "diagQ", html = paste(m$message, "<br>", sep = ""), add = TRUE)
          }  #use java script to display message when running, without it, only display the last message
        )
      }
    if (input$type_out == "SpM"){
      withCallingHandlers({
        shinyjs::html("diagS", "")
        if (input$summa == "sum"){
          MS <- SpectroMinetoMSstatsTMTFormat(main_data(), annot(),
                                              useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                              rmPSM_withfewMea_withinRun =  input$miss_somePSM, filter_with_Qvalue = input$filt_Qv,
                                              qvalue_cutoff = input$qv_cut,
                                              rmProtein_with1Feature = input$pr_onepep)
        }
        else if (input$summa == "max"){
          MS <- SpectroMinetoMSstatsTMTFormat(main_data(), annot(),
                                              useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                              rmPSM_withfewMea_withinRun =  input$miss_somePSM, filter_with_Qvalue = input$filt_Qv,
                                              qvalue_cutoff = input$qv_cut,
                                              rmProtein_with1Feature = input$pr_onepep,
                                              summaryforMultipleRows = max)
        }

      },
      #print the message from the function
      message = function(m) {shinyjs::html(id = "diagS", html = paste(m$message, "<br>", sep = ""), add = TRUE)
      }  #use java script to display message when running, without it, only display the last message
      )
    }
    if (input$type_out == "OpMS"){
      withCallingHandlers({
        shinyjs::html("diagO", "")
        if (input$summa == "sum"){
          MS <- OpenMStoMSstatsTMTFormat(main_data(),
                                         useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                         rmPSM_withfewMea_withinRun =  input$miss_somePSM,
                                         rmProtein_with1Feature = input$pr_onepep)
        }
        else if (input$summa == "max"){
          MS <- OpenMStoMSstatsTMTFormat(main_data(),
                                         useUniquePeptide = input$Unpep, rmPSM_withMissing_withinRun = input$missPSM,
                                         rmPSM_withfewMea_withinRun =  input$miss_somePSM,
                                         rmProtein_with1Feature = input$pr_onepep,
                                         summaryforMultipleRows = max)
        }

      },
      #print the message from the function
      message = function(m) {shinyjs::html(id = "diagO", html = paste(m$message, "<br>", sep = ""), add = TRUE)
      }  #use java script to display message when running, without it, only display the last message
      )
    }

    MS
  }, ignoreInit = TRUE)

  #check if data are ready
  output$MSform_up <- reactive({
    return(!is.null(MSstat_data()))
  })
  outputOptions(output, "MSform_up", suspendWhenHidden = FALSE)

  #print the tab
  output$MSform_out <- DT::renderDataTable({
    DT::datatable(MSstat_data(),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('MSstats file for quantification')
                    ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
                  )
  })

  #download the tab
  output$save_MSform <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_MSstats_format", ".csv", sep = "")
    },
    content = function(file) {
      write.table(MSstat_data(), file, row.names = FALSE, quote = FALSE, sep=",")
    }
  )


  ### QUANTIFICATION

  #import your own file or use the data calculated by the app
  quant_data <- reactive({
    if (input$imp){
      File <- input$forqu
      if (is.null(File))
        return(NULL)
      QU <- import_list(File$datapath)[[1]]
    }
    else{
      QU <- MSstat_data()
    }
    QU
  })

  #check if a file is upload
  output$MSform_fileup <- reactive({
    if(input$imp){
      return(!is.null(quant_data()))
    }
    })
  outputOptions(output, "MSform_fileup", suspendWhenHidden = FALSE)

  #start quantification
  quant_final <- eventReactive(input$goqu, {
    mQuant <- input$maxQu
    if (is.na(mQuant)){
      mQuant <- NULL
    }
    withCallingHandlers({
      shinyjs::html("diagQuant", "")
      proteinSummarization(quant_data(), method = input$norm_meth,
                           global_norm = input$glo_norm, reference_norm = input$ref_norm,
                           remove_norm_channel = input$rm_normch, remove_empty_channel = input$rm_empch,
                           MBimpute = input$MBimp, maxQuantileforCensored = mQuant)
      },
    message = function(m) {
      m_ <- m$message
      if(str_length(m_) > 120)   #some message are very long so not display the whole message if it's the case
        m_ <- paste(str_sub(m_, 1, 120), "...", sep = "")
      shinyjs::html(id = "diagQuant", html = paste(m_, "<br>", sep = ""), add = FALSE)
    }  #use java script to display message when running, without it, only display the last message
    )

  }, ignoreInit = TRUE)


  #check if data are ready
  output$MSQuanti_up <- reactive({
    return(!is.null(quant_final()))
  })
  outputOptions(output, "MSQuanti_up", suspendWhenHidden = FALSE)

  #print the tab
  output$Quanti_output <- DT::renderDataTable({
    DT::datatable(quant_final(),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('MSstats results')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })

  #download the tab
  output$save_MSQuanti <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_MSstats_results", ".csv", sep = "")
    },
    content = function(file) {
      write.table(quant_final(), file, row.names = FALSE, quote = FALSE, sep = ",")
    }
  )


  ### Group comparison

  #import your own file or the results from the app
  compar_data <- reactive({
    if (input$imp_g){
      File <- input$forgrp
      if (is.null(File))
        return(NULL)
      GC <- import_list(File$datapath)[[1]]
    }
    else{
      GC <- quant_final()
    }
    GC
  })

  #check if a file is upload
  output$MSQUANT_fileup <- reactive({
    if(input$imp_g){
      return(!is.null(compar_data()))
    }
  })
  outputOptions(output, "MSQUANT_fileup", suspendWhenHidden = FALSE)


  #print the matrix input for comparison according to the data
  contr_mat <- reactive({
    if (is.null(compar_data())){
      MAT <- NULL
    }
    else{
      MAT <- matrix(rep(1, length(unique(compar_data()$Condition))), nrow = 1)
      colnames(MAT) <- unique(compar_data()$Condition)
      rownames(MAT) <- "Comparison"
    }
    MAT
  })
  output$con_matui <- renderUI({
    matrixInput("con_mat", value = contr_mat(), class = "numeric", cols = list(names = TRUE),
                rows = list(names = TRUE, editableNames = TRUE))
  })

  #if more than one batch, allow to choose specific batch
  output$one_batch <- renderUI({
    if (length(unique(compar_data()$Mixture)) > 1){
      checkboxInput("o_batch", "Perfom group comparison on one batch ?", FALSE)
    }
    else{
      NULL
    }
  })

  #if more than one batch and want to see specific one, print the selectinput according tot the data
  output$sel_batch <- renderUI({
    if (length(unique(compar_data()$Mixture)) > 1){
      if (input$o_batch){
        selectInput("batch", "Select the batch you want to perform the group comparison on",
                    choices = as.character(unique(compar_data()$Mixture)))
      }
      else{
        NULL
      }
    }
    else{
      NULL
    }
  })

  #start the Group comparison
  GrC_final <- eventReactive(input$gogrco, {
    compar_data2 <- compar_data()
    if (length(unique(compar_data()$Mixture)) > 1  & input$y_contrast){
      if (input$o_batch){
        compar_data2 <- compar_data2[which(compar_data2$Mixture == input$batch),]  #filter the batch selected if so
      }
    }
    withCallingHandlers({
      shinyjs::html("diagGrC", "")
      if (input$y_contrast){
        GC <- groupComparisonTMT(compar_data2, moderated = input$moder,
                                 adj.method = input$adj_meth, remove_norm_channel = input$rm_normch_grp,
                                 remove_empty_channel = input$rm_empch_grp, contrast.matrix = input$con_mat)
      }
      else{
        GC <- groupComparisonTMT(compar_data2, moderated = input$moder,
                                 adj.method = input$adj_meth, remove_norm_channel = input$rm_normch_grp,
                                 remove_empty_channel = input$rm_empch_grp)
      }
        },
      message = function(m) {shinyjs::html(id = "diagGrC", html = paste(m$message, "<br>", sep = ""), add = FALSE)
      }  #use java script to display message when running, without it, only display the last message
      )

    GC

  }, ignoreInit = TRUE)

  #check if data are ready
  output$GrC_up <- reactive({
    return(!is.null(GrC_final()))
  })
  outputOptions(output, "GrC_up", suspendWhenHidden = FALSE)

  #print the tab
  output$GrC_output <- DT::renderDataTable({
    DT::datatable(GrC_final(),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    htmltools::strong('MSstats results from group comparison')
                  ),
                  rownames = FALSE,
                  options = list(lengthMenu = c(10,20,30), pageLength = 10)
    )
  })

  #download the tab
  output$save_GrC <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_GrpCompare_results", ".csv", sep = "")
    },
    content = function(file) {
      write.table(GrC_final(), file, row.names = FALSE, quote = FALSE, sep = ",")
    }
  )

  ## VOLCANO

  #import your own file or use results from the app
  VOLC_data <- reactive({
    if (input$imp_volc){
      File <- input$forVOLC
      if (is.null(File))
        return(NULL)
      QU <- import_list(File$datapath)[[1]]
    }
    else{
      QU <- GrC_final()
    }
    QU
  })

  #check if a file is upload
  output$volc_fileup <- reactive({
    return(!is.null(VOLC_data()))
  })
  outputOptions(output, "volc_fileup", suspendWhenHidden = FALSE)

  #update value of argument if some are selected to 'satisfy' argument of the function
  observe({
    if(input$pv_corr){
      updateSelectInput(session, "pv_adjmeth", label = "What was your correction method ?")
      updateCheckboxInput(session, "pv_corrY", value = FALSE)
    }
    else{
      updateSelectInput(session, "pv_adjmeth", selected = "none")
    }
  })
  observe({
    if(input$pv_corrY){
      updateSelectInput(session, "pv_adjmeth", label = "Chooose a correction method")
      updateCheckboxInput(session, "pv_corr", value = FALSE)
    }
    else{
      updateSelectInput(session, "pv_adjmeth", selected = "none")
    }
  })

  #update the selectInput for condition according to the data
  output$compar_ui <- renderUI({
    selectInput("compar_val", "Choose the two condition you want to compare", choices = unique(VOLC_data()$Label))
  })

  #plot the volcano
  Volcano_Plot <- reactive({
    plo_volc(VOLC_data(), lim_pv = 10**-input$pv_thr, lim_dif = input$FC_thr, correction = input$pv_corr, perf_corr = input$pv_corrY,
             comp = input$compar_val, curve = input$curv, curvature = input$curv_val, tit = input$volc_title,
              ytit = "-log10(p-value)", your_corr = input$pv_adjmeth)

  })

  #use a button to display the plot
  volc_plot <- reactiveValues(
    ch = NULL
  )

  observeEvent(input$VOLCANO_button, {
    volc_plot$ch <- Volcano_Plot()
  },
  ignoreInit = TRUE, ignoreNULL = FALSE
  )

  output$p_vol <- renderPlot({
    volc_plot$ch
  })

  ##download the plot
  output$down_volc <- downloadHandler(
    filename = function() {
      paste("Volcano_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file){
      ggsave(file, volc_plot$ch, device = "png", width = 25, height = 25, units = "cm")
    }
  )



  ### PROFILE PLOT

  #import your own or use results from the app (MSstatsTMT Format)
  Prof_MS_data <- reactive({
    if (input$imp_P){
      File <- input$forP_pep
      if (is.null(File))
        return(NULL)
      QU <- import_list(File$datapath)[[1]]
    }
    else{
      QU <- MSstat_data()
    }
    QU
  })

  #check if a file is upload
  output$MSform_prof_fileup <- reactive({
      return(!is.null(Prof_MS_data()))
  })
  outputOptions(output, "MSform_prof_fileup", suspendWhenHidden = FALSE)

  #import your own file or use results from the app (quantification)
  Prof_qu_data <- reactive({
    if (input$imp_P){
      File <- input$forP_pro
      if (is.null(File))
        return(NULL)
      QU <- import_list(File$datapath)[[1]]
    }
    else{
      QU <- quant_final()
    }
    QU
  })

  #check if a file is upload
  output$MSQUANT_prof_fileup <- reactive({
      return(!is.null(Prof_qu_data()))
  })
  outputOptions(output, "MSQUANT_prof_fileup", suspendWhenHidden = FALSE)



  #update the proteins you can selected according to the data
  observe({
    if (is.null(Prof_MS_data()) | is_empty(Prof_MS_data())){
      updateSelectizeInput(session, "PROT", choices = unique(Prof_qu_data()$Protein))
    }
    if (is.null(Prof_qu_data()) | is_empty(Prof_qu_data())){
      updateSelectizeInput(session, "PROT", choices = unique(Prof_MS_data()$ProteinName))
    }
  })


  #plot the protein profile or QC plot, also use a button to display the plots
  Profile_Plot <- reactive({
    if (input$Prof_Plot){
      PlotsTMT_profile(Prof_MS_data(), Prof_qu_data(),
                       which.Protein = input$PROT, originalPlot = input$Orig_Plot,
                       summaryPlot = input$summa_Plot)
    }
    else{
      NULL
    }
    })

  prof_plot <- reactiveValues(
    ch = NULL
  )

  observeEvent(input$prof_button, {
    prof_plot$ch <- Profile_Plot()
    },
    ignoreInit = TRUE, ignoreNULL = FALSE
    )

  output$p_pr <- renderPlot({
    prof_plot$ch
  })



  QC_PLOT <- reactive({
    if (input$QC_Plot){
      PlotsTMT_QC(Prof_MS_data(), Prof_qu_data(),
                  which.Protein = input$PROT)
    }
    else{
      NULL
    }
  })

  qc_plot <- reactiveValues(
    ch = NULL
  )

  observeEvent(input$qc_button, {
    qc_plot$ch <- QC_PLOT()
  },
  ignoreInit = TRUE, ignoreNULL = FALSE
  )

  output$p_qc <- renderPlot({
    qc_plot$ch
  })


  #download the plots
  output$down_prof <- downloadHandler(
    filename = function() {
      paste("profile_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file){
      ggsave(file, prof_plot$ch, device = "png", width = 30, height = 22, units = "cm")
    }
  )

  output$down_qc <- downloadHandler(
    filename = function() {
      paste("QC_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file){
      ggsave(file, qc_plot$ch, device = "png", width = 25, height = 18, units = "cm")
    }
  )

}

shinyApp(ui, server)

