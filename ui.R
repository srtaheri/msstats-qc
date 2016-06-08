library(shiny)
library(plotly)
library(markdown)
shinyUI(fluidPage(
  titlePanel(p(strong("MSstatsQC"),align = "center",style="color:#444444;",style="font-size:150%;",style="font-family:inherit;")),
  navbarPage(h4("System suitability monitoring tools for quantitative mass spectrometry based proteomic experiments"),
              tabPanel("Home", theme = "bootstrap.css",
                         tags$img(src='logo.png', height=220, width=220, style = "float: right"),
                         p("MSstatsQC is an open-source web-based software which provides longitudinal
                           system suitability monitoring tools (control charts) for SRM based proteomic experiments."),
                        h5("Statistical functionalities"),
                        p("MSstatsQC uses control charts to monitor the instrument performance by tracking system
                           suitability metrics including total peak area, retention time and full width at half maximum (FWHM) and peak assymetry."),
                         p("This framework includes simultaneous monitoring tools for mean and dispersion of suitability metrics and presents
                           alternative methods of monitoring such as time weighted control charts to ensure that various types
                           of process disturbances are detected effectively. Simultaneous control charts used in this framework
                           can be classified into two groups: individual-moving range (XmR) control charts and mean and dispersion
                           cumulative sum (CUSUM) control charts. To successfully identify the time of change, change point analysis
                           is also included in this framework. Experiment specific control limits are provided with the control 
                           charts to distinguish between random noise and systematic error."),
                       h5("Using MSstatsQC"),
                       p("The steps for generating results are as follows:"),
                         ("1. Import your QC data "),
                         br(),
                         ("2.	Determine the guide set to estimate metric mean and variance "),
                         br(),
                         ("3.	Select specific precursor(s) or select all"),
                         br(),
                         ("4. Run and generate control charts"),
                         br(),
                         ("5.	Check with change point analysis for better reasoning"),
                         br(),
                         ("6.	Navigate results and download them for your QC reports"),
                         br(),
                         br(),
                         tags$img(src='home.png', height=200, width=500, style = "float: center"),
                         br(),
                         br(),
                         h5 ("Project Team: "),
                         h5("Eralp Dogu,",span("e.dogu@neu.edu",style = "color:blue")),
                         h5("Sara Taheri,",span("mohammadtaheri.s@husky.neu.edu",style = "color:blue")),
                         h5("Olga Vitek,",span("o.vitek@neu.edu",style = "color:blue")),
                         br(),
                         br(),
                         ("Olga Vitek Lab"),
                         br(),
                         ("College of Science"),
                         br(),
                         ("College of Computer and Information Science"),
                         br(),
                         ("360 Huntington Ave"),
                         br(),
                         ("Boston, Massachusetts 02115"),
                         br(),
                         br(),
                         br()

              ), # end mainPanel

#              ),
#############origin/master
              tabPanel("Data Import", 
                       sidebarLayout(

                        sidebarPanel(
                           p("Please upload your data (Comma-separated (*.csv) 8 column QC file format)"),
                           p("To see an acceptable sample data, look at", strong("Help"),"tab"),
                           fileInput("filein", "Upload file"),
                           p("Please select a guide set"),
                           numericInput("L","Lower bound of guide set",value = 1, min = 1, step = 1),
                           numericInput("U","Upper bound of guide set", value = 5, min = 2, step = 1),
                           p("Please select a precursor or select all"),
                           uiOutput("pepSelect"),
                           #p("If you have uploaded your data set, selected the guidset and chosen your precursor type
                             #, click on this button to see the plots"),
                           p("if you have uploaded your data set, click on this button to view it."),
                           #actionButton("act_button", "click to see your data set"),
                           helpText("please select the columns of your data that you need to see"),
                           uiOutput("prodata_column_select"),
                           tags$style("body{background-color:linen; color:black}"),
                           p("If you want to run", strong("MSstatsQC"), "with sample data file, please click this button"),
                           actionButton("act-button", "Run with sample data")
                           ),
                        mainPanel(
                           
                          tabPanel("Data",  
                                      DT::dataTableOutput('prodata_table'))
                        ), 
                           position = "left")
                       ),
              
              tabPanel("Metric Summary",
                       tabsetPanel(
                         tabPanel("Boxplot",
                                  plotlyOutput("box_plot", height = 2000)
                         ),
                         tabPanel("Scatterplot",
                                  selectInput("metric_precursor", "Choose the metric",
                                              choices = c("Retention Time","Total Area","FWHM","Peak Assymetry")),
                                  plotOutput("scatter_plot")
                         )
                       )),
              navbarMenu("Control Charts",
                         tabPanel("XmR",
                                  tabsetPanel(
                                    tabPanel("Retention Time",
                                             plotlyOutput("RT_ZMR")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              ),
                                    tabPanel("Total Peak Area", 
                                             plotlyOutput("TA_ZMR")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              ),
                                    tabPanel("Full Width at Half Maximum (FWHM)", 
                                             plotlyOutput("Max_ZMR")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              ),
                                    tabPanel("Peak Assymetry", 
                                             plotlyOutput("PA_ZMR")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              )
                                    #,
                                    # tabPanel("Mass Accuracy"
                                    #          , textOutput("MA_ZMR_txt")
                                    #          ,plotlyOutput("MA_ZMR")
                                    #          ,tags$head(tags$style(type="text/css", "
                                    #                                #loadmessage 
                                    #                                {
                                    #                                position: fixed;
                                    #                                top: 0px;
                                    #                                left: 0px;
                                    #                                width: 100%;
                                    #                                padding: 5px 0px 5px 0px;
                                    #                                text-align: center;
                                    #                                font-weight: bold;
                                    #                                font-size: 100%;
                                    #                                color: #000000;
                                    #                                background-color: #CCFF66;
                                    #                                z-index: 105;
                                    #                                }
                                    #                                ")),
                                    #          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    #                           tags$div("It may take a while to load the plots, please wait...
                                    #                                    ",id="loadmessage"))
                                    #                           )
                                             )), # End "XMR" tabPanel and it's tabsetPanel
                         
                         tabPanel("CUSUM",
                                  tabsetPanel(
                                    tabPanel("Retention Time", 
                                             #plotOutput("RT_CUSUM"
                                             plotlyOutput("RT_CUSUM"
                                             )
                                             #,verbatimTextOutput("RT_CUSUM_info")  # for interactive plots
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              ),
                                    tabPanel("Total Peak Area", 
                                             #plotOutput("TA_CUSUM")
                                             plotlyOutput("TA_CUSUM")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              ),
                                    tabPanel("Full Width at Half Maximum (FWHM)", 
                                             #plotOutput("Max_CUSUM")
                                             plotlyOutput("Max_CUSUM")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              ),
                                    tabPanel("Peak Assymetry", 
                                             #plotOutput("PA_CUSUM")
                                             plotlyOutput("PA_CUSUM")
                                             ,tags$head(tags$style(type="text/css", "
                                                                   #loadmessage 
                                                                   {
                                                                   position: fixed;
                                                                   top: 0px;
                                                                   left: 0px;
                                                                   width: 100%;
                                                                   padding: 5px 0px 5px 0px;
                                                                   text-align: center;
                                                                   font-weight: bold;
                                                                   font-size: 100%;
                                                                   color: #000000;
                                                                   background-color: #CCFF66;
                                                                   z-index: 105;
                                                                   }
                                                                   ")),
                                             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                              tags$div("It may take a while to load the plots, please wait...
                                                                       ",id="loadmessage"))
                                                              )
                                    #,
                                    # tabPanel("Mass Accuracy" 
                                    #          , textOutput("MA_CUSUM_txt")
                                    #          ,plotlyOutput("MA_CUSUM")
                                    #          ,tags$head(tags$style(type="text/css", "
                                    #                                #loadmessage 
                                    #                                {
                                    #                                position: fixed;
                                    #                                top: 0px;
                                    #                                left: 0px;
                                    #                                width: 100%;
                                    #                                padding: 5px 0px 5px 0px;
                                    #                                text-align: center;
                                    #                                font-weight: bold;
                                    #                                font-size: 100%;
                                    #                                color: #000000;
                                    #                                background-color: #CCFF66;
                                    #                                z-index: 105;
                                    #                                }
                                    #                                ")),
                                    #          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    #                           tags$div("It may take a while to load the plots, please wait...
                                    #                                    ",id="loadmessage"))
                                    #                           )
                                             ))
                         #, # End tabpanle "CUSUM" and tabsetPanel of it
                         #tabPanel("EWMA", textOutput("EWMA_txt")),
                         #tabPanel("Short run SPC", textOutput("Short_run_SPC_txt")),
                       
                        # tabPanel("Multivariate Control Charts"
                        #          , textOutput("Multivariate_Control_Charts_txt"),
                        #          tabsetPanel(
                        #            tabPanel("Retention Time",
                        #                     plotlyOutput("RT_Multi")
                        #                     ,tags$head(tags$style(type="text/css", "
                        #                                           #loadmessage 
                        #                                           {
                        #                                           position: fixed;
                        #                                           top: 0px;
                        #                                           left: 0px;
                        #                                           width: 100%;
                        #                                           padding: 5px 0px 5px 0px;
                        #                                           text-align: center;
                        #                                           font-weight: bold;
                        #                                           font-size: 100%;
                        #                                           color: #000000;
                        #                                           background-color: #CCFF66;
                        #                                           z-index: 105;
                        #                                           }
                        #                                           ")),
                        #                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      ),
                        #            tabPanel("Total Peak Area", 
                        #                     plotlyOutput("TA_Multi")
                        #                     ,tags$head(tags$style(type="text/css", "
                        #                                           #loadmessage 
                        #                                           {
                        #                                           position: fixed;
                        #                                           top: 0px;
                        #                                           left: 0px;
                        #                                           width: 100%;
                        #                                           padding: 5px 0px 5px 0px;
                        #                                           text-align: center;
                        #                                           font-weight: bold;
                        #                                           font-size: 100%;
                        #                                           color: #000000;
                        #                                           background-color: #CCFF66;
                        #                                           z-index: 105;
                        #                                           }
                        #                                           ")),
                        #                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      ),
                        #            tabPanel("Full Width at Half Maximum (FWHM)", 
                        #                     plotlyOutput("Max_Multi")
                        #                     ,tags$head(tags$style(type="text/css", "
                        #                                           #loadmessage 
                        #                                           {
                        #                                           position: fixed;
                        #                                           top: 0px;
                        #                                           left: 0px;
                        #                                           width: 100%;
                        #                                           padding: 5px 0px 5px 0px;
                        #                                           text-align: center;
                        #                                           font-weight: bold;
                        #                                           font-size: 100%;
                        #                                           color: #000000;
                        #                                           background-color: #CCFF66;
                        #                                           z-index: 105;
                        #                                           }
                        #                                           ")),
                        #                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      ),
                        #            tabPanel("Peak Assymetry", 
                        #                     plotlyOutput("PA_Multi")
                        #                     ,tags$head(tags$style(type="text/css", "
                        #                                           #loadmessage 
                        #                                           {
                        #                                           position: fixed;
                        #                                           top: 0px;
                        #                                           left: 0px;
                        #                                           width: 100%;
                        #                                           padding: 5px 0px 5px 0px;
                        #                                           text-align: center;
                        #                                           font-weight: bold;
                        #                                           font-size: 100%;
                        #                                           color: #000000;
                        #                                           background-color: #CCFF66;
                        #                                           z-index: 105;
                        #                                           }
                        #                                           ")),
                        #                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      )
                                   #,
                                   # tabPanel("Mass Accuracy", 
                                   #          plotlyOutput("MA_Multi")
                                   #          ,tags$head(tags$style(type="text/css", "
                                   #                                #loadmessage 
                                   #                                {
                                   #                                position: fixed;
                                   #                                top: 0px;
                                   #                                left: 0px;
                                   #                                width: 100%;
                                   #                                padding: 5px 0px 5px 0px;
                                   #                                text-align: center;
                                   #                                font-weight: bold;
                                   #                                font-size: 100%;
                                   #                                color: #000000;
                                   #                                background-color: #CCFF66;
                                   #                                z-index: 105;
                                   #                                }
                                   #                                ")),
                                   #          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                   #                           tags$div("It may take a while to load the plots, please wait...
                                   #                                    ",id="loadmessage"))
                                   #                           )
                                    #        ))
), # End "Multi" tabPanel and it's tabsetPanel
              
              tabPanel("Change Point Analysis",
                       tabsetPanel(
                         tabPanel("Retention Time"
                                  #,downloadButton(outputId = "down_CP_RT", label = "Download the plots")
                                  , plotlyOutput("RT_CP")
                                  ,tags$head(tags$style(type="text/css", "
                                                        #loadmessage 
                                                        {
                                                        position: fixed;
                                                        top: 0px;
                                                        left: 0px;
                                                        width: 100%;
                                                        padding: 5px 0px 5px 0px;
                                                        text-align: center;
                                                        font-weight: bold;
                                                        font-size: 100%;
                                                        color: #000000;
                                                        background-color: #CCFF66;
                                                        z-index: 105;
                                                        }
                                                        ")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...
                                                            ",id="loadmessage"))
                                                   ),
                         tabPanel("Total Peak Area"
                                  #,downloadButton(outputId = "down_CP_TA", label = "Download the plots")
                                  , plotlyOutput("TA_CP")
                                  ,tags$head(tags$style(type="text/css", "
                                                        #loadmessage 
                                                        {
                                                        position: fixed;
                                                        top: 0px;
                                                        left: 0px;
                                                        width: 100%;
                                                        padding: 5px 0px 5px 0px;
                                                        text-align: center;
                                                        font-weight: bold;
                                                        font-size: 100%;
                                                        color: #000000;
                                                        background-color: #CCFF66;
                                                        z-index: 105;
                                                        }
                                                        ")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...
                                                            ",id="loadmessage"))
                                                   ),
                         tabPanel("Full Width at Half Maximum (FWHM)"
                                  #,downloadButton(outputId = "down_CP_FWHM", label = "Download the plots")
                                  , plotlyOutput("Max_CP")
                                  ,tags$head(tags$style(type="text/css", "
                                                        #loadmessage 
                                                        {
                                                        position: fixed;
                                                        top: 0px;
                                                        left: 0px;
                                                        width: 100%;
                                                        padding: 5px 0px 5px 0px;
                                                        text-align: center;
                                                        font-weight: bold;
                                                        font-size: 100%;
                                                        color: #000000;
                                                        background-color: #CCFF66;
                                                        z-index: 105;
                                                        }
                                                        ")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...
                                                            ",id="loadmessage"))
                                                   ),
                         tabPanel("Peak Assymetry"
                                  #,column(5,radioButtons("file_type","select the file type to download the plots",choices = list("png","pdf"), selected = "pdf"))
                                  #,column(3,downloadButton(outputId = "down_CP_PA", label = "Download the plots"))
                                  ,plotlyOutput("PA_CP")
                                  ,tags$head(tags$style(type="text/css", "
                                                        #loadmessage 
                                                        {
                                                        position: fixed;
                                                        top: 0px;
                                                        left: 0px;
                                                        width: 100%;
                                                        padding: 5px 0px 5px 0px;
                                                        text-align: center;
                                                        font-weight: bold;
                                                        font-size: 100%;
                                                        color: #000000;
                                                        background-color: #CCFF66;
                                                        z-index: 105;
                                                        }
                                                        ")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...
                                                            ",id="loadmessage"))
                                                   )
                         #,
                         # tabPanel("Mass Accuracy"
                         #          , textOutput("CP_MA_txt")
                         #          , plotlyOutput("MA_CP")
                         #          ,tags$head(tags$style(type="text/css", "
                         #                                #loadmessage 
                         #                                {
                         #                                position: fixed;
                         #                                top: 0px;
                         #                                left: 0px;
                         #                                width: 100%;
                         #                                padding: 5px 0px 5px 0px;
                         #                                text-align: center;
                         #                                font-weight: bold;
                         #                                font-size: 100%;
                         #                                color: #000000;
                         #                                background-color: #CCFF66;
                         #                                z-index: 105;
                         #                                }
                         #                                ")),
                         #          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         #                           tags$div("It may take a while to load the plots, please wait...
                         #                                    ",id="loadmessage"))
                         #                           )
                         #          )), # end "Change Point Analysis" tabPanel and it's tabsetPanel
              
              # tabPanel("Capability Analysis",
              #          tabsetPanel(
              #            tabPanel("Retention Time"
              #                     , textOutput("CA_RT_txt")
              #                     , plotlyOutput("RT_CA")
              #                     ,tags$head(tags$style(type="text/css", "
              #                                           #loadmessage 
              #                                           {
              #                                           position: fixed;
              #                                           top: 0px;
              #                                           left: 0px;
              #                                           width: 100%;
              #                                           padding: 5px 0px 5px 0px;
              #                                           text-align: center;
              #                                           font-weight: bold;
              #                                           font-size: 100%;
              #                                           color: #000000;
              #                                           background-color: #CCFF66;
              #                                           z-index: 105;
              #                                           }
              #                                           ")),
              #                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
              #                                      tags$div("It may take a while to load the plots, please wait...
              #                                               ",id="loadmessage"))
              #                                      ),
                         # tabPanel("Mass Accuracy"
                         #          , textOutput("CA_MA_txt")
                         #          , plotlyOutput("MA_CA")
                         #          ,tags$head(tags$style(type="text/css", "
                         #                                #loadmessage 
                         #                                {
                         #                                position: fixed;
                         #                                top: 0px;
                         #                                left: 0px;
                         #                                width: 100%;
                         #                                padding: 5px 0px 5px 0px;
                         #                                text-align: center;
                         #                                font-weight: bold;
                         #                                font-size: 100%;
                         #                                color: #000000;
                         #                                background-color: #CCFF66;
                         #                                z-index: 105;
                         #                                }
                         #                                ")),
                         #          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         #                           tags$div("It may take a while to load the plots, please wait...
                         #                                    ",id="loadmessage"))
                         #                           )
                                  )), # end "Capability Analysis" tabPanel and it's tabsetPanel
              #tabPanel("Overall QC Performance", textOutput("OverallQC_txt")),
              tabPanel("Help",
                       tabsetPanel(
                         tabPanel("Metrics"
                                  ,h5(strong("Retention Time")),
                                  p("Retention time is the time it takes a solute to travel through the column. The retention time is assigned to 
                                    the corresponding solute peak. The retention time is a measure of the amount of time a solute spends in a column. 
                                    It is the sum of the time spent in the stationary phase and the mobile phase."
                                    ,a("visit for more info",href="http://www.britannica.com/science/retention-time")),
                                  
                                  h5(strong("Total Peak Area")),
                                  p("Total Peak Area is the sum of all integrated signals for a certain peptide."),
                                  
                                  h5(strong("Full Width at Half Maximum (FWHM)")),
                                  p("Full width at half maximum 'FWHM' is an expression of the extent of a 
                                    function given by the difference between the two extreme values 
                                    of the independent variable at which the dependent variable is equal 
                                    to half of its maximum value." 
                                    ,a("visit for more info",href="https://en.wikipedia.org/wiki/Full_width_at_half_maximum")),
                                  
                                  h5(strong("Peak Assymetry")),
                                  p("Peak Assymetry is a measure of symetry for a peak. Calculated by taking 2*a/(a+b). Optimal value is around 1 for a Gaussian peak.")
                                  ),
                         
                         tabPanel("Plots"
                                  
                                  ,h5(strong("XmR control charts")),
                                  h5("Can detect=Large shifts and spikes in mean and dispersion of suitability metric."),
                                  h5("Statistics used to construct=using the sequential differences between successive values as a measure 
                                    of dispersion."),
                                  p("When data are collected as individual observations, you cannot calculate the standard 
                                    deviation for each subgroup. The moving range is an alternative way to calculate process 
                                    variation by computing the ranges of two or more consecutive observations.XmR is a cobination 
                                    of a chart for individual observations and 
                                    a chart for moving ranges. I-MR chart consists of 2 charts; 
                                    Individuals (X) chart and Moving Range (mR) chart. 
                                    ", 
                                    a("visit for more info",href="https://en.wikipedia.org/wiki/Shewhart_individuals_control_chart"))
                                  
                                  
                                  ,h5(strong("CUSUMm and CUSUMv control charts")),
                                  h5("Can detect=Large shifts and spikes in mean and dispersion of suitability metric."),
                                  h5("Statistics used to construct=using the sequential differences between successive values as a measure 
                                     of dispersion."),
                                    p("A CUSUM chart is a time-weighted control chart that displays the cumulative sums 
                                    'CUSUMs' of the deviations of each sample value from the target value. Because it is cumulative, 
                                    even minor drifting in the process mean will lead to steadily 
                                    increasing or decreasing cumulative deviation values. In that case, you are often primarily interested 
                                    in detecting small shifts in the product quality (for example, gradual deterioration of quality due to 
                                    machine wear)", 
                                    a("visit for more info",href="https://en.wikipedia.org/wiki/CUSUM"))
                                
                                  
                                  ,h5(strong("Change Point Analysis")),
                                  h5("Can detect=Large shifts and spikes in mean and dispersion of suitability metric."),
                                  h5("Statistics used to construct=using the sequential differences between successive values as a measure 
                                     of dispersion."),
                                  p("A change in the process parameters triggers a control chart to generate an out of
                                    control signal. The QC sample at which the signal is issued is considered as the 
                                    stopping time and after the signal search for an assignable cause is recommended."
                                  , a("visit for more info",href="http://www.eng.fsu.edu/~pigna/pdf/"))
                         ),
                         
                         tabPanel("Documentation",
                                  h5(strong("MSstatsQC webpage")),
                                  p("Source codes, related documents and user manual can be found via our MSstats website"
                                    , a("visit for more info",href="http://www.msstats.org/msstatsqc"))
                                  
                                  ,h5(strong("MSstatsQC Github")),
                                  p("Latest documantation is also available via our Github page"
                                    , a("visit for more info",href="https://github.com/srtaheri/msstats-qc"))
      
                         )
                         
                                  ))
                       )
  ))