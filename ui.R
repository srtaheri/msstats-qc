library(shiny)
library(plotly)

shinyUI(fluidPage(
  titlePanel("MSstatsQC: Quality control tools for LC MS/MS proteomic experiments"),
  navbarPage( "",
              tabPanel("Home", theme = "bootstrap.css",
                         tags$img(src='logo.png', height=200, width=200, style = "float: right"),
                         
                         #h1("MSstatsQC"),
                         #br(),
                         p("MSstatsQC is an open-source web-based software which provides longitudinal
                           system suitability monitoring tools (control charts) for LC MS/MS proteomic experiments."),
                         p("MSstatsQC uses control charts to monitor the chromatographic performance by tracking system
                           suitability metrics including peak area, retention time and full width at half maximum (FWHM)."),
                         p("MSstatsQC includes simultaneous monitoring tools for metric-wise mean and dispersion and presents
                           alternative methods of monitoring such as time weighted control charts to ensure that various types
                           of process disturbances are detected effectively. Simultaneous control charts used in this framework
                           can be classified into two groups: individual-moving range (ZmR) control charts and mean and dispersion
                           cumulative sum (CUSUM) control charts. To successfully identify the time of change, change point analysis
                           is also included in this framework. Experiment specific control limits are provided with the control 
                           charts to distinguish between random noise and systematic error."),
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
                         tags$img(src='XmRchart.png', height=800, width=1000, style = "float: bottom"),
                         p("Project Team: "),
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
              ),
              tabPanel("Data Import", theme = "bootstrap.css",
<<<<<<< HEAD
                       sidebarLayout(
                         
                         sidebarPanel(
=======
                       sidebarLayout( 
                         mainPanel(
                           tabsetPanel('Data',  DT::dataTableOutput('prodata_table'))
                           ), 
                      
                        sidebarPanel(
>>>>>>> origin/master
                           
                           p("Please upload your data (Comma-separated (*.csv) 12 column QC file format)"),
                           p("To see an acceptable sample data, look at", strong("Help"),"tab"),
                           fileInput("filein", "Upload file"),
                           p("Please select a guide set"),
                           numericInput("L","Lower bound of guide set",value = 1, min = 1, step = 1),
                           numericInput("U","Upper bound of guide set", value = 5, min = 2, step = 1),
                           p("Please select a precursor or select all"),
                           uiOutput("pepSelect"),
                           p("If you have uploaded your data set, selected the guidset and chosen your precursor type
                             , click on this button to see the plots"),
                           actionButton("act_button", "click to see plots"),
                           tags$style("body{background-color:linen; color:black}")
                           ),
                         mainPanel(tableOutput("prodata_table")),
                           position = "left")
                       ),
              
              tabPanel("Metric Summary",
                       tabsetPanel(
                         tabPanel("Boxplot",
                                  plotlyOutput("box_plot", height = 2000)
                         ),
                         tabPanel("Scatterplot", #helpText(")
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
                                                              ),
                                    tabPanel("Mass Accuracy"
                                             , textOutput("MA_ZMR_txt")
                                             ,plotlyOutput("MA_ZMR")
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
                                                              ),
                                    tabPanel("Mass Accuracy" 
                                             , textOutput("MA_CUSUM_txt")
                                             ,plotlyOutput("MA_CUSUM")
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
                                             )), # End tabpanle "CUSUM" and tabsetPanel of it
                         tabPanel("EWMA", textOutput("EWMA_txt")),
                         tabPanel("Short run SPC", textOutput("Short_run_SPC_txt")),
                        # end navbarMenu
                        tabPanel("Multivariate Control Charts"
                                 , textOutput("Multivariate_Control_Charts_txt"),
                                 tabsetPanel(
                                   tabPanel("Retention Time",
                                            plotlyOutput("RT_Multi")
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
                                            plotlyOutput("TA_Multi")
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
                                            plotlyOutput("Max_Multi")
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
                                            plotlyOutput("PA_Multi")
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
                                   tabPanel("Mass Accuracy", 
                                            plotlyOutput("MA_Multi")
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
                                            ))), # End "Multi" tabPanel and it's tabsetPanel
              
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
                                                   ),
                         tabPanel("Mass Accuracy"
                                  , textOutput("CP_MA_txt")
                                  , plotlyOutput("MA_CP")
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
                                  )), # end "Change Point Analysis" tabPanel and it's tabsetPanel
              
              tabPanel("Capability Analysis",
                       tabsetPanel(
                         tabPanel("Retention Time"
                                  , textOutput("CA_RT_txt")
                                  , plotlyOutput("RT_CA")
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
                         tabPanel("Mass Accuracy"
                                  , textOutput("CA_MA_txt")
                                  , plotlyOutput("MA_CA")
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
                                  )), # end "Capability Analysis" tabPanel and it's tabsetPanel
              tabPanel("Overall QC Performance", textOutput("OverallQC_txt")),
              tabPanel("Help",
                       tabsetPanel(
                         tabPanel("Metrics"
                                  #, plotOutput("help_metric", hover = "help_metric_hover")
                                  #, verbatimTextOutput("help_metric_info")
                                  ,h5(strong("Retention Time")),
                                  p("Retention time is the time it takes a solute to travel through the column. The retention time is assigned to 
                                    the corresponding solute peak. The retention time is a measure of the amount of time a solute spends in a column. 
                                    It is the sum of the time spent in the stationary phase and the mobile phase."
                                    ,a("visit for more info",href="http://www.britannica.com/science/retention-time")),
                                  
                                  h5(strong("Total Peak Area")),
                                  
                                  h5(strong("Full Width at Half Maximum (FWHM)")),
                                  p("is an expression of the extent of a function given by the difference between the
                                    two extreme values of the independent variable at which the dependent variable is 
                                    equal to half of its maximum value. In other words, it is the width of a spectrum 
                                    curve measured between those points on the y-axis which are half the maximum amplitude."),
                                  
                                  h5(strong("Peak Assymetry")),
                                  p("Calculated by taking 2*a/(a+b). Optimal value is around 1 for a Gaussian peak")
                                  ),
                         
                         tabPanel("Plots"
                                  #, plotlyOutput("Exxx")
                                  ,h5(strong("CUSUMm and CUSUMv control charts")),
                                  
                                  h5(strong("Z and MR control charts")),
                                  p("By using the sequential differences between successive values as a measure 
                                    of dispersion, a chart for standardized individual observations (_𝑧_𝑖_) and 
                                    a moving range chart can be created "),
                                  
                                  h5(strong("Change Point Analysis")),
                                  p("A change in the process parameters triggers a control chart to generate an out of
                                    control signal. The QC sample at which the signal is issued is considered as the 
                                    stopping time and after the signal search for an assignable cause is recommended.")
                                  ),
                         tabPanel("Examples", uiOutput("video"),
                                  h5(strong("AAAA")),
                                  
                                  h5(strong("BBBBB")),
                                  p("By using the sequential differences between successive values as a measure 
                                    of dispersion, a chart for standardized individual observations (_𝑧_𝑖_) and 
                                    a moving range chart can be created "),
                                  
                                  h5(strong("Change Point Analysis")),
                                  p("A change in the process parameters triggers a control chart to generate an out of
                                    control signal. The QC sample at which the signal is issued is considered as the 
                                    stopping time and after the signal search for an assignable cause is recommended.")
                         ),
                         tabPanel("Troubleshooting",
                                  h5(strong("Missing data or poor quality peaks")),
                                  p(" "),
                                  h5(strong("Retention time drift")),
                                  p(" "),
                                  h5(strong("Change Point Analysis"))
                         )
                         
                                  ))
                       )
  ))