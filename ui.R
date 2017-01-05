library(shiny)
library(shinyBS)
library(plotly)
library(markdown)
shinyUI(fluidPage(
  shinyjs::useShinyjs(),
  titlePanel(title=p(strong("MSstatsQC"),align = "center",style="color:#444444;",style="font-size:150%;",
               style="font-family:inherit;"),windowTitle = "MSstatsQC"),
  navbarPage(h4("System suitability monitoring tools for quantitative mass spectrometry based proteomic
                experiments"),
              tabPanel("Home", theme = "bootstrap.css",
                         tags$img(src='logo.png', height=220, width=220, style = "float: right"),
                         p("MSstatsQC is an open-source web-based software which provides longitudinal
                           system suitability monitoring tools (control charts) for SRM based proteomic experiments."),
                        h5(strong("Metrics you can monitor")),
                        p("MSstatsQC uses control charts to monitor the instrument performance by tracking system
                           suitability metrics including total peak area, retention time and full width at half maximum (FWHM) and peak assymetry."),
                       h5(strong("Statistical functionalities")),
                        p("This framework includes simultaneous monitoring tools for mean and dispersion of suitability metrics and presents
                           alternative methods of monitoring such as time weighted control charts to ensure that various types
                           of process disturbances are detected effectively. Simultaneous control charts used in this framework
                           can be classified into two groups: individual-moving range (XmR) control charts and mean and dispersion
                           cumulative sum (CUSUM) control charts. To successfully identify the time of change, change point analysis
                           is also included in this framework. Experiment specific control limits are provided with the control 
                           charts to distinguish between random noise and systematic error."),
                       h5(strong("Using MSstatsQC")),
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
                         h5("Eralp Dogu,",span("eralp.dogu@gmail.com",style = "color:blue")),
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

              tabPanel("Data Import", 
                       sidebarLayout(

                        sidebarPanel(
                          wellPanel(
                            p("Please upload your data (Comma-separated (*.csv) 8 column QC file format)"),
                            p("To see an acceptable sample data, look at", strong("Help"),"tab"),
                            fileInput("filein", "Upload file"),
                            
                            p("Please select your preferred decision rule: "),
                            p(strong("Decision Rule: "),"To go process is when"),
                            fluidRow(
                              column(6,
                                     selectInput('peptideThreshold', '% of peptides', seq(0:100),
                                                 selected = 50)),
                              column(6,
                                     selectInput('metricThreshold', '# of metrics', seq(1:4),
                                                 selected = 2))
                            ),
                            p("are out of controll.")
                                        
                          ),
                           
                           wellPanel(
                             #p("If you want to run", strong("MSstatsQC"), "with sample data file, please click this button"),
                             actionButton("sample_button", "Run with sample data"),
                             bsTooltip("sample_button","If you want to run MSstatsQC with sample data file, please click this button", placement = "bottom", trigger = "hover",
                                       options = NULL)
                           ),
                           
                          wellPanel(
                            actionButton("clear_button", "Clear the data and plots"),
                            bsTooltip("clear_button","click this button to clear your data and all the tables and plots from the system.", placement = "bottom", trigger = "hover",
                                      options = NULL)
                          ),
                           
                          wellPanel(
                            p("Please select a guide set"),
                            numericInput("L","Lower bound of guide set",value = 1, min = 1, step = 1),
                            numericInput("U","Upper bound of guide set", value = 5, min = 2, step = 1),
                            p("Please select a precursor or select all"),
                            uiOutput("pepSelect")
                          ),
                           #br(),
                           wellPanel(
                             helpText("please select the columns of your data that you need to see."),
                             uiOutput("prodata_column_select")
                           ),
                           
                           tags$style("body{background-color:linen; color:black}")
                           
                           
                           ),
                        mainPanel(
                           
                          tabPanel("Data",  
                                   dataTableOutput('prodata_table'))
                        ), 
                           position = "left")
                       ),
              
              tabPanel("Metric Summary",
                       tabsetPanel(
                           tabPanel("Plot Summary",
                                    tags$head(tags$style(type="text/css")),
                                    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                     tags$div("It may take a while to load the plots, please wait...",
                                                              id="loadmessage")),
                                    fluidRow(
                                      column(9,
                                             plotOutput("plot_summary")
                                             ),
                                      column(3,
                                             wellPanel(
                                               textOutput("XmR_summary_decision_txt"),
                                               br(), br(),br(), br(),br(), br(),br(), br(),
                                               br(), br(),br(), br(),br(), br(),br(), br(),
                                               br(), br(),br(), br(),br(), br(),br(), br(),
                                               textOutput("CUSUM_summary_decision_txt"),
                                               br(), br(),br(), br(),br(), br(),br(), br(),
                                               br(), br(),br(), br(),br(), br(),br(), br(),
                                               br(), br(),br(), br(),br()
                                             )
                                             
                                      )
                                    )
                                     
                                    
                                    
                                    
           
                           ),
                         tabPanel("Boxplot",
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  plotlyOutput("box_plot", height = 2000)
           
                         ),
                         tabPanel("Scatterplot",
                                  uiOutput("scatter_plot_metric_selection"),
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  plotOutput("scatter_plot")
                                 ),
                         tabPanel("heat Map",
                                  tags$head(tags$style(type="text/css")),
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                   tags$div("It may take a while to load the plots, please wait...",
                                                            id="loadmessage")),
                                  plotOutput("heat_map")
                                  #,textOutput("heat_map_txt")
                                  
                         )
                         
                       )
                       ),
             
              navbarMenu("Control Charts",
                         tabPanel("XmR",
                                  
                                    sidebarLayout(
                                      sidebarPanel(
                                        uiOutput("XmR_select_metric")
                                      ), # end sidebarPanel
                                      mainPanel(
                                        uiOutput("XmR_tabset")
                                      ) # end mainPanel
                                    ) # end sidebarLayout
                                    
                                    
                                    

                                    #,
                                    # tabPanel("Mass Accuracy"
                                    #          , textOutput("MA_XmR_txt")
                                    #          ,plotlyOutput("MA_XmR")
                                    #          ,tags$head(tags$style(type="text/css"))
                                    #          , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    #                           tags$div("It may take a while to load the plots, please wait...
                                    #                                    ",id="loadmessage"))
                                    #                           )
                                             ), # End "XMR" tabPanel and it's tabsetPanel
                         
                         tabPanel("CUSUM",
                                  sidebarLayout(
                                    sidebarPanel(
                                      uiOutput("CUSUM_select_metric")
                                    ), # end sidebarPanel
                                    mainPanel(
                                      uiOutput("CUSUM_tabset")
                                    ) # end mainPanel
                                  ) # end sidebarLayout
                                    #,
                                    # tabPanel("Mass Accuracy" 
                                    #          , textOutput("MA_CUSUM_txt")
                                    #          ,plotlyOutput("MA_CUSUM")
                                    #          , tags$head(tags$style(type="text/css"))
                                    #          , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    #                           tags$div("It may take a while to load the plots, please wait...
                                    #                                    ",id="loadmessage"))
                                    #                           )
                                             )
                         #, # End tabpanle "CUSUM" and tabsetPanel of it
                         #tabPanel("EWMA", textOutput("EWMA_txt")),
                         #tabPanel("Short run SPC", textOutput("Short_run_SPC_txt")),
                       
                        # tabPanel("Multivariate Control Charts"
                        #          , textOutput("Multivariate_Control_Charts_txt"),
                        #          tabsetPanel(
                        #            tabPanel("Retention Time",
                        #                     plotlyOutput("RT_Multi")
                        #                     , tags$head(tags$style(type="text/css"))
                        #                     , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      ),
                        #            tabPanel("Total Peak Area", 
                        #                     plotlyOutput("TA_Multi")
                        #                     , tags$head(tags$style(type="text/css"))
                        #                     , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      ),
                        #            tabPanel("Full Width at Half Maximum (FWHM)", 
                        #                     plotlyOutput("Max_Multi")
                        #                     , tags$head(tags$style(type="text/css"))
                        #                     , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      ),
                        #            tabPanel("Peak Assymetry", 
                        #                     plotlyOutput("PA_Multi")
                        #                     , tags$head(tags$style(type="text/css"))
                        #                     , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                        #                                      tags$div("It may take a while to load the plots, please wait...
                        #                                               ",id="loadmessage"))
                        #                                      )
                                   #,
                                   # tabPanel("Mass Accuracy", 
                                   #          plotlyOutput("MA_Multi")
                                   #          , tags$head(tags$style(type="text/css"))
                                   #          , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                   #                           tags$div("It may take a while to load the plots, please wait...
                                   #                                    ",id="loadmessage"))
                                   #                           )
                                    #        ))
), # End "Multi" tabPanel and it's tabsetPanel
              
              tabPanel("Change Point Analysis",
                       sidebarLayout(
                         sidebarPanel(
                           uiOutput("CP_select_metric")
                         ), # end sidebarPanel
                         mainPanel(
                           uiOutput("CP_tabset")
                         ) # end mainPanel
                       ) # end sidebarLayout
                         #,
                         # tabPanel("Mass Accuracy"
                         #          , textOutput("CP_MA_txt")
                         #          , plotlyOutput("MA_CP")
                         #          , tags$head(tags$style(type="text/css"))
                         #          , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         #                           tags$div("It may take a while to load the plots, please wait...
                         #                                    ",id="loadmessage"))
                         #                           )
                         #          )), # end "Change Point Analysis" tabPanel and it's tabsetPanel
              
              # tabPanel("Capability Analysis",
              #          tabsetPanel(
              #            tabPanel("Retention Time"
              #                     , textOutput("CA_RT_txt")
              #                     , plotlyOutput("RT_CA")
              #                     , tags$head(tags$style(type="text/css"))
              #                     , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
              #                                      tags$div("It may take a while to load the plots, please wait...
              #                                               ",id="loadmessage"))
              #                                      ),
                         # tabPanel("Mass Accuracy"
                         #          , textOutput("CA_MA_txt")
                         #          , plotlyOutput("MA_CA")
                         #          , tags$head(tags$style(type="text/css"))
                         #          , conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         #                           tags$div("It may take a while to load the plots, please wait...
                         #                                    ",id="loadmessage"))
                         #                           )
                                  ), # end "Capability Analysis" tabPanel and it's tabsetPanel
              #tabPanel("Overall QC Performance", textOutput("OverallQC_txt")),
              tabPanel("Help",
                       tabsetPanel(
                         tabPanel("Metrics"
                                  ,h5(strong("Retention Time")),
                                  p("Retention time is the time it takes a solute to travel through the column. The retention time is assigned to 
                                    the corresponding solute peak. The retention time is a measure of the amount of time a solute spends in a column. 
                                    It is the sum of the time spent in the stationary phase and the mobile phase."
                                    ,a("visit for more info",href="http://www.britannica.com/science/retention-time")),
                                  br(),
                                  h5(strong("Total Peak Area")),
                                  p("Total Peak Area is the sum of all integrated signals for a certain peptide."),
                                  br(),
                                  h5(strong("Full Width at Half Maximum (FWHM)")),
                                  p("Full width at half maximum 'FWHM' is an expression of the extent of a 
                                    function given by the difference between the two extreme values 
                                    of the independent variable at which the dependent variable is equal 
                                    to half of its maximum value." 
                                    ,a("visit for more info",href="https://en.wikipedia.org/wiki/Full_width_at_half_maximum")),
                                  br(),
                                  h5(strong("Peak Assymetry")),
                                  p("Peak Assymetry is a measure of symetry for a peak. Calculated by taking 2*a/(a+b). Optimal value is around 1 for a Gaussian peak.")
                                  ),
                         
                         tabPanel("Plots"
                                  
                                  ,h5(strong("XmR control charts")),
                                  h5("Can detect large shifts and spikes in the mean and dispersion of suitability metric."),
                                  h5("The sequential differences between successive values as a measure 
                                    of dispersion and individual observations are used to construct the plots."),
                                  p("A measure of dispersion can be estimated by computing the ranges of two consecutive observations. This 
                                  approach is used to construct XmR charts. XmR charts plot 
                                    original observations and moving ranges to investigate deviations from random process behaviour.
                                    XmR chart consists of 2 charts; 
                                    Individuals (X) chart and Moving Range (mR) chart. 
                                    ", 
                                  a("visit for more info",href="https://en.wikipedia.org/wiki/Shewhart_individuals_control_chart")),
                                  
                                  br(),
                                  h5(strong("CUSUMm and CUSUMv control charts")),
                                  h5("Can detect small shifts and sustained drifts in the mean and dispersion of suitability metric."),
                                  h5("Time weighted cumulative sums are used to construct the plots."),
                                    p("A CUSUM chart is a time-weighted control chart that displays the cumulative sums 
                                    'CUSUMs'. Because it is cumulative, even minor drifts in the process mean or dispersion and gradual deterioration of quality will cause steadily 
                                    increasing or decreasing cumulative values. We introduce two CUSUM control charts: a mean CUSUM (CUSUMm) and 
                                      a dispersion CUSUM (CUSUMv).", 
                                  a("visit for more info",href="https://en.wikipedia.org/wiki/CUSUM")),
                               
                                  br(),
                                  h5(strong("Change Point Analysis")),
                                  h5("Can identify the exact time of a change in the mean and dispersion of suitability metric. "),
                                  h5("Likelihood functions are plotted and the QC sample which maximizes the functions is considered as a candidate change point."),
                                  p("A change in the process parameters triggers a control chart to generate an out of
                                    control signal. The QC sample at which the signal is issued is considered as the 
                                    stopping time and after the signal search for an assignable cause is recommended. 
                                    However, the signal does not always designate that the special cause actually occurred 
                                    at that certain time. A remedy to this problem is to use follow-up change 
                                    point analysis along with control charts. Change point estimation procedures have a potential 
                                    to save time by narrowing the search window. We introduce 
                                    two change point models: step shift change model for mean and step shift change model for variance.", 
                                   a("visit for more info",href="http://www.eng.fsu.edu/~pigna/pdf/"))
                         ),
                         
                         tabPanel("Documentation",
                                  h5(strong("MSstatsQC webpage")),
                                  p("Source codes, related documents and user manual can be found via our MSstats website"
                                    , a("visit for more info",href="http://www.msstats.org/msstatsqc")),
                                  br(),
                                  h5(strong("MSstatsQC Github")),
                                  p("Latest documantation is also available via our Github page"
                                    , a("visit for more info",href="https://github.com/srtaheri/msstats-qc"))
      
                         )
                         
                                  ))
                       )
  ))