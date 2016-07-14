library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)
library(RecordLinkage)
library(hash)
library(gridExtra)
library(markdown)

source("plot-functions.R")
source("data-validation.R")
source("helper-functions.R")


shinyServer(function(input,output,session) {
  
  #prodata <- data.frame()
  #### Read data  ##################################################################################################
  data <- reactiveValues(df = NULL)
  
  observeEvent(input$filein, {
    file1 <- input$filein
    data$df <- input_checking(read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE))
  })
  
  observeEvent(input$sample_button, {
    data$df <- read.csv("./Datasets/Sampledata_CPTAC_Study_9_1_Site54.csv")
  })
  
  observeEvent(input$clear_button, {
    data$df <- NULL
  })
  ##### Precursor type selection #####################################################################################
  output$pepSelect <- renderUI({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    selectInput("pepSelection","Choose precursor type"
                #,choices = c(levels(prodata$Precursor)
                            ,choices = c(levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
                            ,"all peptides"))
  })
  #### selecting columns to view in Data Import section ##############################################################
  output$prodata_column_select <- renderUI({
    prodata <- data$df
    checkboxGroupInput("show_prodata_columns", "columns of your data", choices = colnames(prodata), selected = colnames(prodata))
  })
  ######Show table of data #####################################################################################################
  output$prodata_table <- renderDataTable({
    data$df[,input$show_prodata_columns, drop = FALSE]
  }, options = list(pageLength = 25))
  ################################################################# plots ###################################################
  
  
  ################################################################################################################
  output$XmR_select_metric <- renderUI({

    checkboxGroupInput("XmR_checkbox_select","choose your prefered metric to view plots",
                       #choices = c("Retention Time" = "RT_XmR","Peak Assymetry" = "PA_XmR",
                         #          "Full Width at Half Maximum (FWHM)" = "Max_XmR","Total Peak Area" = "TA_XmR"),
                       choices = c("Peak Assymetry",find_metrics(data$df)),
                       selected = c("Peak Assymetry","BestRetentionTime",
                            "MaxFWHM",
                            "TotalArea")
                       )

  })
  #################################################################################################################
  map_XmR <- hash(keys= c("BestRetentionTime","Peak Assymetry",
                          "MaxFWHM",
                          "TotalArea","metric1","metric2","meric3","metric4","metric5"),
                     values=c("RT_XmR","PA_XmR","Max_XmR","TA_XmR","m1","m2","m3","m4","m5"))
  
  output$XmR_tabset <- renderUI({
    
    Tabs <- lapply(input$XmR_checkbox_select,
                   function(x) {
                     tabPanel(x,
                              plotlyOutput(map_XmR[[x]]),
                              tags$head(tags$style(type="text/css")),
                              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                               tags$div("It may take a while to load the plots, please wait...",
                                                        id="loadmessage"))
                              )
                     })
    
    do.call(tabsetPanel, Tabs)
  })
  ####################################################### plot XmR for RT ###################################################################################################################
  output$RT_XmR <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Retention Time", plot.method = "XmR", normalization.type = FALSE, y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR for Peak assymetry #####################################################################################################
  output$PA_XmR <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Peak Assymetry", plot.method = "XmR", normalization.type = FALSE, y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR FOR MaxFWHM ############################################################################################################
  output$Max_XmR <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "FWHM", plot.method = "XmR", normalization.type = FALSE, y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR FOR total area #########################################################################################################
  output$TA_XmR <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Total Area", plot.method = "XmR", normalization.type = FALSE, y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot CUSUM_chart  for RT##########################################################################################################
  output$RT_CUSUM <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Retention Time", plot.method = "CUSUM", normalization.type = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################plot CUSUM for Peak assymetry ######################################################################################################
  output$PA_CUSUM <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Peak Assymetry", plot.method = "CUSUM", normalization.type = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################### plot CUSUM FOR MaxFWHM #########################################################################################################
  output$Max_CUSUM <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "FWHM", plot.method = "CUSUM", normalization.type = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")    
  })
  ############################################################ plot CUSUM FOR total area ####################################
  output$TA_CUSUM <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Total Area", plot.method = "CUSUM", normalization.type = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################## plot Change Point for RT #############################
  output$RT_CP <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Retention Time", plot.method = "CP", normalization.type = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ######################################################## plot Change Point for Peak assymetry ###########
  output$PA_CP <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Peak Assymetry", plot.method = "CP", normalization.type = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ########################################################## plot Change Point FOR MaxFWHM ####################################
  output$Max_CP <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "FWHM", plot.method = "CP", normalization.type = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation")    
  })
  ########################################################## plot Change Point FOR total area ##################################
  output$TA_CP <- renderPlotly({
    render.QC.chart(data$df, input$pepSelection, input$L, input$U, normalize.metric = "Total Area", plot.method = "CP", normalization.type = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })

  ########################################################## box plot in Summary tab ##########################################
  output$box_plot <- renderPlotly({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    metrics_box.plot(prodata)
  })
  ########################################################## scatterplot matrix in Summary tab #################################
  output$scatter_plot <- renderPlot({
   
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    metrics_scatter.plot(prodata, input$L, input$U, input$metric_precursor, normalization = TRUE)
  }, height = 700)
  ######################################################### plot_summary in Summary tab ########################################
  # output$plot_summary <- renderPlot({
  # 
  #   prodata <- data$df
  #   validate(
  #     need(!is.null(prodata), "Please upload your data")
  #   )
    #p1 <- CUSUM.Summary.plot(prodata, input$L, input$U)
    #p2 <- CUSUM.Radar.Plot(prodata,input$L,input$U)
    #p3 <- XmR.Summary.plot(prodata, input$L, input$U)
    #p4 <- XmR.Radar.Plot(prodata,input$L,input$U)
    #grid.arrange(p1,p2,p3,p4, ncol = 1)
    #XmR.Radar.Plot(prodata,input$L,input$U, metric = "Retention Time")
    
  #}, height = 500)
  output$plot_summary <- renderPlotly({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
  XmR.Radar.Plot(prodata,input$L,input$U,metric = "Retention Time")
  XmR.Radar.Plot(prodata,input$L,input$U,metric = "Peak Assymetry")
  XmR.Radar.Plot(prodata,input$L,input$U,metric = "FWHM")
  XmR.Radar.Plot(prodata,input$L,input$U,metric = "Total Area")
  })
  
  ###########################################################################################################################
  ###########################################################################################################################
  ########################################################## "help" tab ################################
  
 
  ############################################################################################################################
})