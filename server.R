library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)
library(RecordLinkage)
library(hash)
library(gridExtra)
library(ggExtra)
library(markdown)

source("plot-functions.R")
source("data-validation.R")
source("helper-functions.R")
source("QCMetrics.R")

shinyServer(function(input,output,session) {
  COL.BEST.RET <- "Retention Time"
  COL.FWHM <- "Full Width at Half Maximum"
  COL.TOTAL.AREA <- "Total Peak Area"
  COL.PEAK.ASS <- "Peak Assymetry"
  

  #### Read data  ##################################################################################################
  data <- reactiveValues(df = NULL, metrics = NULL)
  
  
  observeEvent(input$filein, {
    file1 <- input$filein
    data$df <- input_checking(read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE))
    data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
  }, priority = 20)
  
  observeEvent(input$sample_button, {
    data$df <- input_checking(read.csv("./Datasets/Sampledata_CPTAC_Study_9_1_Site54.csv"))
    data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
  }, priority = 20)
  
  observeEvent(input$clear_button, {
    data$df <- NULL
    data$metrics <- NULL
  }, priority = 20)
  ##### Precursor type selection #####################################################################################
  output$pepSelect <- renderUI({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    selectInput("pepSelection","Choose precursor type"
                            ,choices = c(levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
                            ,"all peptides"))
  })
  #### selecting columns to view in Data Import section ##############################################################
  output$prodata_column_select <- renderUI({
    prodata <- data$df
    checkboxGroupInput("show_prodata_columns", "columns of your data", choices = colnames(prodata), selected = colnames(prodata))
  })
  ######Show table of data #####################################################################################################
  output$prodata_table <- renderDataTable({
    data$df[,input$show_prodata_columns, drop = FALSE] # drop = F, is for not considering the last column as arrow and consider it as data frame
  }, options = list(pageLength = 25))
  ################################################################# plots ###################################################
  ###########################################################################################################################
  output$XmR_select_metric <- renderUI({

    checkboxGroupInput("XmR_checkbox_select","choose your prefered metric to view plots",
                       choices = c(data$metrics),
                       selected = c(COL.PEAK.ASS,COL.BEST.RET,
                                    COL.FWHM, COL.TOTAL.AREA)
                       )

  })
  #################################################################################################################
  
  output$XmR_tabset <- renderUI({
    
    Tabs <- lapply(input$XmR_checkbox_select,
                   function(x) {
                       tabPanel(x,
                                tags$head(tags$style(type="text/css")),
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 tags$div("It may take a while to load the plots, please wait...",
                                                          id="loadmessage")),
                                renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L,
                                                             input$U, metric = x,
                                                             plot.method = "XmR", normalization = FALSE,
                                                             y.title1 = "Individual Value", y.title2 = "Moving Range"))
                                
                                )
                   })
    do.call(tabsetPanel, Tabs)
  })
  ################################################################################################################
  ################################################################################################################
  output$CUSUM_select_metric <- renderUI({
    
    checkboxGroupInput("CUSUM_checkbox_select","choose your prefered metric to view plots",
                       choices = c(data$metrics),
                       selected = c(COL.PEAK.ASS,COL.BEST.RET,
                                    COL.FWHM,
                                    COL.TOTAL.AREA)
    )
    
  })
  #################################################################################################################
  output$CUSUM_tabset <- renderUI({
    
    Tabs <- lapply(input$CUSUM_checkbox_select,
                   function(x) {
               
                       tabPanel(x,
                                tags$head(tags$style(type="text/css")),
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 tags$div("It may take a while to load the plots, please wait...",
                                                          id="loadmessage")),
                                renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L, input$U, metric = x, plot.method = "CUSUM", normalization = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation"))
                                )
                   })
    
    do.call(tabsetPanel, Tabs)
  })
  ################################################################################################################
  ################################################################################################################
  output$CP_select_metric <- renderUI({
    
    checkboxGroupInput("CP_checkbox_select","choose your prefered metric to view plots",
                       choices = c(data$metrics),
                       selected = c(COL.PEAK.ASS,COL.BEST.RET,
                                    COL.TOTAL.AREA,
                                    COL.FWHM)
    )
    
  })
  #################################################################################################################
  output$CP_tabset <- renderUI({
    
    Tabs <- lapply(input$CP_checkbox_select,
                   function(x) {
                       tabPanel(x,
                                tags$head(tags$style(type="text/css")),
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 tags$div("It may take a while to load the plots, please wait...",
                                                          id="loadmessage")),
                                renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L, input$U, metric = x, plot.method = "CP", normalization = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation"))
                                )
                   })
    
    do.call(tabsetPanel, Tabs)
  })
  ########################################################## box plot in Summary tab ##########################################
  output$box_plot <- renderPlotly({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    metrics_box.plot(prodata, data.metrics = data$metrics)
  })
  ########################################################## scatterplot matrix in Summary tab #################################
  output$scatter_plot_metric_selection <- renderUI({

     return(selectInput("scatter_metric_select", "Choose the metric",
                 choices = c(data$metrics), selected = COL.FWHM))
    
  })
  output$scatter_plot <- renderPlot({
    prodata <- data$df
      validate(
        need(!is.null(prodata), "Please upload your data")
      )
    metrics_scatter.plot(prodata,input$L, input$U, metric = input$scatter_metric_select, normalization = T)
  }, height = 1500) 
  outputOptions(output, "scatter_plot_metric_selection", priority = 10)
  outputOptions(output, "scatter_plot", priority = 1)
  ######################################################### plot_summary in Summary tab ########################################
  my_height <- reactive({
    my_height <- ceiling(length(data$metrics)/4)*1200
  })
  
  output$plot_summary <- renderPlot({

    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    #Text1 = textGrob("Mean",gp=gpar(col="gray40",fontsize = 11,face="bold"))
    #Text2 = textGrob("Dispersion",gp=gpar(col="gray40",fontsize = 11,face="bold"))
    p1 <- XmR.Summary.plot(prodata, data.metrics = data$metrics, input$L, input$U)
    p2 <- XmR.Radar.Plot(prodata, data.metrics = data$metrics,input$L,input$U)
    p3 <- CUSUM.Summary.plot(prodata, data.metrics = data$metrics, input$L, input$U)
    p4 <- CUSUM.Radar.Plot(prodata, data.metrics = data$metrics, input$L,input$U)
    #p1 <- p1 + annotation_custom(grob = Text1, xmin = 0, xmax = 7, ymin = 0, ymax = 2.6)
    #p1 <- p1 + annotation_custom(grob = Text2, xmin = 0, xmax = 15, ymin = 0, ymax = -2.6)
    #p3 <- p3 + annotation_custom(grob = Text1, xmin = 0, xmax = 7, ymin = 0, ymax = 2.6)
    #p3 <- p3 + annotation_custom(grob = Text2, xmin = 0, xmax = 15, ymin = 0, ymax = -2.6)
    grid.arrange(p1,p2,p3,p4, ncol = 1)

  }, height = my_height )

  #1500
  ############################# heat_map in Summary tab #############################################

  # output$heat_map_metric_selection <- renderUI({
  #   
  #   checkboxGroupInput("heat_map_checkbox_select","choose your prefered metric to view plots",
  #                      choices = c(data$metrics),
  #                      selected = c(COL.PEAK.ASS,COL.BEST.RET,
  #                                   COL.FWHM, COL.TOTAL.AREA)
  #   )
  #   
  # })
###################################################  
  output$heat_map <- renderPlot({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    validate(
      need(!is.null(prodata$AcquiredTime),"To view heatmap, the data set should include AcquiredTime column.")
    )
    if(is.null(prodata$AcquiredTime)) return(NULL)
    metricData <- getMetricData(prodata, precursor = input$pepSelection, input$L, input$U, metric = COL.PEAK.ASS, normalization = FALSE)
    metrics_heat.map(prodata,metricData ,precursorSelection = input$pepSelection, input$L, input$U, type = 1)
  })
 
  ############################################################################################################################
})