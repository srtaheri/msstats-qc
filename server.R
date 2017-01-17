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
    validate(
      need(!is.null(data$df), "Please upload your data"),
      need(is.data.frame(data$df), data$df)
    )
    data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
  }, priority = 20)

  observeEvent(input$sample_button, {
    data$df <- input_checking(read.csv("./Datasets/Sampledata_CPTAC_Study_9_1_Site54.csv"))
    validate(
      need(!is.null(data$df), "Please upload your data"),
      need(is.data.frame(data$df), data$df)
    )
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
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata)
    )
    selectInput("pepSelection","Choose precursor type"
                            ,choices = c(levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
                            ,"all peptides"))
  })
  #### selecting columns to view in Data Import section ##############################################################
  # output$prodata_column_select <- renderUI({
  #   prodata <- data$df
  #   checkboxGroupInput("show_prodata_columns", "columns of your data", choices = colnames(prodata), selected = colnames(prodata))
  # })
  ######Show table of data #####################################################################################################
   output$prodata_table <- renderDataTable({
     validate(
       need(!is.null(data$df), "Please upload your data"),
       need(is.data.frame(data$df), data$df)
     )
     #data$df[,input$show_prodata_columns, drop = FALSE] # drop = F, is for not considering the last column as arrow and consider it as data frame
     data$df
   }, options = list(pageLength = 25))

  ###### Tab for selecting decision rule and metrics ###############################################
  output$metricThresholdGood <- renderUI({
    numOfMetrics <- length(input$user_selected_metrics)
    numericInput('threshold_metric_good', '', value = 1, min = 1, max = numOfMetrics, step = 1)
  })

  output$peptideThresholdWarn <- renderUI({
    threshold_peptide_good <- input$threshold_peptide_good
    numericInput('threshold_peptide_warn', '', value = threshold_peptide_good+1, min = threshold_peptide_good+1, max = 100, step = 1)
  })

  output$metricThresholdWarn <- renderUI({
    validate(
      need(!is.null(input$threshold_metric_good),"loading...")
    )
    numOfMetrics <- length(input$user_selected_metrics)
    threshold_metric_good <- input$threshold_metric_good
    numericInput('threshold_metric_warn', '', value = threshold_metric_good, min = threshold_metric_good, max = numOfMetrics, step = 1)
  })

  output$metricSelection <- renderUI({
    checkboxGroupInput("user_selected_metrics","Please select the metrics",
                       choices = c(data$metrics),
                       selected = c(COL.PEAK.ASS,COL.BEST.RET,
                                    COL.FWHM, COL.TOTAL.AREA),
                       inline = TRUE)
  })


  ###########################################################################

  # output$ToGoPeptidesSelect <- renderUI({
  #   numericInput('peptideThresholdGood', '% of peptides', value = 50, min = 0,
  #                max = 100, step = 1)
  # })
  # output$ToGoPeptides <- renderText({
  #   input$peptideThresholdGood
  # })
  # output$ToGoMetricSelect <- renderUI({
  #   numericInput('metricThresholdGood', '# ', value = 1, min = 1,
  #                max = length(input$user_selected_metrics), step = 1)
  # })
  # output$ToGoMetrics <- renderText({
  #   input$metricThresholdGood
  # })
  ################################################################# plots ###################################################
  ###########################################################################################################################
  # output$XmR_select_metric <- renderUI({
  #
  #   checkboxGroupInput("XmR_checkbox_select","choose your prefered metric to view plots",
  #                      choices = c(data$metrics),
  #                      selected = c(COL.PEAK.ASS,COL.BEST.RET,
  #                                   COL.FWHM, COL.TOTAL.AREA)
  #                      )
  #
  # })
  #################################################################################################################

  output$XmR_tabset <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select your decision rule and metrics from the selection tab")
    )
    Tabs <- lapply(input$user_selected_metrics,
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
  # output$CUSUM_select_metric <- renderUI({
  #
  #   checkboxGroupInput("CUSUM_checkbox_select","choose your prefered metric to view plots",
  #                      choices = c(data$metrics),
  #                      selected = c(COL.PEAK.ASS,COL.BEST.RET,
  #                                   COL.FWHM,
  #                                   COL.TOTAL.AREA)
  #   )
  #
  # })
  #################################################################################################################
  output$CUSUM_tabset <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select your decision rule and metrics from the selection tab")
    )
    Tabs <- lapply(input$user_selected_metrics,
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
  # output$CP_select_metric <- renderUI({
  #
  #   checkboxGroupInput("CP_checkbox_select","choose your prefered metric to view plots",
  #                      choices = c(data$metrics),
  #                      selected = c(COL.PEAK.ASS,COL.BEST.RET,
  #                                   COL.TOTAL.AREA,
  #                                   COL.FWHM)
  #   )
  #
  # })
  #################################################################################################################
  output$CP_tabset <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select your decision rule and metrics from the selection tab")
    )
    Tabs <- lapply(input$user_selected_metrics,
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
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select your decision rule and metrics from the selection tab")
    )
    metrics_box.plot(prodata, data.metrics = data$metrics)
  })

  ######################################################### plot_summary in Summary tab ########################################
  my_height <- reactive({
    my_height <- ceiling(length(input$user_selected_metrics)*length(input$summary_controlChart_select))*200
  })
  heatmap_width <- reactive({
    prodata <- data$df
    heatmap_width <- nrow(prodata[prodata$Precursor == prodata$Precursor[1],])*20
    heatmap_width
  })
  ################ plot the summary and radar plots ############################################################################
  output$plot_summary <- renderPlot({

    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select your decision rule and metrics from the selection tab")
    )

    p1 <- XmR.Summary.plot(prodata, data.metrics = input$user_selected_metrics, input$L, input$U)
    p2 <- XmR.Radar.Plot(prodata, data.metrics = input$user_selected_metrics,input$L,input$U)
    p3 <- CUSUM.Summary.plot(prodata, data.metrics = input$user_selected_metrics, input$L, input$U)
    p4 <- CUSUM.Radar.Plot(prodata, data.metrics = input$user_selected_metrics, input$L,input$U)

    if("XmR" %in% input$summary_controlChart_select && "CUSUM" %in% input$summary_controlChart_select) {
      grid.arrange(p1,p2,p3,p4, ncol = 1)
      print(length(input$summary_controlChart_select))
    }
    else if("XmR" %in% input$summary_controlChart_select) {
      grid.arrange(p1,p2, ncol = 1)
    }else if("CUSUM" %in% input$summary_controlChart_select) {
      grid.arrange(p3,p4,ncol = 1)
    }else {
    }

  }, height = my_height )

  ################## decision message for XmR in summary tab #########
   output$summary_decision_txt <- renderText({
     prodata <- data$df
     validate(
       need(!is.null(prodata), "Please upload your data"),
       need(is.data.frame(prodata), prodata),
       need(!is.null(input$user_selected_metrics),"Please first select the metrics and decision thresholds in the selection tab")
     )
     peptideThresholdGood <- (as.numeric(input$threshold_peptide_good))/100 #For Eralp : the default is set to 50 (look at the selection tab)
     metricThresholdGood <- as.numeric(input$threshold_metric_good) # For Eralp : the default is 1 (look at the selection tab)
     peptideThresholdWarn <- (as.numeric(input$threshold_peptide_warn))/100 #For Eralp : the default is set to 70 (look at the selection tab)
     metricThresholdWarn <- as.numeric(input$threshold_metric_warn) # For Eralp : the default is 2 (look at the selection tab)

     # For Eralp : This gives the number of out of range metrics for XmR "mean"(type = 1) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdGood.[1]
     XmRCounterAboveGood1 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdGood,peptideThresholdWarn,
                                                          input$L, input$U, type = 1)[1]
     # For Eralp : This gives the number of out of range metrics for XmR "mean"(type = 1) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdWarn.[2]
     XmRCounterAboveWarn1 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdGood,peptideThresholdWarn,
                                                                input$L, input$U, type = 1)[2]
     # For Eralp : This gives the number of out of range metrics for XmR "dispersion"(type = 2) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdGood.[1]
     XmRCounterAboveGood2 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdGood,peptideThresholdWarn,
                                                                input$L, input$U, type = 2)[1]
     # For Eralp : This gives the number of out of range metrics for XmR "dispersion"(type = 2) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdWarn.[2]
     XmRCounterAboveWarn2 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdGood,peptideThresholdWarn,
                                                                input$L, input$U, type = 2)[2]

  #   CUSUMCounter1 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "CUSUM", peptideThresholdGood,peptideThresholdWarn,
   #  input$L, input$U, type = 1)[1]
  #   CUSUMCounter2 <- CUSUM.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThreshold,
  #                                                         input$L, input$U, type = 2)

      if(XmRCounterAboveGood1 <= metricThresholdGood && XmRCounterAboveGood2 <= metricThresholdGood) {"System is in-control"}
     # if(XmRCounterAboveGood1 > metricThresholdGood && XmRCounterAboveGood1 <= metricThresholdWarn &&
     #    XmRCounterAboveGood2 <= metricThresholdGood){ "Warning! System is out-of-control (A change in QC metric mean is possible)"}
    # if(XmRCounterAboveGood2 > metricThresholdGood && XmRCounterAboveGood2 <= metricThresholdWarn &&
    #    XmRCounterAboveGood1 <= metricThresholdGood){"Warning! System is out-of-control (A change in QC metric variation is possible)"}
    # if(XmRCounterAboveGood1 > metricThresholdGood && XmRCounterAboveGood1 <= metricThresholdWarn &&
    #    XmRCounterAboveGood2 > metricThresholdGood && XmRCounterAboveGood2 <= metricThresholdWarn) {"Warning! System is out-of-control (A simultaneous change in QC metric mean and variation is possible)"}
    # if(XmRCounterAboveGood1 > metricThresholdWarn &&
    #    XmRCounterAboveGood2 > metricThresholdGood && XmRCounterAboveGood2 <= metricThresholdWarn) {"Bad! Mean is bad and variation is in warning area"}
    # if(XmRCounterAboveGood2 > metricThresholdWarn &&
    #    XmRCounterAboveGood1 > metricThresholdGood && XmRCounterAboveGood1 <= metricThresholdWarn) {"Bad! variation is bad and mean is in warning area"}
    # if(XmRCounterAboveGood1 > metricThresholdWarn && XmRCounterAboveGood2 > metricThresholdWarn) {"Bad! both mean and variation are in bad area"}
  })
  ############################# heat_map in Summary tab #############################################

  output$heat_map <- renderPlot({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select the metrics and decision thresholds in the selection tab"),
      need(!is.null(prodata$AcquiredTime),"To view heatmap, the data set should include AcquiredTime column.")
    )

    peptideThresholdGood <- (as.numeric(input$threshold_peptide_good))/100
    peptideThresholdWarn <- (as.numeric(input$threshold_peptide_warn))/100
    if(is.null(prodata$AcquiredTime)) return(NULL)

    p1 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "XmR",
                           peptideThresholdGood, peptideThresholdWarn,input$L, input$U, type = 1,
                           title = "XmR heat map - Mean")
    p2 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "XmR",
                           peptideThresholdGood, peptideThresholdWarn,input$L, input$U, type = 2,
                           title = "XmR heat map - Dispersion")
    p3 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "CUSUM",
                           peptideThresholdGood, peptideThresholdWarn,input$L, input$U, type = 1,
                           title = "CUSUM heat map - Mean")
    p4 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "CUSUM",
                           peptideThresholdGood, peptideThresholdWarn,input$L, input$U, type = 2,
                           title = "CUSUM heat map - Dispersion")
    if("XmR" %in% input$heatmap_controlChart_select && "CUSUM" %in% input$heatmap_controlChart_select) {
      grid.arrange(p1,p2,p3,p4, ncol = 1)
    }
    else if("XmR" %in% input$heatmap_controlChart_select) {
      grid.arrange(p1,p2, ncol = 1)
    }else if("CUSUM" %in% input$heatmap_controlChart_select) {
      grid.arrange(p3,p4,ncol = 1)
    }else {
    }
  }, height = my_height, width = heatmap_width)

  ############################################################################################################################
})
