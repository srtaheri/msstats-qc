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
    selectInput("pepSelection","Choose peptide"
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
  output$metricThresholdRed <- renderUI({
    numOfMetrics <- length(input$user_selected_metrics)
    numericInput('threshold_metric_red', '', value = 2, min = 0, max = numOfMetrics, step = 1)
  })

  output$peptideThresholdYellow <- renderUI({
    threshold_peptide_red <- input$threshold_peptide_red
    numericInput('threshold_peptide_yellow', '', value = threshold_peptide_red - 1, min = 1, max = threshold_peptide_red, step = 1)
  })

  output$metricThresholdYellow <- renderUI({
    # validate(
    #   need(!is.null(input$threshold_metric_red),"loading...")
    # )
    numOfMetrics <- length(input$user_selected_metrics)
    threshold_metric_red <- input$threshold_metric_red
    numericInput('threshold_metric_yellow', '', value = threshold_metric_red , min = 0, max = threshold_metric_red, step = 1)
  })

  output$metricSelection <- renderUI({
    checkboxGroupInput("user_selected_metrics","",
                       choices = c(data$metrics),
                       selected = c(COL.PEAK.ASS,COL.BEST.RET,
                                    COL.FWHM, COL.TOTAL.AREA),
                       inline = TRUE)
  })

  ################################################################# plots ###################################################
  #################################################################################################################

  output$XmR_tabset <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
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
  #################################################################################################################
  output$CUSUM_tabset <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
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
  #################################################################################################################
  output$CP_tabset <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
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
  ######################################################### height and width in Summary tab ########################################
  my_height <- reactive({
    my_height <- ceiling(length(input$user_selected_metrics)*length(input$summary_controlChart_select))*200
  })

  heatmap_width <- reactive({
    prodata <- data$df
    heatmap_width <- nrow(prodata[prodata$Precursor == prodata$Precursor[1],])*20
    heatmap_width
  })
  ########################################################## box plot in Metric Summary tab ##########################################
  output$box_plot <- renderPlotly({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
    )
    metrics_box.plot(prodata, data.metrics = data$metrics)
  })

  ###############   summary plots and radar plots ############################################################################
  output$plot_summary <- renderPlot({

    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
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
       need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule")
     )
     peptideThresholdRed <- (as.numeric(input$threshold_peptide_red))/100 #For Eralp : this is the percentage of peptide user chooses for red flag
     metricThresholdRed <- as.numeric(input$threshold_metric_red) # For Eralp : this is the number of metric user chooses for red flag
     peptideThresholdYellow <- (as.numeric(input$threshold_peptide_yellow))/100 #For Eralp : this is the percentage of peptide user chooses for yellow flag
     metricThresholdYellow <- as.numeric(input$threshold_metric_yellow) # For Eralp : this is the number of metric user chooses for yellow flag

     # For Eralp : This gives the number of out of range metrics for XmR "mean"(type = 1) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdRed.[1]
     XmRCounterAboveRed1 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,
                                                          input$L, input$U, type = 1)[1]
     # For Eralp : This gives the number of out of range metrics for XmR "mean"(type = 1) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdYellow and below the peptideThresholdRed..[2]
     XmRCounterAboveYellow1 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,
                                                                input$L, input$U, type = 1)[2]
     # For Eralp : This gives the number of out of range metrics for XmR "mean"(type = 2) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdRed.[1]
     XmRCounterAboveRed2 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,
                                                                input$L, input$U, type = 2)[1]
     # For Eralp : This gives the number of out of range metrics for XmR "mean"(type = 2) for which there exists a peptide
     #that its percentage of out of range is above the peptideThresholdYellow and below the peptideThresholdRed..[2]
     XmRCounterAboveYellow2 <- number.Of.Out.Of.Range.Metrics(prodata,data$metrics,method = "XmR", peptideThresholdRed,peptideThresholdYellow,
                                                                input$L, input$U, type = 2)[2]


      if(XmRCounterAboveYellow1 > metricThresholdYellow && XmRCounterAboveYellow2 > metricThresholdYellow) {"Yellow FLAG: System performance is POOR"}
      if(XmRCounterAboveRed1 > metricThresholdRed && XmRCounterAboveRed2 > metricThresholdRed) {paste("RED FLAG: System performance is UNACCEPTABLE","<fontcolor=\"#FF0000\"><b>", input$n, "</b></font>")}
  })
  ############################# heat_map in Summary tab #############################################

  output$heat_map <- renderPlot({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select QC metrics and create a decision rule"),
      need(!is.null(prodata$AcquiredTime),"To view heatmaps, the dataset should include Acquired Time column.")
    )

    peptideThresholdRed <- (as.numeric(input$threshold_peptide_red))/100
    peptideThresholdYellow <- (as.numeric(input$threshold_peptide_yellow))/100
    if(is.null(prodata$AcquiredTime)) return(NULL)

    p1 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "XmR",
                           peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 1,
                           title = "Heatmap (Changes in mean of QC metric-X)")
    p2 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "XmR",
                           peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 2,
                           title = "Heatmap (Changes in variability of QC metric-mR)")
    p3 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "CUSUM",
                           peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 1,
                           title = "Heatmap (Changes in mean of QC metric-CUSUMm)")
    p4 <- metrics_heat.map(prodata,input$pepSelection,
                           data.metrics = input$user_selected_metrics, method = "CUSUM",
                           peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 2,
                           title = "Heatmap (Changes in variability of QC metric-CUSUMv)")
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
