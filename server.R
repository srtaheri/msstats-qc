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
  ###### Tab for selecting decision rule and upper and lower bound decisions ###############################################
  output$decision_rule <- renderUI({
    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df)
    )

    #p("Please select your preferred decision rule: ")

    fluidPage(
      wellPanel(
      fluidRow(

        column(2,
               br(),
               "Good to go is when less than\n",
               br(),
               "Warning area is when, less than"
        ),
        column(2,
               numericInput('peptideThresholdGood', '', value = 50, min = 0,
                            max = 100, step = 1),
               numericInput('peptideThresholdWarn', '', value = 70, min = 1,
                            max = 100, step = 1)
               ),
        column(2,
               br(),
               "percent of peptides and\n",
               br(),br(),
               "and more than", "percentage of peptides and"
        ),
        column(2,
               numericInput('metricThresholdGood', '', value = 1, min = 1,
                            max = as.numeric(data$metrics), step = 1),
               numericInput('metricThresholdWarn', '', value = 2, min = 1,
                            max = as.numeric(data$metrics))
        ),
        column(3,
               br(),
               "of the selected metrics are out of control.\n",
               br(),br(),
               "of the selected metrics are out of control"
        )

      )),
      fluidRow(
        column(12,
               wellPanel(
                 checkboxGroupInput("user_selected_metrics","Please select the metrics",
                                    choices = c(data$metrics),
                                    selected = c(COL.PEAK.ASS,COL.BEST.RET,
                                                 COL.FWHM, COL.TOTAL.AREA),
                                    inline = TRUE)
               )


               )
      ),
      fluidRow(
        column(12,
               checkboxGroupInput("user_selected_type","Please select the type",
                                  choices = c("mean","dispersion"), selected = c("mean","dispersion"))
               )
      )

    )
  })

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
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata)
    )
    metrics_box.plot(prodata, data.metrics = data$metrics)
  })
  ########################################################## scatterplot matrix in Summary tab #################################
  # output$scatter_plot_metric_selection <- renderUI({
  #
  #    return(selectInput("scatter_metric_select", "Choose the metric",
  #                choices = c(data$metrics), selected = COL.FWHM))
  #
  # })
  # output$scatter_plot <- renderPlot({
  #   prodata <- data$df
  #     validate(
  #       need(!is.null(prodata), "Please upload your data"),
  #       need(is.data.frame(prodata), prodata)
  #     )
  #   metrics_scatter.plot(prodata,input$L, input$U, metric = input$scatter_metric_select, normalization = T)
  # }, height = 1500)
  # outputOptions(output, "scatter_plot_metric_selection", priority = 10)
  # outputOptions(output, "scatter_plot", priority = 1)
  ######################################################### plot_summary in Summary tab ########################################
  my_height <- reactive({
    my_height <- ceiling(length(data$metrics)/4)*1200
  })
  ################ plot the summary and radar plots ############################################################################
  output$plot_summary <- renderPlot({

    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata)
    )
    ##Text1 = textGrob("Mean",gp=gpar(col="gray40",fontsize = 11,face="bold"))
    ##Text2 = textGrob("Dispersion",gp=gpar(col="gray40",fontsize = 11,face="bold"))
    p1 <- XmR.Summary.plot(prodata, data.metrics = data$metrics, input$L, input$U)
    p2 <- XmR.Radar.Plot(prodata, data.metrics = data$metrics,input$L,input$U)
    p3 <- CUSUM.Summary.plot(prodata, data.metrics = data$metrics, input$L, input$U)
    p4 <- CUSUM.Radar.Plot(prodata, data.metrics = data$metrics, input$L,input$U)
    ##p1 <- p1 + annotation_custom(grob = Text1, xmin = 0, xmax = 7, ymin = 0, ymax = 2.6)
    ##p1 <- p1 + annotation_custom(grob = Text2, xmin = 0, xmax = 15, ymin = 0, ymax = -2.6)
    ##p3 <- p3 + annotation_custom(grob = Text1, xmin = 0, xmax = 7, ymin = 0, ymax = 2.6)
    ##p3 <- p3 + annotation_custom(grob = Text2, xmin = 0, xmax = 15, ymin = 0, ymax = -2.6)
    grid.arrange(p1,p2,p3,p4, ncol = 1)

  }, height = my_height )

  ################## decision message for XmR in summary tab #########
  # output$XmR_summary_decision_txt <- renderText({
  #   prodata <- data$df
  #   peptideThresholdGood <- (as.numeric(input$peptideThresholdGood))/100
  #   metricThresholdGood <- as.numeric(input$metricThresholdGood)
  #   peptideThresholdWarn <- (as.numeric(input$peptideThresholdWarn))/100
  #   metricThresholdWarn <- as.numeric(input$metricThresholdWarn)

   ####????? dorost konam
    # XmRCounterGood1 <- XmR.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThresholdGood,
    #                                                      input$L, input$U, type = 1)
    # XmRCounterGood2 <- XmR.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThresholdGood,
    #                                                      input$L, input$U, type = 2)
    # XmRCounterWarn1 <- XmR.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThresholdWarn,
    #                                                       input$L, input$U, type = 1)
    # XmRCounterWarn2 <- XmR.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThresholdWarn,
    #                                                       input$L, input$U, type = 2)

    # if(XmRCounterGood1 <= metricThresholdGood && XmRCounterGood2 <= metricThresholdGood) {"System is in-control"}
    # if(XmRCounterGood1 > metricThresholdGood && XmRCounterGood1 <= metricThresholdWarn &&
    #    XmRCounterGood2 <= metricThresholdGood){ "Warning! System is out-of-control (A change in QC metric mean is possible)"}
    # if(XmRCounterGood2 > metricThresholdGood && XmRCounterGood2 <= metricThresholdWarn &&
    #    XmRCounterGood1 <= metricThresholdGood){"Warning! System is out-of-control (A change in QC metric variation is possible)"}
    # if(XmRCounterGood1 > metricThresholdGood && XmRCounterGood1 <= metricThresholdWarn &&
    #    XmRCounterGood2 > metricThresholdGood && XmRCounterGood2 <= metricThresholdWarn) {"Warning! System is out-of-control (A simultaneous change in QC metric mean and variation is possible)"}
    # if(XmRCounterGood1 > metricThresholdWarn &&
    #    XmRCounterGood2 > metricThresholdGood && XmRCounterGood2 <= metricThresholdWarn) {"Bad! Mean is bad and variation is in warning area"}
    # if(XmRCounterGood2 > metricThresholdWarn &&
    #    XmRCounterGood1 > metricThresholdGood && XmRCounterGood1 <= metricThresholdWarn) {"Bad! variation is bad and mean is in warning area"}
    # if(XmRCounterGood1 > metricThresholdWarn && XmRCounterGood2 > metricThresholdWarn) {"Bad! both mean and variation are in bad area"}
  #})
  ###################
  # output$CUSUM_summary_decision_txt <- renderText({
  #   prodata <- data$df
  #   peptideThreshold <- (as.numeric(input$peptideThreshold))/100
  #   metricThreshold <- as.numeric(input$metricThreshold)
  #
  #   CUSUMCounter1 <- CUSUM.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThreshold,
  #                                                         input$L, input$U, type = 1)
  #   CUSUMCounter2 <- CUSUM.number.Of.Out.Of.Range.Metrics(prodata,data$metrics, peptideThreshold,
  #                                                         input$L, input$U, type = 2)
  #
  #   if(CUSUMCounter1 > metricThreshold){ "System is out-of-control (A change in QC metric mean is possible)"}
  #   if(CUSUMCounter2 > metricThreshold){"System is out-of-control (A change in QC metric variation is possible)"}
  #   if(CUSUMCounter1 > metricThreshold && CUSUMCounter2 > metricThreshold) {"System is out-of-control (A simultaneous change in QC metric mean and variation is possible)"}
  #
  # })

  #1500
  ############################# heat_map in Summary tab #############################################

  output$heat_map <- renderPlot({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata)
    )
    validate(
      need(!is.null(prodata$AcquiredTime),"To view heatmap, the data set should include AcquiredTime column.")
    )
    peptideThresholdGood <- (as.numeric(input$peptideThresholdGood))/100
    peptideThresholdWarn <- (as.numeric(input$peptideThresholdWarn))/100
    if(is.null(prodata$AcquiredTime)) return(NULL)
    XmR.heatmap.DataFrame.type1 <- XmR.heatmap.DataFrame(prodata,input$pepSelection,input$user_selected_metrics, peptideThresholdGood, peptideThresholdWarn, input$L, input$U, type = 1)
    #XmR.heatmap.DataFrame.type2 <- XmR.heatmap.DataFrame(prodata,input$pepSelection, input$L, input$U, type = 2)

    p1 <- metrics_heat.map(XmR.heatmap.DataFrame.type1,input$pepSelection, input$L, input$U, type = 1)
    #p2 <- metrics_heat.map(prodata,XmR.heatmap.DataFrame.type2,input$pepSelection, input$L, input$U, type = 2)
    #XmR.Decision.DataFrame.prepare(prodata, metric = "Retention Time", input$L, input$U,type = 1)
    #grid.arrange(p1,p2, ncol = 1)
    p1

  }
  )

  ############################################################################################################################
})
