options(shiny.maxRequestSize=100*1024^2)

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
  # COL.BEST.RET <- "Retention Time"
  # COL.FWHM <- "Full Width at Half Maximum"
  # COL.TOTAL.AREA <- "Total Peak Area"
  # COL.PEAK.ASS <- "Peak Assymetry"


  #### Read data  ##################################################################################################
  data <- reactiveValues(df = NULL, metrics = NULL)


  observeEvent(input$filein, {
    file1 <- input$filein
    data$df <- input_checking(read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE))
    validate(
      need(!is.null(data$df), "Please upload your data"),
      need(is.data.frame(data$df), data$df)
    )
    #data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
    data$metrics <- c(find_custom_metrics(data$df))
  }, priority = 20)

  observeEvent(input$sample_button, {
    data$df <- input_checking(read.csv("./Datasets/Sampledata_CPTAC_Study_9_1_Site54.csv"))
    validate(
      need(!is.null(data$df), "Please upload your data"),
      need(is.data.frame(data$df), data$df)
    )

    #data$metrics <- c(COL.BEST.RET, COL.TOTAL.AREA, COL.FWHM, COL.PEAK.ASS, find_custom_metrics(data$df))
    data$metrics <- c(find_custom_metrics(data$df))
  }, priority = 20)

  observeEvent(input$clear_button, {
    data$df <- NULL
    data$metrics <- NULL
  }, priority = 20)
  ##### Precursor type selection #####################################################################################
  output$pepSelect <- renderUI({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data.\n\n If your data contains min start time and max end time columns,the App will add a peak assymetry column automatically.\n\n Your data should contain a column named Annotation. Put all your metrics after this column.To see an example of a sample data click on the {Run with example data} button."),
      need(is.data.frame(prodata), prodata)
    )
    selectInput("pepSelection","Choose peptide"
                #,choices = c(levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET])),"all peptides")
                ,choices = c(levels(prodata$Precursor),"all peptides")
                )
  })
  ######Show table of data #####################################################################################################
   output$prodata_table <- renderDataTable({  
     validate(
       need(!is.null(data$df), "Please upload your data.\n\n If your data contains min start time and max end time columns, the App will add a peak assymetry column automatically.\n\n Your data should contain a column named Annotation. Put all your metrics after this column.To see an example of a sample data click on the {Run with example data} button."),
       need(is.data.frame(data$df), data$df)
     )
     data$df
   }, options = list(pageLength = 25))
###### selection tab in Data Improt and selection #####################################################
  output$selectMeanSD <- renderUI({
    lapply(input$user_selected_metrics,
           function(x){
             fluidRow(
               column(4,paste(x,":")),
               column(4,
                      numericInput(paste0("selectMean@",x),"mean",value = 1)
               ),
               column(4,
                      numericInput(paste0("selectSD@",x),"standard deviation",value = 1)
               )
             )
           })
  })


  output$selectGuideSet <- renderUI({
    fluidRow(
      column(6,
             numericInput("L","Lower bound of guide set",value = 1, min = 1, step = 1)
      ),
      column(6,
             numericInput("U","Upper bound of guide set", value = 5, min = 2, step = 1)
      )
    )


  })
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
    numOfMetrics <- length(input$user_selected_metrics)
    threshold_metric_red <- input$threshold_metric_red
     validate(
       need(!is.null(numOfMetrics),"loading..."),
       need(!is.null(threshold_metric_red),"loading...")
     )

    numericInput('threshold_metric_yellow', '', value = threshold_metric_red , min = 0, max = threshold_metric_red, step = 1)
  })

  output$metricSelection <- renderUI({
    checkboxGroupInput("user_selected_metrics","",
                       choices = c(data$metrics),
                       inline = TRUE)
  })
  
  output$metricSelectionErrorMsg <- renderUI({
  AfterannoColNum <- (which(colnames(data$df)=="Annotations")) + 1
  null_columns <- c()
  message <- c("This/these metrics contain NA:")
  if(AfterannoColNum < ncol(data$df)) {
    colNames <- colnames(data$df)
    for(i in  AfterannoColNum:ncol(data$df)) {
      if(sum(is.na(data$df[,i])) > 0) {
        null_columns <- c(null_columns,colNames[i])
      }
    }
  }
  message <- HTML(paste(strong(message)))
  for(i in 1:length(null_columns)) {
    message <- HTML(paste(message,null_columns[i], sep = '<br/>'))
  }
  if(length(message) > 0) {
    HTML(paste(message,strong("You cannot select this/these metrics for further analysis."), sep = '<br/>'))
  }
  
  })

  ################################################################# plots ###################################################
  #################################################################################################################
  output$XmR_tabset <- renderUI({

    validate(
      need(!is.null(data$df), "Please upload your data first"),
      need(is.data.frame(data$df), data$df),
      need(!is.null(input$user_selected_metrics),"Please first select metrics and create a decision rule")
    )
    is_guidset_selected <- FALSE
    if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      is_guidset_selected <- TRUE
    }
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
                                                             y.title1 = "Individual Value", y.title2 = "Moving Range",
                                                             selectMean = input[[paste0("selectMean@",x)]],selectSD = input[[paste0("selectSD@",x)]],
                                                             guidset_selected = is_guidset_selected)
                                             )

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
      need(!is.null(input$user_selected_metrics),"Please first select metrics and create a decision rule")
    )
    is_guidset_selected <- FALSE
    if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      is_guidset_selected <- TRUE
    }
    Tabs <- lapply(input$user_selected_metrics,
                   function(x) {

                       tabPanel(x,
                                tags$head(tags$style(type="text/css")),
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 tags$div("It may take a while to load the plots, please wait...",
                                                          id="loadmessage")),
                                renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L, input$U, metric = x, plot.method = "CUSUM", normalization = TRUE, y.title1 = "CUSUM mean", y.title2 = "CUSUM variation",selectMean = input[[paste0("selectMean@",x)]],selectSD = input[[paste0("selectSD@",x)]],guidset_selected = is_guidset_selected))
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
      need(!is.null(input$user_selected_metrics),"Please first select metrics and create a decision rule")
    )
    is_guidset_selected <- FALSE
    if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      is_guidset_selected <- TRUE
    }
    Tabs <- lapply(input$user_selected_metrics,
                   function(x) {
                       tabPanel(x,
                                tags$head(tags$style(type="text/css")),
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 tags$div("It may take a while to load the plots, please wait...",
                                                          id="loadmessage")),
                                renderPlotly(render.QC.chart(data$df, input$pepSelection, input$L, input$U, metric = x, plot.method = "CP", normalization = TRUE, y.title1 = "Change point for mean", y.title2 = "Change point for variation",selectMean = input[[paste0("selectMean@",x)]],selectSD = input[[paste0("selectSD@",x)]],guidset_selected = is_guidset_selected))
                                )
                   })

    do.call(tabsetPanel, Tabs)
  })
  ######################################################### height and width in Summary tab ########################################
  my_height <- reactive({
    l <- length(input$user_selected_metrics)
    k <- length(input$summary_controlChart_select)
    if(l < 5) {
      my_height <- ceiling(k)*700
    }else if(l < 9) {
      my_height <- ceiling(k)*1300
    }else if(l <15) {
      my_height <- ceiling(k)*1700
    }else if(l < 20) {
      my_height <- ceiling(k)*2000
    }else if(l < 25) {
      my_height <- ceiling(k)*2300
    }else {
      my_height <- ceiling(k)*2600
    }
  })
  
  my_width <- reactive({
    l = length(input$user_selected_metrics)
    if(l == 1) {
      my_width = 300
    }else if(l == 2) {
      my_width = 800
    }else if(l == 3) {
      my_width = 1200
    }
    my_width <- 1500
    
  })
  
  heatmap_height <- reactive({
    l <- length(input$user_selected_metrics)
    k <- length(input$heatmap_controlChart_select)
    if(l == 1) {
      heatmap_height <- ceiling(k)*300
    }else {
      heatmap_height <- ceiling(k)*ceiling(l)*200
    }
  })
  
  heatmap_width <- reactive({
    heatmap_width <- 1000
    print(heatmap_width)
  })
  ########################################################## box plot in Metric Summary tab ##########################################
  output$box_plot <- renderPlotly({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select metrics and create a decision rule")
    )
    metrics_box.plot(prodata, data.metrics = input$user_selected_metrics)
  })

  ###############   summary plots and radar plots ############################################################################
  output$plot_summary <- renderPlot({

    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select metrics and create a decision rule")
    )

    is_guidset_selected <- FALSE
    if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      is_guidset_selected <- TRUE
    }
    listMean <- list()
    listSD <- list()
    for(metric in input$user_selected_metrics){
      listMean[[metric]] <- input[[paste0("selectMean@",metric)]]
      listSD[[metric]] <- input[[paste0("selectSD@",metric)]]
    }

    plots <- list()
    i <- 1
    for(method in input$summary_controlChart_select) {
      p1 <- NULL
      p2 <- NULL
      if(method == "XmR") {
        p1 <- XmR.Summary.plot(prodata, data.metrics = input$user_selected_metrics, input$L, input$U, listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
        p2 <- XmR.Radar.Plot(prodata, data.metrics = input$user_selected_metrics,input$L,input$U,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
      } else if(method == "CUSUM") {
        p1 <- CUSUM.Summary.plot(prodata, data.metrics = input$user_selected_metrics, input$L, input$U,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
        p2 <- CUSUM.Radar.Plot(prodata, data.metrics = input$user_selected_metrics, input$L,input$U,listMean = listMean,listSD = listSD, guidset_selected = is_guidset_selected)
      }
      plots[[i]]   <- p1
      plots[[i+1]] <- p2

      i <- i+2
    }
    if(length(plots) > 0)
      do.call("grid.arrange", c(plots, ncol = 1))
  }, height = my_height, width = my_width)
  ############################# heat_map in Summary tab #############################################
  output$heat_map <- renderPlot({
    prodata <- data$df

    validate(
      need(!is.null(prodata), "Please upload your data"),
      need(is.data.frame(prodata), prodata),
      need(!is.null(input$user_selected_metrics),"Please first select metrics and create a decision rule"),
      need(!is.null(prodata$AcquiredTime),"To view heatmaps, the dataset should include Acquired Time column.")
    )

    peptideThresholdRed <- (as.numeric(input$threshold_peptide_red))/100
    peptideThresholdYellow <- (as.numeric(input$threshold_peptide_yellow))/100
    if(is.null(prodata$AcquiredTime)) return(NULL)

    is_guidset_selected <- FALSE
    if(input$selectGuideSetOrMeanSD == "Mean and standard deviation estimated from guide set") {
      is_guidset_selected <- TRUE
    }

    listMean <- list()
    listSD <- list()
    for(metric in input$user_selected_metrics){
     listMean[[metric]] <- input[[paste0("selectMean@",metric)]]
     listSD[[metric]] <- input[[paste0("selectSD@",metric)]]
    }

    plots <- list()
    i <- 1
    for(method in input$heatmap_controlChart_select) {
      p1 <- metrics_heat.map(prodata,
                             data.metrics = input$user_selected_metrics, method = method,
                             peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 1,
                             title = "Decision-map : mean",
                             listMean = listMean, listSD = listSD, guidset_selected = is_guidset_selected)
      p2 <- metrics_heat.map(prodata,
                             data.metrics = input$user_selected_metrics, method = method,
                             peptideThresholdRed, peptideThresholdYellow,input$L, input$U, type = 2,
                             title = "Decision-map : variability",
                             listMean = listMean, listSD = listSD, guidset_selected = is_guidset_selected)
      plots[[i]]   <- p1
      plots[[i+1]] <- p2

      i <- i+2
    }
    if(length(plots) > 0)
      do.call("grid.arrange", c(plots, ncol = 1))

  }, height = heatmap_height, width = heatmap_width)

  ############################################################################################################################
})
