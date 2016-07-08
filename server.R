library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)
library(RecordLinkage)
library(hash)

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
  
  render.tab <- function(normalize.metric, plot.method, normalization.type, main.title, y.title1, y.title2){
    
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    plots <- list()
    
    if(input$pepSelection == "all peptides") {
      
      results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
        z <- prepare_column(prodata, j, input$L, input$U, metric = normalize.metric, normalization = normalization.type)
        plots[[2*j-1]] <<- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title1, 1)
        plots[[2*j]] <<- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title2, 2)
      })
      
      do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
        layout(autosize = F, width = 1400, height = nlevels(prodata$Precursor)*200)
    }
    
    else {
      j = which(levels(reorder(prodata$Precursor,prodata$BestRetentionTime)) == input$pepSelection)
      z <- prepare_column(prodata, j, input$L, input$U, metric = normalize.metric, normalization = normalization.type)
      
      plot1 <- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title1, 1)
      plot2 <- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title2, 2)
      
      subplot(plot1,plot2)
    }
  }
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
    render.tab(normalize.metric = "Retention Time", plot.method = "XmR", normalization.type = FALSE, main.title = "Retention Time", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR for Peak assymetry #####################################################################################################
  output$PA_XmR <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "XmR", normalization.type = FALSE, main.title = "Peak Assymetry", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR FOR MaxFWHM ############################################################################################################
  output$Max_XmR <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "XmR", normalization.type = FALSE, main.title = "FWHM", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR FOR total area #########################################################################################################
  output$TA_XmR <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "XmR", normalization.type = FALSE, main.title = "Total Area", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR for metric1 #############################################################################################################
  output$m1 <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "XmR", normalization.type = FALSE, main.title = "Total Area", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot XmR for metric2 #############################################################################################################
  output$m2 <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "XmR", normalization.type = FALSE, main.title = "Total Area", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## plot CUSUM_chart  for RT##########################################################################################################
  output$RT_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "CUSUM", normalization.type = TRUE, main.title = "Retention Time", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################plot CUSUM for Peak assymetry ######################################################################################################
  output$PA_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "CUSUM", normalization.type = TRUE, main.title = "Peak Assymetry", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################### plot CUSUM FOR MaxFWHM #########################################################################################################
  output$Max_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "CUSUM", normalization.type = TRUE, main.title = "FWHM", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")    
  })
  ############################################################ plot CUSUM FOR total area ####################################
  output$TA_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "CUSUM", normalization.type = TRUE, main.title = "Total Area", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################## plot Change Point for RT #############################
  output$RT_CP <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "CP", normalization.type = TRUE, main.title = "Retention Time", y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ######################################################## plot Change Point for Peak assymetry ###########
  output$PA_CP <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "CP", normalization.type = TRUE, main.title = "Peak Assymetry", y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ########################################################## plot Change Point FOR MaxFWHM ####################################
  output$Max_CP <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "CP", normalization.type = TRUE, main.title = "FWHM", y.title1 = "Change point for mean", y.title2 = "Change point for variation")    
  })
  ########################################################## plot Change Point FOR total area ##################################
  output$TA_CP <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "CP", normalization.type = TRUE, main.title = "Total Area", y.title1 = "Change point for mean", y.title2 = "Change point for variation")
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
  output$plot_summary <- renderPlot({
    
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    #CUSUM.summary.plot(prodata, input$L, input$U,type = 1)
    #XmR.Summary.plot(prodata,input$L,input$U)
    CUSUM.Summary.plot(prodata, input$L, input$U)

  })
  
  output$plot_summary1 <- renderPlot({
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    multiplot(XmR.Radar.Plot(prodata,input$L,input$U),CUSUM.Radar.Plot(prodata,input$L,input$U))
    

  }, width = 1200, height = 1000)
  ###########################################################################################################################
  ###########################################################################################################################
  ########################################################## "help" tab ################################
  
  #### Text messages in empty places - This part will be removed in future, when the metrics codes is complete ######
  # output$EWMA_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$Short_run_SPC_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$Multivariate_Control_Charts_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  
  # output$Capability_Analysis_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  
  # output$MA_XmR_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$MA_CUSUM_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$CP_MA_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$CA_RT_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$CA_MA_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  # 
  # output$OverallQC_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  ############################################################################################################################
})