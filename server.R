library(shiny)
library(plotly)
library(RecordLinkage)
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
  
  observeEvent(input$act_button, {
    data$df <- read.csv("./Datasets/Sampledata_CPTAC_Study_9_1_Site54.csv")
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
  # output$prodata_column_select <- renderUI({
  #   prodata <- data$df
  #   checkboxGroupInput("show_prodata_columns", "columns of your data", choices = colnames(prodata), selected = colnames(prodata))
  # })
  ######Show data#####################################################################################################
  output$prodata_table <- DT::renderDataTable(
    #DT::datatable(data$df[,ifelse(input$show_prodata_columns==NA,1, input$show_prodata_columns)], options = list(pageLength = 25))
    DT::datatable(data$df, options = list(pageLength = 25))
   
  )
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
  ########################################################## plot CUSUM_chart  for RT####################
  output$RT_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "CUSUM", normalization.type = TRUE, main.title = "Retention Time", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################## plot CP for RT #############################
  output$RT_CP <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "CP", normalization.type = TRUE, main.title = "Retention Time", y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ####################################################### plot ZMR for RT ##############################
  output$RT_ZMR <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "ZMR", normalization.type = FALSE, main.title = "Retention Time", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################plot CUSUM for Peak assymetry ################
  output$PA_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "CUSUM", normalization.type = TRUE, main.title = "Peak Assymetry", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ######################################################## plot Change Point for Peak assymetry ###########
  output$PA_CP <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "CP", normalization.type = TRUE, main.title = "Peak Assymetry", y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ########################################################## plot ZMR for Peak assymetry ###################################
  output$PA_ZMR <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "ZMR", normalization.type = FALSE, main.title = "Peak Assymetry", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################### plot CUSUM FOR MaxFWHM ####################################
  output$Max_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "CUSUM", normalization.type = TRUE, main.title = "FWHM", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")    
  })
  ########################################################## plot Change Point FOR MaxFWHM ####################################
  output$Max_CP <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "CP", normalization.type = TRUE, main.title = "FWHM", y.title1 = "Change point for mean", y.title2 = "Change point for variation")    
  })
  ########################################################## plot ZMR FOR MaxFWHM ####################################
  output$Max_ZMR <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "ZMR", normalization.type = FALSE, main.title = "FWHM", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ############################################################ plot CUSUM FOR total area ####################################
  output$TA_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "CUSUM", normalization.type = TRUE, main.title = "Total Area", y.title1 = "CUSUM mean", y.title2 = "CUSUM variation")
  })
  ########################################################## plot Change Point FOR total area ##################################
  output$TA_CP <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "CP", normalization.type = TRUE, main.title = "Total Area", y.title1 = "Change point for mean", y.title2 = "Change point for variation")
  })
  ########################################################## plot ZMR FOR total area ##########################################
  output$TA_ZMR <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "ZMR", normalization.type = FALSE, main.title = "Total Area", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## box plot in Summary tab ##########################################
  output$box_plot <- renderPlotly({
    prodata <- data$df
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
  output$plot_summary <- renderPlotly({
   
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    CUSUM.summary.plot(prodata, input$L, input$U,type = 1, ytitle = "probability of out of range points for CUSUM mean")
    
  })
  output$plot_summaryy <- renderPlotly({
    
    prodata <- data$df
    validate(
      need(!is.null(prodata), "Please upload your data")
    )
    CUSUM.summary.plot.version2(prodata, input$L, input$U,type = 1)
  })
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
  
  # output$MA_ZMR_txt <- renderText({
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