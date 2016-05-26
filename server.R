library(shiny)
library(plotly)
library(RecordLinkage)

source("plot-functions.R")
source("data-validation.R")
source("helper-functions.R")

shinyServer(function(input,output,session) {
  
  #### Read data  ##################################################################################################
  prodata <- reactive({ # data is what user, upload in the app
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    file1 <- input$filein
    if(is.null(file1)){return()} 
    prodata <- read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE)
    input_checking(prodata)
    #prodata
    # read.xlsx2(file1$path , sheetName = "data")
  })
  
  ##### Precursor type selection #####################################################################################
  output$pepSelect <- renderUI({
    prodata <- prodata()
    selectInput("pepSelection","Choose precursor type", choices = c(levels(prodata$Precursor),"all peptides"))
  })
  ######Show data############################################
  
  output$prodata_table <- DT::renderDataTable(
    
    DT::datatable(prodata(), options = list(pageLength = 25))
  )
  ################################################################# plots ###################################################
  render.tab <- function(normalize.metric, plot.method, main.title, y.title1, y.title2){
    prodata <- prodata()
    plots <- list()
    
    if(input$pepSelection == "all peptides") {
      
      results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
        z <- normalize(prodata, j, input$L, input$U, metric = normalize.metric)
        plots[[2*j-1]] <<- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title1, 1)
        plots[[2*j]] <<- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title2, 2)
      })
      
      do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
        layout(autosize = F, width = 1000, height = nlevels(prodata$Precursor)*200)
    }
    
    else {
      j = which(levels(prodata$Precursor) == input$pepSelection)
      z <- normalize(prodata, j, input$L, input$U, metric = normalize.metric)
      
      plot1 <- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title1, 1)
      plot2 <- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title2, 2)
      
      subplot(plot1,plot2)
    }
  }
  ########################################################## plot CUSUM_chart  for RT####################
  output$RT_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "CUSUM", main.title = "Retention Time", y.title1 = "CUSUMm", y.title2 = "CUSUMv")
  })
  ########################################################## plot CP for RT #############################
  output$RT_CP <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "CP", main.title = "Retention Time", y.title1 = "Ci", y.title2 = "Di")
  })
  ####################################################### plot ZMR for RT ##############################
  output$RT_ZMR <- renderPlotly({
    render.tab(normalize.metric = "Retention Time", plot.method = "ZMR", main.title = "Retention Time", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################plot CUSUM for Peak assymetry ################
  output$PA_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "CUSUM", main.title = "Peak Assymetry", y.title1 = "CUSUMm", y.title2 = "CUSUMv")
  })
  ######################################################## plot Change Point for Peak assymetry ###########
  output$PA_CP <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "CP", main.title = "Peak Assymetry", y.title1 = "Ci", y.title2 = "Di")
  })
  ########################################################## plot ZMR for Peak assymetry ###################################
  output$PA_ZMR <- renderPlotly({
    render.tab(normalize.metric = "Peak Assymetry", plot.method = "ZMR", main.title = "Peak Assymetry", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################### plot CUSUM FOR Max.FWHM ####################################
  output$Max_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "CUSUM", main.title = "FWHM", y.title1 = "CUSUMm", y.title2 = "CUSUMv")    
  })
  ########################################################## plot Change Point FOR Max.FWHM ####################################
  output$Max_CP <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "CP", main.title = "FWHM", y.title1 = "Ci", y.title2 = "Di")    
  })
  ########################################################## plot ZMR FOR Max.FWHM ####################################
  output$Max_ZMR <- renderPlotly({
    render.tab(normalize.metric = "FWHM", plot.method = "ZMR", main.title = "FWHM", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ############################################################ plot CUSUM FOR total area ####################################
  output$TA_CUSUM <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "CUSUM", main.title = "Total Area", y.title1 = "CUSUMm", y.title2 = "CUSUMv")
  })
  ########################################################## plot Change Point FOR total area ##################################
  output$TA_CP <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "CP", main.title = "Total Area", y.title1 = "Ci", y.title2 = "Di")
  })
  ########################################################## plot ZMR FOR total area ##########################################
  output$TA_ZMR <- renderPlotly({
    render.tab(normalize.metric = "Total Area", plot.method = "ZMR", main.title = "Total Area", y.title1 = "Individual Value", y.title2 = "Moving Range")
  })
  ########################################################## box plot in Summary tab ##########################################
  output$box_plot <- renderPlotly({
    prodata <- prodata()
    metrics_box.plot(prodata)
  })
  ########################################################## scatterplot matrix in Summary tab #################################
  output$scatter_plot <- renderPlot({
    prodata <- prodata()
    metrics_scatter.plot(prodata, input$L, input$U, input$metric_precursor)
  }, height = 700)
  
  ###########################################################################################################################
  ###########################################################################################################################
  ########################################################## "help" tab ################################
  
  #### Text messages in empty places - This part will be removed in future, when the metrics codes is complete ######
  output$EWMA_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$Short_run_SPC_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$Multivariate_Control_Charts_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  # output$Capability_Analysis_txt <- renderText({
  #   paste0("This part is not complete yet, we will complete it in near future.")
  # })
  
  output$MA_ZMR_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$MA_CUSUM_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$CP_MA_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$CA_RT_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$CA_MA_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  
  output$OverallQC_txt <- renderText({
    paste0("This part is not complete yet, we will complete it in near future.")
  })
  ############################################################################################################################
})