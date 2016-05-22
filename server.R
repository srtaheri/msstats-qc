library(ggplot2)
library(shiny)
library(geometry)
library(car)
library(qcc)
library(cowplot) # for plot_grid()
library(Rmisc) # for multiplot() function
library(rsconnect)
library(GGally)
library(RecordLinkage)
library(plotly)
library(gridExtra)
source("plot-functions.R")
source("data-validation.R")


shinyServer(function(input,output,session) {
  
#### Read data  ##################################################################################################
  prodata <- reactive({ # data is what user, upload in the app
     validate(
       need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
     )
    file1 <- input$filein
    if(is.null(file1)){return()} 
    input_checking(read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE))
    
    # read.xlsx2(file1$path , sheetName = "data")
  })
### Read sample data in "Functions and Sample Data" tab ###########################################################
  # sample_data <- reactive({
  #   read.csv('/Users/sarataheri/Desktop/updated\ ui\ and\ server/AutoQCdata.csv') # put your own path 
  # })
##### Precursor type selection #####################################################################################
  output$pep1 <- renderUI({
    prodata <- prodata()
    #selectInput("pep1.1","choose precursor type", choices = 1:nlevels(prodata$Precursor))
    selectInput("pep1.1","Choose precursor type", choices = c(levels(prodata$Precursor),"all peptides"))
  })
  
#### plot height ##################################################################################################
   my_height <- reactive({
     prodata <- prodata()
     if(input$pep1.1 == "all peptides") { 
       if(nlevels(prodata$Precursor) > 0 && nlevels(prodata$Precursor) < 6) {
         my_height <- 1000
        
       }
       else if(nlevels(prodata$Precursor) > 5 && nlevels(prodata$Precursor) <11) {
         my_height <- 2000
         
       }
       else if(nlevels(prodata$Precursor) > 10 && nlevels(prodata$Precursor) <16) {
         my_height <- 4000
        
       }
       else if(nlevels(prodata$Precursor) > 15 && nlevels(prodata$Precursor) < 21) {
         my_height <- 4800
         
       }
       else if(nlevels(prodata$Precursor) > 20 && nlevels(prodata$Precursor) < 26) {
         my_height <- 6000
       } 
       else {
         print("hello")
        
       }
       
     } # end first if
     
     else {  # if all peptide is not chosen
       my_height <- 300
       
     }
     
   })
#### title ########################################################################################################
  title <- reactive({
    prodata <- prodata()
    levels(prodata$Precursor) # names of precursors
  })
 #############################################################################################################################
  #################################################### Function  #############################################################
  ############################################################################################################################



#### Panel.cor for drawing scatter plot matrix ############################################################################
  panel.cor <- function(x, y, digits = 2, cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    # correlation coefficient
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste("r= ", txt, sep = "")
    text(0.5, 0.6, txt)
    
    # p-value calculation
    p <- cor.test(x, y)$p.value
    txt2 <- format(c(p, 0.123456789), digits = digits)[1]
    txt2 <- paste("p= ", txt2, sep = "")
    if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
    text(0.5, 0.4, txt2)
  }

###########################################################################################################################
###########################################################################################################################
################################################################# plots ###################################################
   normalize <- function(prodata, j, L, U, method) {
     precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
     x <- 0
     if(method == "Best.RT"){
       x = precursdata$Best.RT # raw data for retention time
     } else if(method == "Peak Assymetry") {
       x = precursdata$Max.End.Time-precursdata$Min.Start.Time # raw data for peak assymetry
     } else if(method == "FWHM") {
       x=precursdata$Max.FWHM
     } else if(method == "Total Area") {
       x=precursdata$Total.Area # raw data for total area
     } else {
       print("Error")
     }
     
     mu=mean(x[L:U]) # in-control process mean
     sd=sd(x[L:U]) # in-control process variance
     z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
     return(z)
   }
   
   render.tab <- function(normalize.method, plot.method, main.title, y.title1, y.title2){
     validate(
       need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
     )
     prodata <- prodata()
     plots <- list()
     
     if(input$pep1.1 == "all peptides") {
       
       results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
         z <- normalize(prodata, j, input$L, input$U, method = normalize.method)
         plots[[2*j-1]] <<- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title1, 1)
         plots[[2*j]] <<- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title2, 2)
       })
       
       do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
         layout(autosize = F, width = 1500, height = my_height())
     }
     
     else {
       j = which(levels(prodata$Precursor) == input$pep1.1)
       z <- normalize(prodata, j, input$L, input$U, method = normalize.method)
       
       plot1 <- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title1, 1)
       plot2 <- do.plot(prodata, z,j,input$L,input$U, method=plot.method, main.title, y.title2, 2)
       
       subplot(plot1,plot2) %>% layout(title = levels(prodata$Precursor)[j])
     }
   }
  ########################################################## plot CUSUM_chart  for RT####################
   output$RT_CUSUM <- renderPlotly({
     render.tab(normalize.method = "Best.RT", plot.method = "CUSUM", main.title = "Retention Time", y.title1 = "CUSUMm", y.title2 = "CUSUMv")
   })
   ########################################################## plot CP for RT #############################
   output$RT_CP <- renderPlotly({
     render.tab(normalize.method = "Best.RT", plot.method = "CP", main.title = "Retention Time", y.title1 = "Ci", y.title2 = "Di")
   })
  ####################################################### plot ZMR for RT ##############################
   output$RT_ZMR <- renderPlotly({
     render.tab(normalize.method = "Best.RT", plot.method = "ZMR", main.title = "Retention Time", y.title1 = "Individual Value", y.title2 = "Moving Range")
   })
  ########################################################plot CUSUM for Peak assymetry ################
   output$PA_CUSUM <- renderPlotly({
     render.tab(normalize.method = "Peak Assymetry", plot.method = "CUSUM", main.title = "Peak Assymetry", y.title1 = "CUSUMm", y.title2 = "CUSUMv")
   })
 ######################################################## plot Change Point for Peak assymetry ###########
   output$PA_CP <- renderPlotly({
     render.tab(normalize.method = "Peak Assymetry", plot.method = "CP", main.title = "Peak Assymetry", y.title1 = "Ci", y.title2 = "Di")
   })
 ########################################################## plot ZMR for Peak assymetry ###################################
   output$PA_ZMR <- renderPlotly({
     render.tab(normalize.method = "Peak Assymetry", plot.method = "ZMR", main.title = "Peak Assymetry", y.title1 = "Individual Value", y.title2 = "Moving Range")
   })
  ########################################################### plot CUSUM FOR Max.FWHM ####################################
   output$Max_CUSUM <- renderPlotly({
     render.tab(normalize.method = "FWHM", plot.method = "CUSUM", main.title = "FWHM", y.title1 = "CUSUMm", y.title2 = "CUSUMv")    
   })
########################################################## plot Change Point FOR Max.FWHM ####################################
   output$Max_CP <- renderPlotly({
     render.tab(normalize.method = "FWHM", plot.method = "CP", main.title = "FWHM", y.title1 = "Ci", y.title2 = "Di")    
   })
########################################################## plot ZMR FOR Max.FWHM ####################################
   output$Max_ZMR <- renderPlotly({
     render.tab(normalize.method = "FWHM", plot.method = "ZMR", main.title = "FWHM", y.title1 = "Individual Value", y.title2 = "Moving Range")
   })
############################################################ plot CUSUM FOR total area ####################################
   output$TA_CUSUM <- renderPlotly({
     render.tab(normalize.method = "Total Area", plot.method = "CUSUM", main.title = "Total Area", y.title1 = "CUSUMm", y.title2 = "CUSUMv")
   })
########################################################## plot Change Point FOR total area ##################################
   output$TA_CP <- renderPlotly({
     render.tab(normalize.method = "Total Area", plot.method = "CP", main.title = "Total Area", y.title1 = "Ci", y.title2 = "Di")
   })
########################################################## plot ZMR FOR total area ##########################################
   output$TA_ZMR <- renderPlotly({
     render.tab(normalize.method = "Total Area", plot.method = "ZMR", main.title = "Total Area", y.title1 = "Individual Value", y.title2 = "Moving Range")
   })
########################################################## box plot in Summary tab ##########################################
  output$box_plot <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    prodata <- prodata()
    input$act_button
    if (input$act_button == 0)
      return()
    
    prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$Best.RT) # to plot boxplots in decreasing order
    RT <- plot_ly(prodata, y = Best.RT, color = PrecursorRT, type = "box") %>% layout(showlegend = FALSE)
    
    prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$Max.End.Time - prodata$Min.Start.Time) # to plot boxplots in increasing order
    PA <- plot_ly(prodata, y = (Max.End.Time-Min.Start.Time), color = PrecursorPA, type = "box") %>% layout(showlegend = FALSE)
      #ylab("Peak Assymetry")+
                              
    prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
    prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$Total.Area) # to plot boxplots in decreasing order
    TPA <- plot_ly(prodata, y = Total.Area, color = PrecursorTA, type = "box") %>% layout(showlegend = FALSE)
      #ylab("Total Peak Area")+
      
    
    prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
    prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$Max.FWHM) 
    FWHM <- plot_ly(prodata, y = Max.FWHM, color = PrecursorFWHM, type = "box") %>% layout(showlegend = FALSE)
      #ylab("FWHM")+
      
    
    
  
    #isolate(plot_grid(RT, TPA, FWHM,PA,  ncol = 1, nrow = 4))
    subplot(RT, PA, TPA, FWHM, nrows = 4) 
 
    
  }
  #,height = 1500    # I set the height in ui.R for plotly
  )
########################################################## scatterplot matrix in Summary tab #################################
  output$scatter_plot <- renderPlot({
    
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
    prodata <- prodata()
    input$act_button
    if (input$act_button == 0)
      return()
    prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
    prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
    
    if(input$metric_precursor == "Peak Assymetry") {
      multidata<-matrix(0,length(prodata$Precursor),nlevels(prodata$Precursor))
      
      for (j in 1:nlevels(prodata$Precursor)) {
        precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
        x= precursdata$Max.End.Time-precursdata$Min.Start.Time # raw data for peak assymetry
        mu=mean(x[input$L:input$U]) # in-control process mean
        sd=sd(x[input$L:input$U]) # in-control process variance
        z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
        multidata[1:length(z),j]<-z
      }
      colnames(multidata) <- title()
      multidata=data.frame(multidata)
      isolate(pairs(multidata, upper.panel = panel.cor, col = "blue"))
    } else {   
      multidata<-matrix(0,length(prodata$Precursor),nlevels(prodata$Precursor))
      
      for (j in 1:nlevels(prodata$Precursor)) {
        precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
        x=precursdata[,which(colnames(prodata) == input$metric_precursor)]  # raw data for total area
        mu=mean(x[input$L:input$U]) # in-control process mean
        sd=sd(x[input$L:input$U]) # in-control process variance
        z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
        multidata[1:length(z),j]<-z
      }
      colnames(multidata) <- title()
      multidata=data.frame(multidata)
      isolate(pairs(multidata, upper.panel = panel.cor, col = "blue"))
      #pairs(multidata, upper.panel = panel.cor, col = "blue")
    }
  }
  , height = 1000
  )

###########################################################################################################################
###########################################################################################################################
########################################################## "help" tab ################################


########################################################## sample data set in "help" tab ###############
  
########################################################## upload video ####################################
output$video <- renderUI({
  tags$video(src='reactive.mp4', type="video/mp4", width="800px", 
             height="800px", controls='controls')
})
##############################################################################################################   

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
   
   output$Capability_Analysis_txt <- renderText({
     paste0("This part is not complete yet, we will complete it in near future.")
   })
   
   output$MA_ZMR_txt <- renderText({
     paste0("This part is not complete yet, we will complete it in near future.")
   })
   
   output$MA_CUSUM_txt <- renderText({
     paste0("This part is not complete yet, we will complete it in near future.")
   })
   
   output$CP_MA_txt <- renderText({
     paste0("This part is not complete yet, we will complete it in near future.")
   })
   
   output$CA_MA_txt <- renderText({
     paste0("This part is not complete yet, we will complete it in near future.")
   })

############################################################################################################################
})