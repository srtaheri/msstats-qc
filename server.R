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
source("CUSUM-Plot-functions.R")
source("CP-Plot-functions.R")

shinyServer(function(input,output,session){
  
#### Read data  ##################################################################################################
  prodata <- reactive({ # data is what user, upload in the app
     validate(
       need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
     )
    file1 <- input$filein
    if(is.null(file1)){return()} 
    Input_checking(read.csv(file=file1$datapath, sep=",", header=TRUE, stringsAsFactors=TRUE))
    
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
#### plotly trace titles ##########################################################################################
  #trace_title 
#### Best selection for column names of uploaded data sets ########################################################
   best_colnames <- reactive({
     # here we put a selection of most column names that users use. The first element of each vector should be the best name that
     # we suggest users to use and  which our code is based on. for example "Best.RT" and "Max FWHM" which are the first element
     # of each vector in the list, are our suggestion so we wrote them in the fisrt place.
     list(
          c("Best.RT","best retention time", "retention time","rt","best ret time","intensity"),
          c("Max.FWHM","fwhm")

          )
   })
########################################
   output$activeTab <- reactive({
     return(input$tab)
   })
   outputOptions(output, 'activeTab', suspendWhenHidden=FALSE)

 #############################################################################################################################
  #################################################### Function  #############################################################
  ############################################################################################################################

#### IMR plot 1 ###########################################################################################################
  IMR_plot1 <- function(z,j,L,U,Main.title) {
    
    #title <- title()
    prodata <- prodata()
    MR<- numeric(length(z)-1)
    ###################
    
    ## Calculate X chart statistics and limits 
    UCLI=3
    LCLI=-3
    #UCLI=MRmean+2.66*sd(MR[L:U])
    #LCLI=MRmean-2.66*sd(MR[L:U])
    QCno=1:length(z)
    plot.data=data.frame(QCno,z,UCLI,LCLI)
    y.max=ifelse(max(plot.data$z)>=UCLI,(max(plot.data$z)),UCLI)
    y.min=ifelse(min(plot.data$z)<=LCLI,(min(plot.data$z)),LCLI)
    #plot.subtitle=title[j]
    Main=Main.title
    #####################
    x <- list(
      title = paste("QCno - ", levels(prodata$Precursor)[j])
    )
    y <- list(
      title = "Individual Value"
    )
    
    plot_ly(plot.data, x = QCno, y = z
            ,type = "scatter"
            ,line = list(shape = "linear")
            ,marker=list(color="dodgerblue" , size=4 , opacity=0.5)
            ,showlegend = FALSE
            ) %>%
      layout(xaxis = x,yaxis = y) %>%
      add_trace( y = UCLI, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines",showlegend = FALSE) %>%
      add_trace(y = LCLI, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines",showlegend = FALSE) %>%
      add_trace(x = plot.data[z <= LCLI, ]$QCno, y = plot.data[z <= LCLI, ]$z
                , mode = "markers"
                , marker=list(color="red" , size=8 , opacity=0.5)
                ,showlegend = FALSE
                ) %>%
      add_trace(x = plot.data[z >= UCLI, ]$QCno, y = plot.data[z >= UCLI, ]$z
                , mode = "markers"
                , marker=list(color="red" , size=8 , opacity=0.5)
                ,showlegend = FALSE
                ) %>%
      add_trace(x = plot.data[z > LCLI & z < UCLI, ]$QCno, y = plot.data[z > LCLI & z < UCLI, ]$z
                , mode = "markers"
                , marker=list(color="blue" , size=8 , opacity=0.5)
                ,showlegend = FALSE
                )
  }
#### IMR plot 2 ###########################################################################################################
  IMR_plot2 <- function(z,j,L,U,Main.title) {
    
    #title <- title()
    prodata <- prodata()
    MR<- numeric(length(z)-1)
    ###################
    
    ## Calculate X chart statistics and limits 
    
    ##UCLI=MRmean+2.66*sd(MR[L:U])
    ##LCLI=MRmean-2.66*sd(MR[L:U])

    #plot.subtitle=title[j]
    Main=Main.title
    #####################
    #####################
    ## Calculate MR chart statistics and limits
    for(i in 2:length(z))
    {
      MR[i]=abs(z[i]-z[i-1]) # Compute moving range of z
    }
    UCLMR=3.267*sd(MR[1:L-U])
    LCLMR=0
    QCno=1:length(z)
    plot.data=data.frame(QCno,z,MR,UCLMR,LCLMR)
    ymax=ifelse(max(plot.data$MR)>=UCLMR,(max(plot.data$MR)),UCLMR)
    ymin=ifelse(min(plot.data$MR)<=LCLMR,(min(plot.data$MR)),LCLMR)
    ######################
    
    x <- list(
      title = paste("QCno - ", levels(prodata$Precursor)[j])
    )
    y <- list(
      title = "Moving Range"
    )
    
    plot_ly(plot.data, x = QCno, y = MR, type = "scatter",
            name = "linear",  line = list(shape = "linear"),
            marker=list(color="dodgerblue" , size=4 , opacity=0.5)) %>%
      layout(xaxis = x,yaxis = y) %>%
      add_trace( y = UCLMR, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines") %>%
      add_trace(y = LCLMR, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines") %>%
      add_trace(x = plot.data[MR <= LCLMR, ]$QCno, y = plot.data[MR <= LCLMR, ]$MR, mode = "markers", marker=list(color="red" , size=8 , opacity=0.5)) %>%
      add_trace(x = plot.data[MR >= UCLMR, ]$QCno, y = plot.data[MR >= UCLMR, ]$MR, mode = "markers", marker=list(color="red" , size=8 , opacity=0.5)) %>%
      add_trace(x = plot.data[MR > LCLMR & MR < UCLMR, ]$QCno, y = plot.data[MR > LCLMR & MR < UCLMR, ]$MR, mode = "markers", marker=list(color="blue" , size=8 , opacity=0.5))
  }
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
#### camelCaseSplit function ##############################################################################################
  camelCaseSplit <- function(x) {
    # This function get a camelCase word and splits it.
    # Ex : camelCaseSplit("myComputerIsHere") ---> my Computer Is Here
    return(gsub("([a-z])([A-Z])", "\\1 \\L\\2", x, perl = TRUE))
  }
#### punc_remove function #################################################################################################
   punc_remove <- function(x){
     # This function removes any existing punctuation in your sentence or word and transfer it to space.
     # Ex1: punc_remove(Best.RT) --> Best RT     #Ex2: punc_remove(Best_RT) --> Best RT 
    return(gsub("[[:punct:]///' ]", " ", x))
   }
#### clearString function ###############################################################################################
    clearString <- function(x){
      # This function, gets a word or setence, Splits it (if it is a camelCase), removes any existing punctuations, and transfer
      # all Upper Case letters to lower case letters.
      # Ex: clearString("myName_isSara.Taheri") --> my name is sara taheri
     return(tolower(punc_remove(camelCaseSplit(x))))
   }
#### guessColumnName function ###########################################################################################
    guessColumnName <- function(x){
      best_colnames <- best_colnames()
      # This function receives the data and check the column names of data and changes the column names if it is not the
      # same names as our suggested sample data to fit our suggested sample data.
      a <- clearString(x)
      
      max_index <- 0
      max <- -1
      for(i in 1:length(best_colnames)){
        col <- best_colnames[[i]]
        for(j in 1:length(col)){
          sim <- levenshteinSim(a,col[j])
          if(sim > max){
            max <- sim
            max_index <- i
          }
        }
      }
      if (max > 0.6) {
        return(best_colnames[[max_index]][1])
      } 
      else {
        return(x)
      }
    }
### Input_checking function #########################################################################################
    Input_checking <- function(data){
      data[data==""] <- NA
      colnames(data) <- unlist(lapply(colnames(data), function(x)guessColumnName(x)))
      #print(is.logical(colnames(data)[c(1,2)] == c("Precursor","Sample.File")) )
      #print(colnames(sample_data()))
      #print(colnames(data))
      return(data)
    }
    ####
###########################################################################################################################
###########################################################################################################################
################################################################# plots ###################################################
   
  ########################################################## plot CUSUM_chart  for RT####################
  output$RT_CUSUM <- renderPlotly({


     validate(
       need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
     )

     prodata <- prodata()
     plots <- list()
     #input$act_button   # if user doesn't press "click to see plots" button, nothing is shown
     #prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$Best.RT) # order precursors in increasing order based on RT

     #if (input$act_button == 0)
        #return()

     
     if(input$pep1.1 == "all peptides") {
       #input$act_button
       
       results <- lapply(c(1:nlevels(prodata$Precursor)), function(j){
         precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
         x=precursdata$Best.RT # raw data for retention time
         mu=mean(x[input$L:input$U]) # in-control process mean
         sd=sd(x[input$L:input$U]) # in-control process variance
         z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
         plots[[2*j-1]] <<- CUSUM_plot(prodata, z,j,input$L,input$U,"Retention Time", "CUSUMm", 1) #CUSUM_plot1(prodata, z,j,input$L,input$U,"Retention Time")
         plots[[2*j]] <<- CUSUM_plot(prodata, z,j,input$L,input$U,"Retention Time", "CUSUMv", 2)  #CUSUM_plot2(prodata, z,j,input$L,input$U,"Retention Time")
       })
      
       do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
         layout(autosize = F, width = 1500, height = my_height())
     }

     else{
       #input$act_button
     j = which(levels(prodata$Precursor) == input$pep1.1)
     precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
     x=precursdata$Best.RT # raw data for retention time
     mu=mean(x[input$L:input$U]) # in-control process mean
     sd=sd(x[input$L:input$U]) # in-control process variance
     z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
    
      plot1 <- CUSUM_plot(prodata, z,j,input$L,input$U,"Retention Time", "CUSUMm", 1) #CUSUM_plot1(prodata, z,j,input$L,input$U, "Retention Time")
      plot2 <- CUSUM_plot(prodata, z,j,input$L,input$U,"Retention Time", "CUSUMv", 2) #CUSUM_plot2(prodata, z,j,input$L,input$U, "Retention Time")
      
      subplot(plot1,plot2) %>% layout(title = levels(prodata$Precursor)[j])
     }

   }
  
   )
  ########################################################## plot CP for RT #############################
  output$RT_CP <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
    prodata <- prodata()
    plots <- list()
    #prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$Best.RT) # order precursors in increasing order based on RT
    #input$act_button
    #if (input$act_button == 0)
    #  return()
    
    if(input$pep1.1 == "all peptides") {
      
      results <- lapply(c(1:nlevels(prodata$Precursor)), function(j){
        precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
        x=precursdata$Best.RT # raw data for retention time
        mu=mean(x[input$L:input$U]) # in-control process mean
        sd=sd(x[input$L:input$U]) # in-control process variance
        z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
        plots[[2*j-1]] <<- CP_plot(prodata, z,j,"Retention Time",1,"Ci")
        plots[[2*j]] <<- CP_plot(prodata, z,j,"Retention Time", 2, "Di")
      })

      do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
        layout(autosize = F, width = 1500, height = my_height())
    }
    else{
    
    j = which(levels(prodata$Precursor) == input$pep1.1)
    precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
    x=precursdata$Best.RT # raw data for retention time
    mu=mean(x[input$L:input$U]) # in-control process mean
    sd=sd(x[input$L:input$U]) # in-control process variance
    z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
    
    plot1 <- CP_plot(prodata, z,j,"Retention Time",1,"Ci")
    plot2 <- CP_plot(prodata, z,j,"Retention Time", 2, "Di")
    subplot(plot1,plot2)
    }
    }
    )
  ####################################################### plot ZMR for RT ##############################
  output$RT_ZMR <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
        prodata <- prodata()
        plots <- list()
        #prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$Best.RT) # order precursors in increasing order based on RT
        input$act_button
        if (input$act_button == 0)
          return()
        
        if(input$pep1.1 == "all peptides") {
          
          
          
          # for (j in 1:nlevels(prodata$Precursor)) {
          #   precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorRT)[j],] # subset for a particular precursor
          #   x=precursdata$Best.RT # raw data for retention time
          #   mu=mean(x[input$L:input$U]) # in-control process mean
          #   sd=sd(x[input$L:input$U]) # in-control process variance
          #   z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
          #   plots[[2*j-1]] <- IMR_plot1(z,j,input$L,input$U,"Retention Time") 
          #   #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorRT)[j]), "")))
          #   plots[[2*j]] <- IMR_plot2(z,j,input$L,input$U,"Retention Time") 
          #   #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorRT)[j]), "")))
          #   
          # }
          results <- lapply(c(1:nlevels(prodata$Precursor)), function(j){
            precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
            x=precursdata$Best.RT # raw data for retention time
            mu=mean(x[input$L:input$U]) # in-control process mean
            sd=sd(x[input$L:input$U]) # in-control process variance
            z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
            plots[[2*j-1]] <<- IMR_plot1(z,j,input$L,input$U,"Retention Time") 
            plots[[2*j]] <<- IMR_plot2(z,j,input$L,input$U,"Retention Time") 
          })
          #number_of_plots <- 2*nlevels(prodata$Precursor)
          #layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
          #isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
          do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
            layout(autosize = F, width = 1500, height = my_height())
        }
        else{
        #j = as.numeric(input$pep1.1)
        j = which(levels(prodata$Precursor) == input$pep1.1)
        precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
        x=precursdata$Best.RT # raw data for retention time
        mu=mean(x[input$L:input$U]) # in-control process mean
        sd=sd(x[input$L:input$U]) # in-control process variance
        z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
        
        plot1 <- IMR_plot1(z,j,input$L,input$U, "Retention Time")
        plot2 <- IMR_plot2(z,j,input$L,input$U, "Retention Time")
        subplot(plot1,plot2)
        }
      }
    #, height = my_height
    )
  ########################################################plot CUSUM for Peak assymetry ################
  output$PA_CUSUM <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
         prodata <- prodata()
         plots <- list()
         #prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$Max.End.Time - prodata$Min.Start.Time)
         input$act_button
         if (input$act_button == 0)
           return()
         
         if(input$pep1.1 == "all peptides") {
           
           results <- lapply(c(1:nlevels(prodata$Precursor)), function(j){
             precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
             x=precursdata$Max.End.Time-precursdata$Min.Start.Time # raw data for peak assymetry
             mu=mean(x[input$L:input$U]) # in-control process mean
             sd=sd(x[input$L:input$U]) # in-control process variance
             z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
             plots[[2*j-1]] <<- CUSUM_plot1(z,j,input$L,input$U,"Peak Assymetry")
             plots[[2*j]] <<- CUSUM_plot2(z,j,input$L,input$U,"Peak Assymetry")
           })

           do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% layout(autosize = F, width = 1500, height = my_height())
         }
         else{
         j = which(levels(prodata$Precursor) == input$pep1.1)
         precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
         x=precursdata$Max.End.Time-precursdata$Min.Start.Time  # raw data for peak assymetry               
         mu=mean(x[input$L:input$U]);
         sd=sd(x[input$L:input$U]);
         z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
         
         plot1 <- CUSUM_plot1(z,j,input$L,input$U, "Peak Assymetry")
         plot2 <- CUSUM_plot2(z,j,input$L,input$U, "Peak Assymetry")
         subplot(plot1,plot2)
         }
       }
    )
 ######################################################## plot Change Point for Peak assymetry ###########
  output$PA_CP <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
        prodata <- prodata()
        prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$Max.End.Time - prodata$Min.Start.Time)
        input$act_button
        if (input$act_button == 0)
          return()
        
        if(input$pep1.1 == "all peptides") {
          
          plots <- list()
          
          for (j in 1:nlevels(prodata$Precursor)) {
            precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorPA)[j],] # subset for a particular precursor
            x=precursdata$Max.End.Time-precursdata$Min.Start.Time # raw data for peak assymetry 
            mu=mean(x[input$L:input$U]) # in-control process mean
            sd=sd(x[input$L:input$U]) # in-control process variance
            z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
            plots[[2*j-1]] <- CP_plot1(z,j,"Peak Assymetry") 
            #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorPA)[j]), "")))
            plots[[2*j]] <- CP_plot2(z,j,"Peak Assymetry") 
            #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorPA)[j]), "")))
            
          }
          # number_of_plots <- 2*nlevels(prodata$Precursor)
          # layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
          # isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
          do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
            layout(autosize = F, width = 1500, height = my_height())
          
        }
        else{
        #j = as.numeric(input$pep1.1)
        j = which(levels(prodata$Precursor) == input$pep1.1)
        precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
        x=precursdata$Max.End.Time-precursdata$Min.Start.Time  # raw data for peak assymetry             
        mu=mean(x[input$L:input$U]);
        sd=sd(x[input$L:input$U]);
        z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
         
        plot1 <- CP_plot1(z,j, "Peak Assymetry")
        plot2 <- CP_plot2(z,j, "Peak Assymetry")
        subplot(plot1,plot2)
        }
       }
    #, height = my_height
    )
 ########################################################## plot ZMR for Peak assymetry ###################################
  output$PA_ZMR <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
         prodata <- prodata()
         prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$Max.End.Time - prodata$Min.Start.Time)
         input$act_button
         if (input$act_button == 0)
           return()
         
         if(input$pep1.1 == "all peptides") {
           
           plots <- list()
           
           for (j in 1:nlevels(prodata$Precursor)) {
             precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorPA)[j],] # subset for a particular precursor
             x=precursdata$Max.End.Time-precursdata$Min.Start.Time  # raw data for peak assymetry
             mu=mean(x[input$L:input$U]) # in-control process mean
             sd=sd(x[input$L:input$U]) # in-control process variance
             z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
             plots[[2*j-1]] <- IMR_plot1(z,j,input$L,input$U,"Peak Assymetry") 
             #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorPA)[j]), "")))
             plots[[2*j]] <- IMR_plot2(z,j,input$L,input$U,"Peak Assymetry") 
             #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorPA)[j]), "")))
             
           }
           # number_of_plots <- 2*nlevels(prodata$Precursor)
           # layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
           # isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
           do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
             layout(autosize = F, width = 1500, height = my_height())
         }
         else{
         #j = as.numeric(input$pep1.1)
         j = which(levels(prodata$Precursor) == input$pep1.1)
         precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
         x=precursdata$Max.End.Time-precursdata$Min.Start.Time  # raw data for peak assymetry             
         mu=mean(x[input$L:input$U]);
         sd=sd(x[input$L:input$U]);
         z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
         
         plot1 <- IMR_plot1(z,j,input$L,input$U, "Peak Assymetry")
         plot2 <- IMR_plot2(z,j,input$L,input$U, "Peak Assymetry")
         subplot(plot1,plot2)
         }
        }
    #, height = my_height
    )
  ########################################################### plot CUSUM FOR Max.FWHM ####################################
  output$Max_CUSUM <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
        
          prodata <- prodata()
          plots <- list()
          prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
          #prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$Max.FWHM) 
          input$act_button
          
          if (input$act_button == 0)
            return()
          
          if(input$pep1.1 == "all peptides") {
            
            results <- lapply(c(1:nlevels(prodata$Precursor)), function(j){
              precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
              x=precursdata$Max.FWHM # raw data for fwhm
              mu=mean(x[input$L:input$U]) # in-control process mean
              sd=sd(x[input$L:input$U]) # in-control process variance
              z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
              plots[[2*j-1]] <<- CUSUM_plot1(z,j,input$L,input$U,"FWHM")
              plots[[2*j]] <<- CUSUM_plot2(z,j,input$L,input$U,"FWHM")
            })
            
            do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% layout(autosize = F, width = 1500, height = my_height())
          }
          else{
          
          j = which(levels(prodata$Precursor) == input$pep1.1)
          precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
          x=precursdata$Max.FWHM; # raw data for fwhm
          mu=mean(x[input$L:input$U]); # in-control process mean
          sd=sd(x[input$L:input$U]); # in-control process variance
          z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
          
          plot1 <- CUSUM_plot1(z,j,input$L,input$U, "Peak Assymetry")
          plot2 <- CUSUM_plot2(z,j,input$L,input$U, "Peak Assymetry")
          subplot(plot1,plot2)
          }
        }
    )
########################################################## plot Change Point FOR Max.FWHM ####################################
  output$Max_CP <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
          prodata <- prodata()
          prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
          prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$Max.FWHM) 
          input$act_button
          if (input$act_button == 0)
            return()
          
          if(input$pep1.1 == "all peptides") {
            
            plots <- list()
            
            for (j in 1:nlevels(prodata$Precursor)) {
              precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorFWHM)[j],] # subset for a particular precursor
              x=precursdata$Max.FWHM # raw data for fwhm
              mu=mean(x[input$L:input$U]) # in-control process mean
              sd=sd(x[input$L:input$U]) # in-control process variance
              z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
              plots[[2*j-1]] <- CP_plot1(z,j,"FWHM") 
              #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorFWHM)[j]), "")))
              plots[[2*j]] <- CP_plot2(z,j,"FWHM") 
              #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorFWHM)[j]), "")))
              
            }
            # number_of_plots <- 2*nlevels(prodata$Precursor)
            # layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
            # isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
            do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
              layout(autosize = F, width = 1500, height = my_height())
          }
          else{
          #j = as.numeric(input$pep1.1)
          j = which(levels(prodata$Precursor) == input$pep1.1)
          #prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
          precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
          x=precursdata$Max.FWHM; # raw data for fwhm
          mu=mean(x[input$L:input$U]); # in-control process mean
          sd=sd(x[input$L:input$U]); # in-control process variance
          z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
          
          plot1 <- CP_plot1(z,j, "FWHM")
          plot2 <- CP_plot2(z,j, "FWHM")
          subplot(plot1,plot2)
          }
        }
    #, height = my_height
    )
########################################################## plot ZMR FOR Max.FWHM ####################################
  output$Max_ZMR <- renderPlotly({
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
    
          prodata <- prodata()
          prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
          prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$Max.FWHM) 
          input$act_button
          if (input$act_button == 0)
            return()
         
          
          if(input$pep1.1 == "all peptides") {
            
            
            plots <- list()
            for (j in 1:nlevels(prodata$Precursor)) {
              precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorFWHM)[j],] # subset for a particular precursor
              x=precursdata$Max.FWHM  # raw data for fwhm
              mu=mean(x[input$L:input$U]) # in-control process mean
              sd=sd(x[input$L:input$U]) # in-control process variance
              z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
              plots[[2*j-1]] <- IMR_plot1(z,j,input$L,input$U,"FWHM") 
              #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorFWHM)[j]), "")))
              plots[[2*j]] <- IMR_plot2(z,j,input$L,input$U,"FWHM") 
              #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorFWHM)[j]), "")))
              
            }
            # number_of_plots <- 2*nlevels(prodata$Precursor)
            # layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
            # isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
            do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
              layout(autosize = F, width = 1500, height = my_height())
          }
          else{
          #j = as.numeric(input$pep1.1)
          j = which(levels(prodata$Precursor) == input$pep1.1)
          #prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
          precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
          x=precursdata$Max.FWHM; # raw data for fwhm
          mu=mean(x[input$L:input$U]); # in-control process mean
          sd=sd(x[input$L:input$U]); # in-control process variance
          z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
          
          plot1 <- IMR_plot1(z,j,input$L,input$U, "FWHM")
          plot2 <- IMR_plot2(z,j,input$L,input$U, "FWHM")
          subplot(plot1,plot2)
          }
        }
    #, height = my_height
    )
############################################################ plot CUSUM FOR total area ####################################
  output$TA_CUSUM <- renderPlotly({
    
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
           prodata <- prodata()
           plots <- list()
           prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
           #prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$Total.Area) # order precursors in increasing order based on Total Area
           input$act_button
           
           if (input$act_button == 0)
             return()
           
           if(input$pep1.1 == "all peptides") {
             
             results <- lapply(c(1:nlevels(prodata$Precursor)), function(j){
               precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
               x=precursdata$Total.Area # raw data for total area
               mu=mean(x[input$L:input$U]) # in-control process mean
               sd=sd(x[input$L:input$U]) # in-control process variance
               z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
               plots[[2*j-1]] <<- CUSUM_plot1(z,j,input$L,input$U,"Total Area")
               plots[[2*j]] <<- CUSUM_plot2(z,j,input$L,input$U,"Total Area")
             })

             do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% layout(autosize = F, width = 1500, height = my_height())
             
           }
           else{
           j = which(levels(prodata$Precursor) == input$pep1.1)
           precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
           x=precursdata$Total.Area; # raw data for total area
           mu=mean(x[input$L:input$U]); # in-control process mean
           sd=sd(x[input$L:input$U]); # in-control process variance
           z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
           
           plot1 <- CUSUM_plot1(z,j,input$L,input$U, "Peak Assymetry")
           plot2 <- CUSUM_plot2(z,j,input$L,input$U, "Peak Assymetry")
           subplot(plot1,plot2)
           }
         }
    )
########################################################## plot Change Point FOR total area ##################################
  output$TA_CP <- renderPlotly({
    
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
           prodata <- prodata()
           input$act_button
           prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
           prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$Total.Area) # order precursors in increasing order based on Total Area
           
           if (input$act_button == 0)
             return()
           
          
            if(input$pep1.1 == "all peptides") {
             
             plots <- list()
             #prodata$Max.FWHM <- as.numeric(gsub(",","",prodata$Max.FWHM))
             for (j in 1:nlevels(prodata$Precursor)) {
               precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorTA)[j],] # subset for a particular precursor
               x=precursdata$Total.Area # raw data for total area
               mu=mean(x[input$L:input$U]) # in-control process mean
               sd=sd(x[input$L:input$U]) # in-control process variance
               z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
               plots[[2*j-1]] <- CP_plot1(z,j,"Total Area") 
               #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorTA)[j]), "")))
               plots[[2*j]] <- CP_plot2(z,j,"Total Area") 
               #+  ggtitle(bquote(atop(.(levels(prodata$PrecursorTA)[j]), "")))
               
             }
             # number_of_plots <- 2*nlevels(prodata$Precursor)
             # layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
             # isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
             do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
               layout(autosize = F, width = 1500, height = my_height())
           }
           else{
           #j = as.numeric(input$pep1.1)
           j = which(levels(prodata$Precursor) == input$pep1.1)
           #prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
           precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
           x=precursdata$Total.Area; # raw data for total area
           #Main.title="Total Area" # main title of each plot
           mu=mean(x[input$L:input$U]); # in-control process mean
           sd=sd(x[input$L:input$U]); # in-control process variance
           z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
        
           plot1 <- CP_plot1(z,j, "Total Area")
           plot2 <- CP_plot2(z,j, "Total Area")
           subplot(plot1,plot2)
           }
         }
    #, height = my_height
    )
########################################################## plot ZMR FOR total area ##########################################
  output$TA_ZMR <- renderPlotly({
    
    validate(
      need(!is.null(input$filein), "No plot is shown here, because you have not uploaded your data")
    )
           prodata <- prodata()
           prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
           input$act_button
           prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$Total.Area) # order precursors in increasing order based on Total Area
           
           if (input$act_button == 0)
             return()
           
           if(input$pep1.1 == "all peptides") {
             
             plots <- list()
             
             for (j in 1:nlevels(prodata$Precursor)) {
               precursdata<-prodata[prodata$Precursor==levels(prodata$PrecursorTA)[j],] # subset for a particular precursor
               x=precursdata$Total.Area  # raw data for total area
               mu=mean(x[input$L:input$U]) # in-control process mean
               sd=sd(x[input$L:input$U]) # in-control process variance
               z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
               plots[[2*j-1]] <- IMR_plot1(z,j,input$L,input$U,"Total Area") +  ggtitle(bquote(atop(.(levels(prodata$PrecursorTA)[j]), "")))
               plots[[2*j]] <- IMR_plot2(z,j,input$L,input$U,"Total Area") +  ggtitle(bquote(atop(.(levels(prodata$PrecursorTA)[j]), "")))
               
             }
             # number_of_plots <- 2*nlevels(prodata$Precursor)
             # layout <- matrix(1:number_of_plots, ncol=2, byrow = TRUE)
             # isolate(multiplot(plotlist = plots, cols = 2, layout = layout))
             do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
               layout(autosize = F, width = 1500, height = my_height())
             
           }
           else{
           #j = as.numeric(input$pep1.1)
           j = which(levels(prodata$Precursor) == input$pep1.1)
           #prodata$Total.Area <- as.numeric(gsub(",","",prodata$Total.Area))
           precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],] # subset for a particular precursor
           x=precursdata$Total.Area; # raw data for total area
           #Main.title="Total Area" # main title of each plot
           mu=mean(x[input$L:input$U]); # in-control process mean
           sd=sd(x[input$L:input$U]); # in-control process variance
           z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
           
           plot1 <- IMR_plot1(z,j,input$L,input$U, "Total Area")
           plot2 <- IMR_plot2(z,j,input$L,input$U, "Total Area")
           subplot(plot1,plot2)
           
           }
         }
    #, height = my_height
    )
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