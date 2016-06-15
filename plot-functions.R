source("QCMetrics.R")

CUSUM_plot <- function(prodata, z, j, L, U, Main.title, ytitle, type) {
  h <- 5 
  precursor.level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
  plot.data <- CUSUM.data.prepare(prodata, precursor.level, z, L, U, type)
  #print(plot.data[1,1])
  #ymax=ifelse(max(plot.data$CUSUM)>=h,(max(plot.data$CUSUM)),h)
  #ymin=ifelse(min(plot.data$CUSUM)<=-h,(min(plot.data$CUSUM)),-h)
  
  Main=Main.title
  
  x <- list(
    title =  paste("QCno - ", precursor.level),
    range = c(0, max(plot.data$QCno))
  )
  y <- list(
    title = ytitle
  )
  
  p <- plot_ly(plot.data
               , x = QCno
               , y = CUSUM.poz
               , line = list(color = "dodgerblue")
               , name = "CUSUM+"
               , showlegend = FALSE
               , text = Annotations
  ) %>%
    add_trace(  x = QCno
              , y = CUSUM.neg
              , line = list(color = "blue")
              , name = "CUSUM-"
              , showlegend = FALSE
              , text = Annotations
    ) %>%
    layout(xaxis = x,yaxis = y, showlegend = FALSE) %>%
    add_trace(x=c(0, max(plot.data$QCno)),y = c(h,h), marker=list(color="red" , size=4 , opacity=0.5), name = "UCL",showlegend = FALSE) %>%
    add_trace(x=c(0, max(plot.data$QCno)),y = c(-h,-h), marker=list(color="red" , size=4 , opacity=0.5), name = "LCL",showlegend = FALSE) %>%
    add_trace(  x = QCno
              , y = CUSUM.neg
              , mode = "markers"
              , marker=list(color="blue" , size=5 , opacity=0.5)
              , showlegend = FALSE,name=""
    ) %>%
    add_trace(  x = QCno
                , y = CUSUM.poz
                , mode = "markers"
                , marker=list(color="blue" , size=5 , opacity=0.5)
                , showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.poz <= -h, ]$QCno,
              y = plot.data[CUSUM.poz <= -h, ]$CUSUM.poz,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.poz >= h, ]$QCno,
              y = plot.data[CUSUM.poz >= h, ]$CUSUM.poz,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    )%>%
    add_trace(x = plot.data[CUSUM.neg <= -h, ]$QCno,
              y = plot.data[CUSUM.neg <= -h, ]$CUSUM.neg,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.neg >= h, ]$QCno,
              y = plot.data[CUSUM.neg >= h, ]$CUSUM.neg,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    )
  
  return(p)
}
##########################################################################################################
CUSUM.summary.plot <- function(prodata, L, U,type, ytitle) {
  h <- 5
  plot.data.ret.time <- CUSUM.Summary.prepare(prodata, metric = "Retention Time", L, U,type)
  plot.data.peak.assymetry <- CUSUM.Summary.prepare(prodata, metric = "Peak Assymetry", L, U,type)
  plot.data.fwhm <- CUSUM.Summary.prepare(prodata, metric = "FWHM", L, U,type)
  plot.data.total.area <- CUSUM.Summary.prepare(prodata, metric = "Total Area", L, U,type)

  p <- plot_ly( 
                 x = plot.data.ret.time$QCno
               , y = plot.data.ret.time$pr.y.poz
               #,  mode = "markers"
               #, marker=list(color="dodgerblue" , size=8 , opacity=0.5)
               , name = "CUSUM+RT"
               , line = list(shape = "linear", color="dodgerblue")
               , showlegend = FALSE
  ) %>%
    add_trace(  
                  x = plot.data.ret.time$QCno
                , y = plot.data.ret.time$pr.y.neg
                #,  mode = "markers"
                #, marker=list(color="blue" , size=8 , opacity=0.5)
                , name = "CUSUM-RT"
                , line = list(shape = "linear", color="blue")
                , showlegend = FALSE
  ) %>%
    add_trace(  
               x = plot.data.peak.assymetry$QCno
              , y = plot.data.peak.assymetry$pr.y.poz
              #,  mode = "markers"
              #, marker=list(color="rgb(128, 42, 42)" , size=8 , opacity=0.5)
              , name = "CUSUM+PA"
              , line = list(shape = "linear", color="rgb(128, 42, 42)")
              , showlegend = FALSE
              ) %>%
    add_trace(
        x = plot.data.peak.assymetry$QCno
      , y = plot.data.peak.assymetry$pr.y.neg
      #,  mode = "markers"
      #, marker=list(color="rgb(205, 92, 92)" , size=8 , opacity=0.5)
      , name = "CUSUM+PA"
      , line = list(shape = "linear", color="rgb(205, 92, 92)")
      , showlegend = FALSE
    ) %>%
    add_trace(
      x = plot.data.fwhm$QCno
      , y = plot.data.fwhm$pr.y.poz
      #,  mode = "markers"
      #, marker=list(color="rgb(248, 117, 49)" , size=8 , opacity=0.5)
      , name = "CUSUM+FWHM"
      , line = list(shape = "linear", color="rgb(248, 117, 49)")
      , showlegend = FALSE
    ) %>%
  add_trace(
    x = plot.data.fwhm$QCno
    , y = plot.data.fwhm$pr.y.neg
    #,  mode = "markers"
    #, marker=list(color="rgb(94, 38, 5)" , size=8 , opacity=0.5)
    , name = "CUSUM+FWHM"
    , line = list(shape = "linear", color="rgb(94, 38, 5)")
    , showlegend = FALSE
  ) %>%
    add_trace(
      x = plot.data.total.area$QCno
      , y = plot.data.total.area$pr.y.poz
      #,  mode = "markers"
      #, marker=list(color="rgb(84, 99, 44)" , size=8 , opacity=0.5)
      , name = "CUSUM+TA"
      , line = list(shape = "linear", color="rgb(84, 99, 44)")
      , showlegend = FALSE
    ) %>%
    add_trace(
      x = plot.data.total.area$QCno
      , y = plot.data.total.area$pr.y.neg
      #,  mode = "markers"
      #, marker=list(color="rgb(156, 203, 25)" , size=8 , opacity=0.5)
      , name = "CUSUM-TA"
      , line = list(shape = "linear", color="rgb(156, 203, 25)")
      , showlegend = FALSE
    )
  
  return(p)
  
}
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
CP_plot <- function(prodata,z,j,Main.title,type, ytitle) {
  ## Create variables 
  precursor_level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
  prodata_grouped_by_precursor <- prodata[prodata$Precursor==precursor_level,]
  Et <-  numeric(length(z)-1) # this is Ct in type 1, and Dt in type 2.
  SS<- numeric(length(z)-1)
  SST<- numeric(length(z)-1)
  tho.hat <- 0
  
  Main = Main.title
  
  if(type == 1) {
    ## Change point analysis for mean (Single step change model)
    for(i in 1:length(z)-1) {
      Et[i]=(length(z)-i)*(((1/(length(z)-i))*sum(z[(i+1):length(z)]))-0)^2 #change point function
    }
    QCno=1:(length(z)-1) 
  } else if(type == 2) {
    ## Change point analysis for variance (Single step change model)  
    for(i in 1:length(z)) {
      SS[i]=z[i]^2
    }
    for(i in 1:length(z)) {
      SST[i]=sum(SS[i:length(z)])
      Et[i]=((SST[i]/2)-((length(z)-i+1)/2)*log(SST[i]/(length(z)-i+1))-(length(z)-i+1)/2) #change point function
    }
    QCno=1:length(z)
  }
  
  tho.hat = which(Et==max(Et)) # change point estimate
  plot.data=data.frame(QCno,Et,tho.hat) # dataframe for change point plot
  y.max=max(plot.data$Et) # y axis upper limit
  y.min=0 # y axis lower limit
  
  x <- list(
    title = paste("QCno - ", levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j])
  )
  y <- list(
    title = ytitle
  )
  
  plot_ly(plot.data, x = QCno, y = Et
          ,type = "scatter"
          ,line = list(shape = "linear")
          ,showlegend = FALSE,name=""
          , text=prodata_grouped_by_precursor$Annotations
  ) %>%
    layout(xaxis = x,yaxis = y) %>%
    add_trace( x = c(tho.hat,tho.hat), y = c(0, (max(Et)+2)) 
               ,marker=list(color="red", size=4, opacity=0.5)
               , mode = "lines"
               ,showlegend = FALSE,name=""
    ) %>%
    add_trace(x = QCno, y =  Et
              ,mode = "markers"
              , marker=list(color="blue" , size=8 , opacity=0.5)
              ,showlegend = FALSE,name=""
    )
}

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
IMR_plot <- function(prodata,z,j,L,U,Main.title, type, ytitle) {
  
  t <- numeric(length(z)-1) # z in plot 1, MR in plot 2
  precursor_level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
  prodata_grouped_by_precursor <- prodata[prodata$Precursor==precursor_level,]
  
  ## Calculate X chart statistics and limits 
  #UCL = 0
  #LCL = 0
  for(i in 2:length(z)) {
    t[i]=abs(z[i]-z[i-1]) # Compute moving range of z
  }
  Main=Main.title
  
  QCno=1:length(z)
  
  if(type == 1) {
    UCL=mean(z[L:U])+2.66*sd(t[L:U])
    LCL=mean(z[L:U])-2.66*sd(t[L:U])
    #UCL=3
    #LCL=-3
    t <- z
  } else if(type == 2) {
    ## Calculate MR chart statistics and limits

    UCL=3.267*sd(t[1:L-U])
    LCL=0
  }
  plot.data=data.frame(QCno,z,t,UCL,LCL)
  
  y.max=ifelse(max(plot.data$t)>=UCL,(max(plot.data$t)),UCL)
  y.min=ifelse(min(plot.data$t)<=LCL,(min(plot.data$t)),LCL)
  
  x <- list(
    title = paste("QCno - ", levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j])
  )
  y <- list(
    title = ytitle
  )
  plot_ly(plot.data, x = QCno, y = t, type = "scatter",
          name = "",  line = list(shape = "linear"),
          marker=list(color="dodgerblue" , size=4 , opacity=0.5)
          ,showlegend = FALSE
          , text=prodata_grouped_by_precursor$Annotations
  ) %>%
    layout(xaxis = x,yaxis = y) %>%
    add_trace(y = UCL, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines",showlegend = FALSE,name="UCL") %>%
    add_trace(y = LCL, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines",showlegend = FALSE,name="LCL") %>%
    add_trace(x = plot.data[t <= LCL, ]$QCno, y = plot.data[t <= LCL, ]$t
              , mode = "markers"
              , marker=list(color="red" , size=8 , opacity=0.5)
              ,showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[t >= UCL, ]$QCno, y = plot.data[t >= UCL, ]$t
              , mode = "markers"
              , marker=list(color="red" , size=8 , opacity=0.5)
              ,showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[t > LCL & t < UCL, ]$QCno, y = plot.data[t > LCL & t < UCL, ]$t
              , mode = "markers"
              , marker=list(color="blue" , size=8 , opacity=0.5)
              ,showlegend = FALSE,name=""
    )
}
#################################################################################################################
#################################################################################################################
#################################################################################################################

#################################################################################################################
#################################################################################################################
#################################################################################################################

do.plot <- function(prodata, z, j, L, U, method, main.title, y.title, type) {
  if(method=="CUSUM") {
    CUSUM_plot(prodata, z, j, L, U, main.title, y.title, type)
  } else if(method=="CP") {
    CP_plot(prodata, z, j, main.title, type, y.title)
  } else if(method=="ZMR") {
    IMR_plot(prodata, z, j, L, U,main.title, type, y.title)
  }
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
#########################################################################################################################
metrics_scatter.plot <- function(prodata, L, U, metric, normalization) {
  multidata<-matrix(0,length(prodata$Precursor),nlevels(prodata$Precursor))
  for (j in 1:nlevels(prodata$Precursor)) {
    z <- prepare_column(prodata, j, L, U, metric, normalization)
    multidata[1:length(z),j]<-z
  }
  colnames(multidata) <- levels(prodata$Precursor)
  multidata=data.frame(multidata)
  pairs(multidata, upper.panel = panel.cor, col = "blue")
}
#########################################################################################################################
metrics_box.plot <- function(prodata) {
  prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$BestRetentionTime) # to plot boxplots y axis (Retention Time) in decreasing order
  RT <- plot_ly(prodata, y = BestRetentionTime, color = PrecursorRT, type = "box") %>% 
    layout(yaxis = list(title = "Retention Time"),showlegend = FALSE)
  
##########HEAD
  prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$MaxEndTime - prodata$MinStartTime) # to plot boxplots in increasing order
  PA <- plot_ly(prodata, y = (MaxEndTime-MinStartTime), color = PrecursorPA, type = "box") %>%
  layout(yaxis = list(title = "Peak Assymetry"),showlegend = FALSE)

##########origin/master
  
  prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$TotalArea) # to plot boxplots in decreasing order
  TPA <- plot_ly(prodata, y = TotalArea, color = PrecursorTA, type = "box") %>% 
    layout(yaxis = list(title = "Total Peak Area"),showlegend = FALSE)
  
  prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$MaxFWHM) 
  FWHM <- plot_ly(prodata, y = MaxFWHM, color = PrecursorFWHM, type = "box") %>% 
    layout(yaxis = list(title = "FWHM"),showlegend = FALSE)
  
  return(subplot(RT, PA, TPA, FWHM, nrows = 4) %>%
           layout(autosize = F, width = 700, height = 1000))
}