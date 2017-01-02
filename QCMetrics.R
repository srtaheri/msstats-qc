COL.BEST.RET <- "Retention Time"
COL.FWHM <- "Full Width at Half Maximum"
COL.TOTAL.AREA <- "Total Peak Area"
COL.PEAK.ASS <- "Peak Assymetry"
#############################################################################################
#INPUTS : "prodata" is the data user uploads.
#         "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#         "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#         "metric" is one of these metrics: COL.BEST.RET,COL.FWHM, COL.TOTAL.AREA,COL.PEAK.ASS or a metric that user defines in his data set
#         "normalization" is either TRUE or FALSE
#Description of function : it gets the metric column only for the precursor chosen and either return the column as it is or normalize it and then return it
getMetricData <- function(prodata, precursorSelection, L, U, metric, normalization) {
  precursor.data<-prodata[prodata$Precursor==precursorSelection,] #"Precursor" is one of the columns in data that shows the name of peptides
  metricData <- 0

  if(is.null(metric)){
    return(NULL)  
  }
  
  metricData = precursor.data[,metric]
  if(normalization == TRUE) {
    mu=mean(metricData[L:U]) # in-control process mean
    sd=sd(metricData[L:U]) # in-control process variance
    if(sd == 0) {sd <- 0.0001}
    metricData=scale(metricData[1:length(metricData)],mu,sd) # transformation for N(0,1) )
    return(metricData)
  } else if(normalization == FALSE){
    return(metricData)
  }
  
}
#########################################################################################################
find_custom_metrics <- function(prodata) {
    prodata <- prodata[, which(colnames(prodata)=="Annotations"):ncol(prodata),drop = FALSE]
    nums <- sapply(prodata, is.numeric)
    other.metrics <- colnames(prodata[,nums])[1:ifelse(length(colnames(prodata[,nums]))<11,
                                                       length(colnames(prodata[,nums])),
                                                       10)
                                              ] # limiting custom metrics up to 10 metrics and not more
    if(any(is.na(other.metrics))) {
      return(c())
    }
    return(other.metrics)
}
################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads. 
#        "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column 
#        "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame for CUSUM that contains all the information needed to plot CUSUM
CUSUM.data.prepare <- function(prodata, metricData, precursorSelection, L, U, type) {

  k=0.5 
  #precursor.data only gets the data for the selected precursor
  precursor.data <- prodata[prodata$Precursor==precursorSelection,] #"Precursor" is one of the columns in data that shows the name of peptides
  
  v <- numeric(length(metricData))
  
  Cpoz <- numeric(length(metricData))
  Cneg <- numeric(length(metricData))
  
  for(i in 2:length(metricData)) {
    Cpoz[i]=max(0,(metricData[i]-(k)+Cpoz[i-1]))
    Cneg[i]=max(0,((-k)-metricData[i]+Cneg[i-1]))
  }
  
  if(type==2) {
    for(i in 2:length(metricData)) {
      v[i]=(sqrt(abs(metricData[i]))-0.822)/0.349
    }
    for(i in 2:length(metricData)) {
      Cpoz[i]=max(0,(v[i]-(k)+Cpoz[i-1]))
      Cneg[i]=max(0,((-k)-v[i]+Cneg[i-1]))
    }
  }

  QCno = 1:length(metricData)
  
  plot.data = 
    data.frame(QCno = QCno
               ,CUSUM.poz = Cpoz
               ,CUSUM.neg = -Cneg 
               ,Annotations=precursor.data$Annotations
               )
  
  return(plot.data)
}
###################################################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads. 
#        "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column 
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame for CP that contains all the information needed to plot Change Point
CP.data.prepare <- function(prodata, metricData, type) {
  
  Et <-  numeric(length(metricData)-1) # this is Ct in type 1, and Dt in type 2.
  SS<- numeric(length(metricData)-1)
  SST<- numeric(length(metricData)-1)
  tho.hat <- 0
  
  if(type == 1) {
    ## Change point analysis for mean (Single step change model)
    for(i in 1:(length(metricData)-1)) {
      Et[i]=(length(metricData)-i)*(((1/(length(metricData)-i))*sum(metricData[(i+1):length(metricData)]))-0)^2 #change point function
    }
    QCno=1:(length(metricData)-1) 
  } else if(type == 2) {
    ## Change point analysis for variance (Single step change model)  
    for(i in 1:length(metricData)) {
      SS[i]=metricData[i]^2
    }
    for(i in 1:length(metricData)) {
      SST[i]=sum(SS[i:length(metricData)])
      Et[i]=((SST[i]/2)-((length(metricData)-i+1)/2)*log(SST[i]/(length(metricData)-i+1))-(length(metricData)-i+1)/2) #change point function
    }
    QCno=1:length(metricData)
  }
  tho.hat = which(Et==max(Et)) # change point estimate
  plot.data <- data.frame(QCno,Et,tho.hat)

  return(plot.data)
}
###################################################################################################
#INPUT : "prodata" is the data user uploads.
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "data.metrics" is all the available metrics. It is defined in server.R
get_CP_tho.hat <- function(prodata, L, U, data.metrics) {
  tho.hat <- data.frame(tho.hat = c(), metric = c(), group = c(), y=c())
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  for(metric in data.metrics) {
    for (j in 1:nlevels(prodata$Precursor)) {
      metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = TRUE)
      mix <- rbind(
        data.frame(tho.hat = CP.data.prepare(prodata, metricData, type = 1)$tho.hat[1], metric = metric, group = "Individual Value", y=1.1),
        data.frame(tho.hat = CP.data.prepare(prodata, metricData, type = 2)$tho.hat[1], metric = metric, group = "Moving Range", y=-1.1)
      )
      tho.hat <- rbind(tho.hat, mix)
      
    }
  }
  
  return(tho.hat)
}
###################################################################################################
#z = metricData
#INPUT : "prodata" is the data user uploads.
#        "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column 
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame for XmR that contains all the information needed to plot XmR
XmR.data.prepare <- function(prodata, metricData, L,U, type) {
  t <- numeric(length(metricData)-1) 

  for(i in 2:length(metricData)) {
    t[i]=abs(metricData[i]-metricData[i-1]) # Compute moving range of metricData
  }
  #Main=Main.title
  
  QCno=1:length(metricData)
  
  if(type == 1) {
    UCL=mean(metricData[L:U])+2.66*sd(t[L:U])
    LCL=mean(metricData[L:U])-2.66*sd(t[L:U])
    t <- metricData
  } else if(type == 2) {
    ## Calculate MR chart statistics and limits
    
    UCL=3.267*sd(t[1:L-U])
    LCL=0
  }
  plot.data=data.frame(QCno,metricData,t,UCL,LCL)
  return(plot.data)
}
############################################################################################
#INPUTS : "prodata" is the data user uploads.
#         "metric" is one of these metrics: COL.BEST.RET,COL.FWHM, COL.TOTAL.AREA,COL.PEAK.ASS or a metric that user defines in his data set
#         "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : returns a data frame that is used in CUSUM.Summary.DataFrame function below to use for plotting summary plot for CUSUM in Summary Tab of shiny app
CUSUM.Summary.prepare <- function(prodata, metric, L, U,type) {
  h <- 5

  QCno <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = T)
    counter[1:length(metricData)] <- counter[1:length(metricData)]+1
    plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], L, U, type)

    sub.poz <- plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]
    sub.neg <- plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]

    y.poz[sub.poz$QCno] <- y.poz[sub.poz$QCno] + 1
    y.neg[sub.neg$QCno] <- y.neg[sub.neg$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno]
  pr.y.neg = y.neg[1:max_QCno]/counter[1:max_QCno]

  plot.data <- data.frame(QCno = rep(1:max_QCno,2),
                          pr.y = c(pr.y.poz, pr.y.neg),
                          group = ifelse(rep(type==1,2*max_QCno), 
                                         c(rep("Metric mean increase",max_QCno),
                                           rep("Metric mean decrease",max_QCno)),
                                         c(rep("Metric dispersion increase",max_QCno),
                                           rep("Metric dispersion decrease",max_QCno))),
                          metric = rep(metric,max_QCno*2)
     )
  return(plot.data)
}
############################################################################################
CUSUM.Summary.DataFrame <- function(prodata, data.metrics, L, U) {
  dat <- data.frame(QCno = c(),
                    pr.y = c(),
                    group = c(),
                    metric = c())
  for (metric in data.metrics) {
    data.1   <- CUSUM.Summary.prepare(prodata, metric = metric, L, U,type = 1)
    data.2   <- CUSUM.Summary.prepare(prodata, metric = metric, L, U,type = 2)
    data.2$pr.y <- -(data.2$pr.y)
    dat <- rbind(dat,data.1,data.2)
  }
  return(dat)
}
############################################################################################
XmR.Summary.prepare <- function(prodata, metric, L, U,type) {
  QCno    <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))
  
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    counter[1:length(metricData)] <- counter[1:length(metricData)]+1
    plot.data <- XmR.data.prepare(prodata, metricData , L , U , type)
    
    sub.poz <- plot.data[plot.data$t >= plot.data$UCL, ]
    sub.neg <- plot.data[plot.data$t <= plot.data$LCL, ]
    
    y.poz[sub.poz$QCno] <- y.poz[sub.poz$QCno] + 1
    y.neg[sub.neg$QCno] <- y.neg[sub.neg$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno]
  pr.y.neg = y.neg[1:max_QCno]/counter[1:max_QCno]
  
  plot.data <- data.frame(QCno = rep(1:max_QCno,2),
                          pr.y = c(pr.y.poz, pr.y.neg),
                          group = ifelse(rep(type==1,2*max_QCno), 
                                         c(rep("Metric mean increase",max_QCno),
                                           rep("Metric mean decrease",max_QCno)),
                                         c(rep("Metric dispersion increase",max_QCno),
                                           rep("Metric dispersion decrease",max_QCno))),
                          metric = rep(metric,max_QCno*2))
  return(plot.data)
}
###########################################################################################
XmR.Summary.DataFrame <- function(prodata, data.metrics, L, U) {
  dat <- data.frame(QCno = c(),
                    pr.y = c(),
                    group = c(),
                    metric = c())
  for (metric in data.metrics) {
    data.1   <- XmR.Summary.prepare(prodata, metric = metric, L, U,type = 1) 
    data.2   <- XmR.Summary.prepare(prodata, metric = metric, L, U,type = 2)
    data.2$pr.y <- -(data.2$pr.y)
    dat <- rbind(dat, data.1, data.2)
  }
  return(dat)
}
############################################################################################
XmR.heatmap.DataFrame <- function(prodata,precursorSelection, L, U, type) {

  metricDataRT <- getMetricData(prodata, precursorSelection, L = L, U = U, metric = COL.BEST.RET, normalization = F)
  RT <- XmR.data.prepare(prodata, metricDataRT, L, U, type)$t
  
  metricDataPA <- getMetricData(prodata, precursorSelection, L = L, U = U, metric = COL.PEAK.ASS, normalization = F)
  PA <- XmR.data.prepare(prodata, metricDataPA, L, U, type)$t
  
  metricDataFWHM <- getMetricData(prodata, precursorSelection, L = L, U = U, metric = COL.FWHM, normalization = F)
  FWHM <- XmR.data.prepare(prodata, metricDataFWHM, L, U, type)$t
  
  metricDataTPA <- getMetricData(prodata, precursorSelection, L = L, U = U, metric = COL.TOTAL.AREA, normalization = F)
  TPA <- XmR.data.prepare(prodata, metricDataTPA, L, U, type)$t
  
  datMean <- data.frame(RT = RT,
                    PA = PA,
                    FWHM = FWHM,
                    TPA = TPA)
  return(datMean)
}
############################################################################################
Compute.QCno.OutOfRangePeptide.XmR <- function(prodata,L,U,metric,type, XmR.type) {
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    plot.data <- XmR.data.prepare(prodata, metricData , L = L, U = U, type = type)
    if(XmR.type == "poz")
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t >= plot.data$UCL, ]$QCno))
    else
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t <= plot.data$LCL, ]$QCno))
  }
  return(QCno.out.range)
}
#############################################################################################
Compute.QCno.OutOfRangePeptide.CUSUM <- function(prodata,L,U,metric,type, CUSUM.type) {
  h <- 5
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = T)
    plot.data <- CUSUM.data.prepare(prodata, metricData, precursors[j], L, U, type)
    if(CUSUM.type == "poz")
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    else
    QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
  }
  return(QCno.out.range)
}
###############################################################################################################
XmR.Radar.Plot.prepare <- function(prodata,L,U, metric, type,group, XmR.type) {
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  precursors2 <- substring(precursors, first = 1, last = 3)
  QCno.length <- c()
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    QCno.length <- c(QCno.length,length(metricData))
  }
  dat <- data.frame(peptides = precursors2,
             OutRangeQCno  = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type, XmR.type),
             group         = rep(group,length(precursors)),
             orderby       = seq(1:length(precursors)),
             metric        = rep(metric, length(precursors)),
             tool          = rep("XmR",length(precursors)),
             probability   = (Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type, XmR.type)/QCno.length)
             )

  return(dat)
}
################################################################################################
XmR.Radar.Plot.DataFrame <- function(prodata, data.metrics, L,U) {
  dat <- data.frame(peptides = c(), OutRangeQCno = c(), group = c(),
                    orderby = c(), metric = c(), tool = c(),
                    probability   = c()
                    )
  for (metric in data.metrics) {
    data.1 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 1,group = "Metric mean increase", XmR.type = "poz")
    data.2 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 1,group = "Metric mean decrease", XmR.type = "neg")
    data.3 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 2,group = "Metric dispersion increase", XmR.type = "poz")
    data.4 <- XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 2,group = "Metric dispersion decrease", XmR.type = "neg")
    dat <- rbind(dat, data.1, data.2, data.3, data.4)
  }
  return(dat)
}
#################################################################################################################
CUSUM.Radar.Plot.prepare <- function(prodata,L,U, metric,type,group, CUSUM.type) {
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  precursors2 <- substring(precursors, first = 1, last = 3)
  QCno.length <- c()
  for(j in 1:length(precursors)) {
    metricData <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    QCno.length <- c(QCno.length,length(metricData))
  }
  dat <- data.frame(peptides = precursors2,
                    OutRangeQCno  = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = metric,type = type, CUSUM.type),
                    group         = rep(group,length(precursors)),
                    orderby       = seq(1:length(precursors)),
                    metric        = rep(metric, length(precursors)),
                    tool          = rep("XmR",length(precursors)),
                    probability   = (Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = metric,type = type, CUSUM.type)/QCno.length)
  )
  return(dat)
}
#################################################################################################
CUSUM.Radar.Plot.DataFrame <- function(prodata, data.metrics, L,U) {
  dat <- data.frame(peptides = c(), OutRangeQCno = c(), group = c(),
                    orderby = c(), metric = c(), tool = c(),
                    probability   = c()
  )
  for (metric in data.metrics) {
   data.1 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 1, group = "Metric mean increase", CUSUM.type = "poz")
   data.2 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 1, group = "Metric mean decrease", CUSUM.type = "neg")
   data.3 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 2, group = "Metric dispersion increase", CUSUM.type = "poz")
   data.4 <- CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 2, group = "Metric dispersion decrease", CUSUM.type = "neg")
   dat <- rbind(dat, data.1, data.2, data.3, data.4)
  }
  return(dat)
}
#######################################################################################################
