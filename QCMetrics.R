getMetricData <- function(prodata, precursor, L, U, metric, normalization) {
  precursor.data<-prodata[prodata$Precursor==precursor,]
  z <- 0
  if(metric == "Retention Time"){
    z = precursor.data$BestRetentionTime # raw data for retention time
    #paste("bestsss")
  } else if(metric == "Peak Assymetry") {
    z = 2*precursor.data$MinStartTime/(precursor.data$MaxEndTime+precursor.data$MinStartTime) # raw data for peak assymetry
  } else if(metric == "FWHM") {
    z = precursor.data$MaxFWHM
  } else if(metric == "Total Area") {
    z = precursor.data$TotalArea # raw data for total area
  } else {
    z = precursor.data[,metric]
  }
  if(normalization == TRUE) {
    mu=mean(z[L:U]) # in-control process mean
    sd=sd(z[L:U]) # in-control process variance
    z=scale(z[1:length(z)],mu,sd) # transformation for N(0,1) )
    return(z)
  } else if(normalization == FALSE){
    return(z)
  }
  
}
#########################################################################################################
find_metrics <- function(prodata) {
  prodata <- prodata[, (which(colnames(prodata)=="Annotations")+1):ncol(prodata)]
  
  nums <- sapply(prodata, is.numeric)
  
  other.metrics <- colnames(prodata[,nums])
  
  return(other.metrics)
  
}
################################################################
CUSUM.data.prepare <- function(prodata, z, precursor.level, L, U, type) {

  k=0.5 
  
  prodata_grouped_by_precursor <- prodata[prodata$Precursor==precursor.level,]
  
  v <- numeric(length(z))
  
  Cpoz <- numeric(length(z))
  Cneg <- numeric(length(z))
  
  for(i in 2:length(z)) {
    Cpoz[i]=max(0,(z[i]-(k)+Cpoz[i-1]))
    Cneg[i]=max(0,((-k)-z[i]+Cneg[i-1]))
  }
  
  if(type==2) {
    for(i in 2:length(z)) {
      v[i]=(sqrt(abs(z[i]))-0.822)/0.349
    }
    for(i in 2:length(z)) {
      Cpoz[i]=max(0,(v[i]-(k)+Cpoz[i-1]))
      Cneg[i]=max(0,((-k)-v[i]+Cneg[i-1]))
    }
  }

  QCno = 1:length(z)

  plot.data = 
    data.frame(QCno = QCno
               ,CUSUM.poz = Cpoz
               ,CUSUM.neg = -Cneg 
               ,Annotations=prodata_grouped_by_precursor$Annotations
               )
  
  return(plot.data)
}
###################################################################################################
CP.data.prepare <- function(prodata, z, type) {
  
  Et <-  numeric(length(z)-1) # this is Ct in type 1, and Dt in type 2.
  SS<- numeric(length(z)-1)
  SST<- numeric(length(z)-1)
  tho.hat <- 0
  
  #Main = Main.title
  
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

  return(data.frame(QCno,Et,tho.hat)) # dataframe for change point plot
}
###################################################################################################
XmR.data.prepare <- function(prodata, z, L,U, type) {
  t <- numeric(length(z)-1) # z in plot 1, MR in plot 2

  for(i in 2:length(z)) {
    t[i]=abs(z[i]-z[i-1]) # Compute moving range of z
  }
  #Main=Main.title
  
  QCno=1:length(z)
  
  if(type == 1) {
    UCL=mean(z[L:U])+2.66*sd(t[L:U])
    LCL=mean(z[L:U])-2.66*sd(t[L:U])
    t <- z
  } else if(type == 2) {
    ## Calculate MR chart statistics and limits
    
    UCL=3.267*sd(t[1:L-U])
    LCL=0
  }
  plot.data=data.frame(QCno,z,t,UCL,LCL)
  return(plot.data)
}
############################################################################################
CUSUM.Summary.prepare <- function(prodata, metric, L, U,type) {
  h <- 5

  QCno <- 1:nrow(prodata)
  y.poz <- rep(0,nrow(prodata))
  y.neg <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))

  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))

  for(j in 1:length(precursors)) {
    z <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = T)
    counter[1:length(z)] <- counter[1:length(z)]+1
    plot.data <- CUSUM.data.prepare(prodata, z, precursors[j], L, U, type)

    sub.poz <- plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]
    sub.neg <- plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]

    y.poz[sub.poz$QCno] <- y.poz[sub.poz$QCno] + 1
    y.neg[sub.neg$QCno] <- y.neg[sub.neg$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno]
  pr.y.neg = y.neg[1:max_QCno]/counter[1:max_QCno]
  # plot.data <- data.frame(QCno = QCno[1:max_QCno],
  #                         pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno],
  #                         pr.y.neg = y.poz[1:max_QCno]/counter[1:max_QCno]
  #  )
  plot.data <- data.frame(QCno = rep(1:max_QCno,2),
                          pr.y = c(pr.y.poz, pr.y.neg),
                          group = ifelse(rep(type==1,2*max_QCno), 
                                         c(rep("Individual Value CUSUM+",max_QCno),
                                           rep("Individual Value CUSUM-",max_QCno)),
                                         c(rep("Moving Range CUSUM+",max_QCno),
                                           rep("Moving Range CUSUM-",max_QCno))),
                          metric = rep(metric,max_QCno*2)
     )
  return(plot.data)
}
############################################################################################
XmR.Summary.prepare <- function(prodata, metric, L, U,type) {
  QCno    <- 1:nrow(prodata)
  y       <- rep(0,nrow(prodata))
  counter <- rep(0,nrow(prodata))
  
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  
  for(j in 1:length(precursors)) {
    z <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    counter[1:length(z)] <- counter[1:length(z)]+1
    plot.data <- XmR.data.prepare(prodata, z = z, L = L, U = U, type)
    
    sub <- plot.data[plot.data$t >= plot.data$UCL | plot.data$t <= plot.data$LCL, ]
    
    y[sub$QCno] <- y[sub$QCno] + 1
  }
  max_QCno <- max(which(counter!=0))
  
  plot.data <- data.frame(QCno = QCno[1:max_QCno],
                          pr.y = y[1:max_QCno]/counter[1:max_QCno]
  )
  return(plot.data)
}
############################################################################################
Compute.QCno.OutOfRangePeptide.XmR <- function(prodata,L,U,metric,type) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    z <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    plot.data <- XmR.data.prepare(prodata, z = z, L = L, U = U, type = type)
    QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t >= plot.data$UCL | plot.data$t <= plot.data$LCL, ]$QCno))
  }
  return(QCno.out.range)
}
#############################################################################################
Compute.QCno.OutOfRangePeptide.CUSUM <- function(prodata,L,U,metric,type, CUSUM.type) {
  h <- 5
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    z <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = T)
    plot.data <- CUSUM.data.prepare(prodata, z, precursors[j], L, U, type)
    if(CUSUM.type == "poz")
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    else
    QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
    # QCno.out.range1 <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    # QCno.out.range2 <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
    # QCno.out.range <- c(QCno.out.range1,QCno.out.range2)
    
  }
  
  return(QCno.out.range)
}

###############################################################################################################
XmR.Radar.Plot.prepare <- function(prodata,L,U, metric, type,group) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.length <- c()
  for(j in 1:length(precursors)) {
    z <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    QCno.length <- c(QCno.length,length(z))
  }
  dat <- data.frame(peptides = precursors,
             OutRangeQCno  = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type),
             group         = rep(group,length(precursors)),
             orderby       = seq(1:length(precursors)),
             metric        = rep(metric, length(precursors)),
             tool          = rep("XmR",length(precursors)),
             probability   = (Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type)/QCno.length)
             )

  return(dat)
}

#################################################################################################################

CUSUM.Radar.Plot.prepare <- function(prodata,L,U, metric,type,group, CUSUM.type) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.length <- c()
  for(j in 1:length(precursors)) {
    z <- getMetricData(prodata, precursors[j], L = L, U = U, metric = metric, normalization = T)
    QCno.length <- c(QCno.length,length(z))
  }
  dat <- data.frame(peptides = precursors,
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
CUSUM.Radar.Plot.DataFrame <- function(prodata,L,U) {
  dat <- rbind(CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Retention Time",type = 1,group = "Individual Value CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Retention Time",type = 1,group = "Individual Value CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Retention Time",type = 2,group = "Moving Range CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Retention Time",type = 2,group = "Moving Range CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Peak Assymetry",type = 1,group = "Individual Value CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Peak Assymetry",type = 1,group = "Individual Value CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Peak Assymetry",type = 2,group = "Moving Range CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Peak Assymetry",type = 2,group = "Moving Range CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "FWHM",type = 1,group = "Individual Value CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "FWHM",type = 1,group = "Individual Value CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "FWHM",type = 2,group = "Moving Range CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "FWHM",type = 2,group = "Moving Range CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Total Area",type = 1,group = "Individual Value CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Total Area",type = 1,group = "Individual Value CUSUM-", CUSUM.type = "neg"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Total Area",type = 2,group = "Moving Range CUSUM+", CUSUM.type = "poz"),
               CUSUM.Radar.Plot.prepare(prodata,L,U, metric = "Total Area",type = 2,group = "Moving Range CUSUM-", CUSUM.type = "neg")
          )
  return(dat)
}