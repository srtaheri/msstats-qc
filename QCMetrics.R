
CUSUM.data.prepare <- function(prodata, precursor.level, z, L, U, type) {
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
    data.frame(QCno
               ,CUSUM.poz = Cpoz
               ,CUSUM.neg = -Cneg 
               ,Annotations=prodata_grouped_by_precursor$Annotations
               )
  
  return(plot.data)
}
###################################################################################################
CP.data.prepare <- function(prodata,z,j, type) {
  
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
XmR.data.prepare <- function(prodata, z, j,L,U, type) {
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
Compute.QCno.OutOfRange.XmR <- function(prodata,L,U, metric, type) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    z <- prepare_column(prodata, j, L, U, metric = metric, normalization = T)
    plot.data <- XmR.data.prepare(prodata, z, j,L,U, type)
    QCno.out.range <- c(QCno.out.range,plot.data[plot.data$t >= plot.data$UCL | plot.data$t <= plot.data$LCL, ]$QCno)
  }
  return(QCno.out.range)
}
############################################################################################
Compute.QCno.OutOfRangePeptide.XmR <- function(prodata,L,U,metric,type) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    z <- prepare_column(prodata, j, L, U, metric = metric, normalization = T)
    plot.data <- XmR.data.prepare(prodata, z, j,L,U, type)
    QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$t >= plot.data$UCL | plot.data$t <= plot.data$LCL, ]$QCno))
  }
  return(QCno.out.range)
}
#############################################################################################
Compute.QCno.OutOfRangePeptide.CUSUM <- function(prodata,L,U,metric,type) {
  h <- 5
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  QCno.out.range <- c()
  
  for(j in 1:length(precursors)) {
    z <- prepare_column(prodata, j, L, U, metric = metric, normalization = T)
    precursor.level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
    plot.data <- CUSUM.data.prepare(prodata, precursor.level, z, L, U, type)
    # if(CUSUM.type == "poz")    # I shoud add " , CUSUM.type " to the function
    # QCno.out.range <- c(QCno.out.range,plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno)
    # else
    # QCno.out.range <- c(QCno.out.range,plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno)
    QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno)) +
                      c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
    
  }
  print(QCno.out.range)
  return(QCno.out.range)
}
#########################################################################################################
XmR.Radar.Plot.prepare <- function(prodata,L,U) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  
  dat1.ret = data.frame(peptides = precursors,
                        OutRangeQCno  = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Retention Time",type = 1),
                        group         = rep("individual \n value",length(precursors)),
                        orderby       = seq(1:length(precursors)),
                        metric        = rep("Retention Time - XmR", length(precursors)),
                        tool          = rep("XmR",length(precursors)) 
  )
  dat2.ret = data.frame(peptides     = precursors,
                        OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Retention Time",type = 2),
                        group        = rep("moving \n range",length(precursors)),
                        orderby      = seq(1:length(precursors)),
                        metric       = rep("Retention Time - XmR", length(precursors)),
                        tool         = rep("XmR",length(precursors)) 
  )
  
  
  dat1.pa = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Peak Assymetry",type = 1),
                       group        = rep("individual \n value",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Peak Assymetry - XmR", length(precursors)),
                       tool         = rep("XmR",length(precursors)) 
  )
  dat2.pa = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Peak Assymetry",type = 2),
                       group        = rep("moving \n range",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Peak Assymetry - XmR", length(precursors)),
                       tool         = rep("XmR",length(precursors)) 
  )
  
  
  dat1.fwhm = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "FWHM",type = 1),
                         group        = rep("individual \n value",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - XmR", length(precursors)),
                         tool         = rep("XmR",length(precursors)) 
  )
  dat2.fwhm = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "FWHM",type = 2),
                         group        = rep("moving \n range",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - XmR", length(precursors)),
                         tool         = rep("XmR",length(precursors)) 
  )
  
  
  dat1.ta = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Total Area",type = 1),
                       group        = rep("individual \n value",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - XmR", length(precursors)),
                       tool         = rep("XmR",length(precursors)) 
  )
  dat2.ta = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Total Area",type = 2),
                       group        = rep("moving \n range",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - XmR", length(precursors)),
                       tool         = rep("XmR",length(precursors)) 
  )
  
  rbind(dat1.ret,dat2.ret,dat1.pa,dat2.pa,dat1.fwhm,dat2.fwhm,dat1.ta,dat2.ta)
  

  
}
#################################################################################################################
CUSUM.Radar.Plot.prepare <- function(prodata,L,U) {
  
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  
  dat1.ret = data.frame(peptides = precursors,
                        OutRangeQCno  = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric= "Retention Time",type = 1),
                        group         = rep("individual \n value",length(precursors)),
                        orderby       = seq(1:length(precursors)),
                        metric        = rep("Retention Time - CUSUM", length(precursors)),
                        tool          = rep("CUSUM",length(precursors)) 
  )
  dat2.ret = data.frame(peptides     = precursors,
                        OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Retention Time",type = 2),
                        group        = rep("moving \n range",length(precursors)),
                        orderby      = seq(1:length(precursors)),
                        metric       = rep("Retention Time - CUSUM", length(precursors)),
                        tool         = rep("CUSUM",length(precursors))
  )
  
  
  dat1.pa = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Peak Assymetry",type = 1),
                       group        = rep("individual \n value",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Peak Assymetry - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  dat2.pa = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Peak Assymetry",type = 2),
                       group        = rep("moving \n range",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Peak Assymetry - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  
  
  dat1.fwhm = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "FWHM",type = 1),
                         group        = rep("individual \n value",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - CUSUM", length(precursors)),
                         tool         = rep("CUSUM",length(precursors))
  )
  dat2.fwhm = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "FWHM",type = 2),
                         group        = rep("moving \n range",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - CUSUM", length(precursors)),
                         tool         = rep("CUSUM",length(precursors))
  )
  
  
  dat1.ta = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Total Area",type = 1),
                       group        = rep("individual \n value",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  dat2.ta = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Total Area",type = 2),
                       group        = rep("moving \n range",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  
  rbind(dat1.ret,dat2.ret,dat1.pa,dat2.pa,dat1.fwhm,dat2.fwhm,dat1.ta,dat2.ta)
  
}