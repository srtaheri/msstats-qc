
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
    data.frame(QCno = QCno
               ,CUSUM.poz = Cpoz
               ,CUSUM.neg = -Cneg 
               ,Annotations=prodata_grouped_by_precursor$Annotations
               )
  
  return(plot.data)
}
###################################################################################################
CP.data.prepare <- function(prodata,j,z, type) {
  
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
    z <- prepare_column(prodata, j, L, U, metric = metric, normalization = T)
    counter[1:length(z)] <- counter[1:length(z)]+1
    plot.data <- CUSUM.data.prepare(prodata, precursors[j], z, L, U, type)

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
    z <- prepare_column(prodata, j = j, L = L, U = U, metric = metric, normalization = T)
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
    z <- prepare_column(prodata, j = j, L = L, U = U, metric = metric, normalization = T)
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
    z <- prepare_column(prodata, j, L, U, metric = metric, normalization = T)
    precursor.level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
    plot.data <- CUSUM.data.prepare(prodata, precursor.level, z, L, U, type)
    if(CUSUM.type == "poz")
      QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    else
    QCno.out.range <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
    # QCno.out.range1 <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.poz >= h | plot.data$CUSUM.poz <= -h, ]$QCno))
    # QCno.out.range2 <- c(QCno.out.range,length(plot.data[plot.data$CUSUM.neg >= h | plot.data$CUSUM.neg <= -h, ]$QCno))
    # QCno.out.range <- c(QCno.out.range1,QCno.out.range2)
    
  }
  #print(QCno.out.range)
  return(QCno.out.range)
}
#########################################################################################################
# XmR.Radar.Plot.prepare <- function(prodata,L,U) {
#   precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  
  # dat1.ret = data.frame(peptides = precursors,
  #                       OutRangeQCno  = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Retention Time",type = 1),
  #                       group         = rep("individual \n value",length(precursors)),
  #                       orderby       = seq(1:length(precursors)),
  #                       metric        = rep("Retention Time - XmR", length(precursors)),
  #                       tool          = rep("XmR",length(precursors))
  # )
  # 
  # dat2.ret = data.frame(peptides     = precursors,
  #                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Retention Time",type = 2),
  #                       group        = rep("moving \n range",length(precursors)),
  #                       orderby      = seq(1:length(precursors)),
  #                       metric       = rep("Retention Time - XmR", length(precursors)),
  #                       tool         = rep("XmR",length(precursors))
  # )

  # dat1.pa = data.frame(peptides     = precursors,
  #                      OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Peak Assymetry",type = 1),
  #                      group        = rep("individual \n value",length(precursors)),
  #                      orderby      = seq(1:length(precursors)),
  #                      metric       = rep("Peak Assymetry - XmR", length(precursors)),
  #                      tool         = rep("XmR",length(precursors))
  # )
  # dat2.pa = data.frame(peptides     = precursors,
  #                      OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Peak Assymetry",type = 2),
  #                      group        = rep("moving \n range",length(precursors)),
  #                      orderby      = seq(1:length(precursors)),
  #                      metric       = rep("Peak Assymetry - XmR", length(precursors)),
  #                      tool         = rep("XmR",length(precursors))
  # )
  # 
  # 
  # dat1.fwhm = data.frame(peptides     = precursors,
  #                        OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "FWHM",type = 1),
  #                        group        = rep("individual \n value",length(precursors)),
  #                        orderby      = seq(1:length(precursors)),
  #                        metric       = rep("FWHM - XmR", length(precursors)),
  #                        tool         = rep("XmR",length(precursors))
  # )
  # dat2.fwhm = data.frame(peptides     = precursors,
  #                        OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "FWHM",type = 2),
  #                        group        = rep("moving \n range",length(precursors)),
  #                        orderby      = seq(1:length(precursors)),
  #                        metric       = rep("FWHM - XmR", length(precursors)),
  #                        tool         = rep("XmR",length(precursors))
  # )
  # 
  # 
  # dat1.ta = data.frame(peptides     = precursors,
  #                      OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Total Area",type = 1),
  #                      group        = rep("individual \n value",length(precursors)),
  #                      orderby      = seq(1:length(precursors)),
  #                      metric       = rep("Total Area - XmR", length(precursors)),
  #                      tool         = rep("XmR",length(precursors))
  # )
  # dat2.ta = data.frame(peptides     = precursors,
  #                      OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Total Area",type = 2),
  #                      group        = rep("moving \n range",length(precursors)),
  #                      orderby      = seq(1:length(precursors)),
  #                      metric       = rep("Total Area - XmR", length(precursors)),
  #                      tool         = rep("XmR",length(precursors))
  # )
  # 
  # rbind(dat1.ret,dat2.ret,dat1.pa,dat2.pa,dat1.fwhm,dat2.fwhm,dat1.ta,dat2.ta)
  
# rbind(dat1.ret,dat2.ret)
#   
# }
##############################################
# XmR.Radar.Plot.prepare <- function(prodata,L,U, metric, type,group) {
#   precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
#   QCno.length <- c()
#   for(j in 1:length(precursors)) {
#     z <- prepare_column(prodata, j = j, L = L, U = U, metric = metric, normalization = T)
#     QCno.length <- c(QCno.length,length(z))
#   }
#   dat <- data.frame(peptides = precursors,
#              OutRangeQCno  = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type),
#              group         = rep(group,length(precursors)),
#              orderby       = seq(1:length(precursors)),
#              metric        = rep("Retention Time - XmR", length(precursors)),
#              tool          = rep("XmR",length(precursors)),
#              probability   = (Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = type)/QCno.length)
#              )
#   
#   return(dat)
# }
#######################################################################
XmR.Radar.Plot.prepare <- function(prodata,L,U,metric) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  m2 <- matrix(c(Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = 1),
                 Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = metric,type = 2)),
               nrow = 2, byrow = TRUE)
  
  group.names <- c("Individual Value","Moving Range")
  df2 <- data.frame(group = group.names,m2)
  colnames(df2)[2:7] <- precursors
  return(df2)
}
#################################################################################################################
CUSUM.Radar.Plot.prepare <- function(prodata,L,U) {
  
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))

  dat1.ret.poz = data.frame(peptides = precursors,
                        OutRangeQCno  = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric= "Retention Time",
                                                                             type = 1, CUSUM.type = "poz"),
                        group         = rep("Individual Value CUSUM+",length(precursors)),
                        orderby       = seq(1:length(precursors)),
                        metric        = rep("Retention Time - CUSUM", length(precursors)),
                        tool          = rep("CUSUM",length(precursors)) 
  )
  
  dat1.ret.neg = data.frame(peptides = precursors,
                            OutRangeQCno  = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric= "Retention Time",
                                                                                 type = 1, CUSUM.type = "neg"),
                            group         = rep("Individual Value CUSUM-",length(precursors)),
                            orderby       = seq(1:length(precursors)),
                            metric        = rep("Retention Time - CUSUM", length(precursors)),
                            tool          = rep("CUSUM",length(precursors)) 
  )
  
  dat2.ret.poz = data.frame(peptides     = precursors,
                        OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Retention Time",
                                                                            type = 2, CUSUM.type = "poz"),
                        group        = rep("Moving Range CUSUM+",length(precursors)),
                        orderby      = seq(1:length(precursors)),
                        metric       = rep("Retention Time - CUSUM", length(precursors)),
                        tool         = rep("CUSUM",length(precursors))
  )
  dat2.ret.neg = data.frame(peptides     = precursors,
                            OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Retention Time",
                                                                                type = 2, CUSUM.type = "neg"),
                            group        = rep("Moving Range CUSUM-",length(precursors)),
                            orderby      = seq(1:length(precursors)),
                            metric       = rep("Retention Time - CUSUM", length(precursors)),
                            tool         = rep("CUSUM",length(precursors))
  )
  
  dat1.pa.poz = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Peak Assymetry",
                                                                           type = 1, CUSUM.type = "poz"),
                       group        = rep("Individual Value CUSUM+",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Peak Assymetry - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
                       
  )
  dat1.pa.neg = data.frame(peptides     = precursors,
                           OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Peak Assymetry",
                                                                               type = 1, CUSUM.type = "neg"),
                           group        = rep("Individual Value CUSUM-",length(precursors)),
                           orderby      = seq(1:length(precursors)),
                           metric       = rep("Peak Assymetry - CUSUM", length(precursors)),
                           tool         = rep("CUSUM",length(precursors))
  )
  dat2.pa.poz = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Peak Assymetry",
                                                                           type = 2, CUSUM.type = "poz"),
                       group        = rep("Moving Range CUSUM+",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Peak Assymetry - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  dat2.pa.neg = data.frame(peptides     = precursors,
                           OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Peak Assymetry",
                                                                               type = 2, CUSUM.type = "neg"),
                           group        = rep("Moving Range CUSUM-",length(precursors)),
                           orderby      = seq(1:length(precursors)),
                           metric       = rep("Peak Assymetry - CUSUM", length(precursors)),
                           tool         = rep("CUSUM",length(precursors))
  )
  
  dat1.fwhm.poz = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "FWHM",
                                                                             type = 1, CUSUM.type = "poz"),
                         group        = rep("Individual Value CUSUM+",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - CUSUM", length(precursors)),
                         tool         = rep("CUSUM",length(precursors))
  )
  dat1.fwhm.neg = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "FWHM",
                                                                             type = 1, CUSUM.type = "neg"),
                         group        = rep("Individual Value CUSUM-",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - CUSUM", length(precursors)),
                         tool         = rep("CUSUM",length(precursors))
  )
  dat2.fwhm.poz = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "FWHM",
                                                                             type = 2, CUSUM.type = "poz"),
                         group        = rep("Moving Range CUSUM+",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - CUSUM", length(precursors)),
                         tool         = rep("CUSUM",length(precursors))
  )
  dat2.fwhm.neg = data.frame(peptides     = precursors,
                         OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "FWHM",
                                                                             type = 2, CUSUM.type = "neg"),
                         group        = rep("Moving Range CUSUM-",length(precursors)),
                         orderby      = seq(1:length(precursors)),
                         metric       = rep("FWHM - CUSUM", length(precursors)),
                         tool         = rep("CUSUM",length(precursors))
  )
  
  dat1.ta.poz = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Total Area",
                                                                           type = 1, CUSUM.type = "poz"),
                       group        = rep("Individual Value CUSUM+",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  dat1.ta.neg = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Total Area",
                                                                           type = 1, CUSUM.type = "neg"),
                       group        = rep("Individual Value CUSUM-",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  dat2.ta.poz = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Total Area",
                                                                           type = 2, CUSUM.type = "poz"),
                       group        = rep("Moving Range CUSUM+",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  dat2.ta.neg = data.frame(peptides     = precursors,
                       OutRangeQCno = Compute.QCno.OutOfRangePeptide.CUSUM(prodata,L,U,metric = "Total Area",
                                                                           type = 2, CUSUM.type = "neg"),
                       group        = rep("Moving Range CUSUM-",length(precursors)),
                       orderby      = seq(1:length(precursors)),
                       metric       = rep("Total Area - CUSUM", length(precursors)),
                       tool         = rep("CUSUM",length(precursors))
  )
  
  rbind(dat1.ret.poz,dat1.ret.neg,dat2.ret.poz,dat2.ret.neg,
        dat1.pa.poz,dat1.pa.neg,dat2.pa.poz,dat2.pa.neg,
        dat1.fwhm.poz,dat1.fwhm.neg,dat2.fwhm.poz,dat2.fwhm.neg,
        dat1.ta.poz,dat1.ta.neg,dat2.ta.poz,dat2.ta.neg)
  
}