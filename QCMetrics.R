
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
  #print(plot.data)
  return(plot.data)
}
###################################################################################################
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
  plot.data <- data.frame(QCno[1:max_QCno], 
                          pr.y.poz = y.poz[1:max_QCno]/counter[1:max_QCno], 
                          pr.y.neg = y.neg[1:max_QCno]/counter[max_QCno])
  return(plot.data)
}