#### CUSUM_chart plot1 function#############################################################################################
CUSUM_plot1 <- function(prodata, z,j,L,U,Main.title) {
  
  k=0.5
  h=5
  #prodata <- prodata()
  ##title <- title()
  
  Cpoz<- numeric(length(z))
  Cneg<- numeric(length(z))
  v <- numeric(length(z))
  
  for(i in 2:length(z))
  {
    Cpoz[i]=max(0,(z[i]-(k)+Cpoz[i-1]))
    Cneg[i]=max(0,((-k)-z[i]+Cneg[i-1]))
  }  
  QCno=rep(1:length(z),2)
  group=c(rep('CUSUM+',length(z)),rep('CUSUM-',length(z)))
  plot.data=data.frame(QCno,group,CUSUM=c(Cpoz,-Cneg),h)
  ymax=ifelse(max(plot.data$CUSUM)>=h,(max(plot.data$CUSUM)),h)
  ymin=ifelse(min(plot.data$CUSUM)<=-h,(min(plot.data$CUSUM)),-h)
  ##plot.subtitle=title[j]
  Main=Main.title
  
  x <- list(
    title = paste("QCno - ", levels(prodata$Precursor)[j])
    
  )
  y <- list(
    title = "CUSUMm"
  )
  
  p <- plot_ly(plot.data,
               x = plot.data[which(group == "CUSUM+"),"QCno"], 
               y = plot.data[which(group == "CUSUM+"),"CUSUM"], 
               line = list(color = "dodgerblue")
               , name = "CUSUM+"
               ,showlegend = FALSE
  ) %>%
    add_trace(x = plot.data[which(group == "CUSUM-"),"QCno"], 
              y = plot.data[which(group == "CUSUM-"),"CUSUM"], 
              line = list(color = "blue")
              , name = "CUSUM-"
              ,showlegend = FALSE
    ) %>%
    layout(xaxis = x,yaxis = y ) %>%
    
    add_trace(y = h, marker=list(color="red" , size=4 , opacity=0.5), name = "UCL",showlegend = FALSE) %>%
    add_trace(y = -h, marker=list(color="red" , size=4 , opacity=0.5), name = "LCL",showlegend = FALSE) %>%
    add_trace(x = plot.data[CUSUM < h & CUSUM > -h, ]$QCno,
              y = plot.data[CUSUM < h & CUSUM > -h, ]$CUSUM, 
              mode = "markers",
              marker=list(color="blue" , size=5 , opacity=0.5)
              #,name = levels(prodata$Precursor)[j]
              ,showlegend = FALSE
    ) %>%
    add_trace(x = plot.data[CUSUM <= -h, ]$QCno,
              y = plot.data[CUSUM <= -h, ]$CUSUM,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE
    ) %>%
    add_trace(x = plot.data[CUSUM >= h, ]$QCno,
              y = plot.data[CUSUM >= h, ]$CUSUM,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5)
              ,showlegend = FALSE
    ) 
  
  
  p
  
  ##rm(plot.data)
}

#### CUSUM_chart plot2 function ############################################################################################# 
CUSUM_plot2 <- function(prodata, z,j,L,U,Main.title){
  
  k=0.5
  h=5
  #prodata <- prodata()
  #title <- title()
  
  Cpoz<- numeric(length(z))
  Cneg<- numeric(length(z))
  v <- numeric(length(z))
  #plot.subtitle=title[j]
  
  for(i in 2:length(z))
  {
    Cpoz[i]=max(0,(z[i]-(k)+Cpoz[i-1]))
    Cneg[i]=max(0,((-k)-z[i]+Cneg[i-1]))
  }  
  QCno=rep(1:length(z),2)
  group=c(rep('CUSUM+',length(z)),rep('CUSUM-',length(z)))
  plot.data=data.frame(QCno,group,CUSUM=c(Cpoz,-Cneg),h)
  ymax=ifelse(max(plot.data$CUSUM)>=h,(max(plot.data$CUSUM)),h)
  ymin=ifelse(min(plot.data$CUSUM)<=-h,(min(plot.data$CUSUM)),-h)
  #plot.subtitle=title[j]
  Main=Main.title
  
  for(i in 2:length(z))
  {
    v[i]=(sqrt(abs(z[i]))-0.822)/0.349
  }
  vmean=mean(v[L:U]) 
  for(i in 2:length(z))
  {
    Cpoz[i]=max(0,(v[i]-(k)+Cpoz[i-1]));
    Cneg[i]=max(0,((-k)-v[i]+Cneg[i-1]))
  }
  plot.data=data.frame(QCno,group,CUSUM=c(Cpoz,-Cneg),h)
  ymax=ifelse(max(plot.data$CUSUM)>=h,(max(plot.data$CUSUM)),h)
  ymin=ifelse(min(plot.data$CUSUM)<=-h,(min(plot.data$CUSUM)),-h)
  
  x <- list(
    title =  paste("QCno - ", levels(prodata$Precursor)[j])
  )
  y <- list(
    title = "CUSUMv"
  )
  p <- plot_ly(plot.data,
               x = plot.data[which(group == "CUSUM+"),"QCno"], 
               y = plot.data[which(group == "CUSUM+"),"CUSUM"], 
               line = list(color = "dodgerblue")
               , name = "CUSUM+"
               ,showlegend = FALSE
  ) %>%
    add_trace(x = plot.data[which(group == "CUSUM-"),"QCno"], 
              y = plot.data[which(group == "CUSUM-"),"CUSUM"], 
              line = list(color = "blue")
              , name = "CUSUM-"
              ,showlegend = FALSE
    ) %>%
    layout(xaxis = x,yaxis = y) %>%
    add_trace(y = h, marker=list(color="red" , size=4 , opacity=0.5), name = "UCL",showlegend = FALSE) %>%
    add_trace(y = -h, marker=list(color="red" , size=4 , opacity=0.5), name = "LCL",showlegend = FALSE) %>%
    add_trace(x = plot.data[CUSUM < h & CUSUM > -h, ]$QCno,
              y = plot.data[CUSUM < h & CUSUM > -h, ]$CUSUM, 
              mode = "markers",
              marker=list(color="blue" , size=5 , opacity=0.5)
              #,name = levels(prodata$Precursor)[j]
              ,showlegend = FALSE
    ) %>%
    add_trace(x = plot.data[CUSUM <= -h, ]$QCno,
              y = plot.data[CUSUM <= -h, ]$CUSUM,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE
    ) %>%
    add_trace(x = plot.data[CUSUM >= h, ]$QCno,
              y = plot.data[CUSUM >= h, ]$CUSUM,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE
    ) 
  p
  
  
  #rm(plot.data)
  #rm(QCno)
}