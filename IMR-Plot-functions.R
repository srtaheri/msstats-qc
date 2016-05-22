IMR_plot <- function(prodata,z,j,L,U,Main.title, type, ytitle) {

  t <- numeric(length(z)-1) # z in plot 1, MR in plot 2
  
  ## Calculate X chart statistics and limits 
  UCL = 0
  LCL = 0
  #UCLI=MRmean+2.66*sd(MR[L:U])
  #LCLI=MRmean-2.66*sd(MR[L:U])
  
  Main=Main.title
  
  QCno=1:length(z)
  
  if(type == 1) {
    UCL=3
    LCL=-3
    t <- z
  } else if(type == 2) {
    ## Calculate MR chart statistics and limits
    for(i in 2:length(z)) {
      t[i]=abs(z[i]-z[i-1]) # Compute moving range of z
    }
    UCL=3.267*sd(t[1:L-U])
    LCL=0
  }
  plot.data=data.frame(QCno,z,t,UCL,LCL)

  y.max=ifelse(max(plot.data$t)>=UCL,(max(plot.data$t)),UCL)
  y.min=ifelse(min(plot.data$t)<=LCL,(min(plot.data$t)),LCL)

  x <- list(
    title = paste("QCno - ", levels(prodata$Precursor)[j])
  )
  y <- list(
    title = ytitle
  )
  plot_ly(plot.data, x = QCno, y = t, type = "scatter",
          name = "linear",  line = list(shape = "linear"),
          marker=list(color="dodgerblue" , size=4 , opacity=0.5)
          ,showlegend = FALSE
          ) %>%
    layout(xaxis = x,yaxis = y) %>%
    add_trace( y = UCL, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines",showlegend = FALSE) %>%
    add_trace(y = LCL, marker=list(color="red" , size=4 , opacity=0.5), mode = "lines",showlegend = FALSE) %>%
    add_trace(x = plot.data[t <= LCL, ]$QCno, y = plot.data[t <= LCL, ]$t
              , mode = "markers"
              , marker=list(color="red" , size=8 , opacity=0.5)
              ,showlegend = FALSE
    ) %>%
    add_trace(x = plot.data[t >= UCL, ]$QCno, y = plot.data[t >= UCL, ]$t
              , mode = "markers"
              , marker=list(color="red" , size=8 , opacity=0.5)
              ,showlegend = FALSE
    ) %>%
    add_trace(x = plot.data[t > LCL & t < UCL, ]$QCno, y = plot.data[t > LCL & t < UCL, ]$t
              , mode = "markers"
              , marker=list(color="blue" , size=8 , opacity=0.5)
              ,showlegend = FALSE
    )
}

