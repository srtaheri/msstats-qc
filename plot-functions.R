source("QCMetrics.R")
source("http://pcwww.liv.ac.uk/~william/Geodemographic%20Classifiability/func%20CreateRadialPlot.r")
source('ggradar.R')
source('helper-functions.R')
library(dplyr)
library(ggplot2)
library(scales)

CUSUM.outrange.thld <- 5
#################################################################################################################
render.QC.chart <- function(prodata, precursorSelection, L, U, normalize.metric, plot.method, normalization.type, y.title1, y.title2){
  validate(
    need(!is.null(prodata), "Please upload your data")
  )
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  plots <- list()
  
  if(precursorSelection == "all peptides") {
    results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
      metricData <- getMetricData(prodata, precursors[j], L, U, metric = normalize.metric, normalization = normalization.type)
      plots[[2*j-1]] <<- do.plot(prodata, metricData, precursors[j],L,U, method=plot.method, y.title1, 1)
      plots[[2*j]] <<- do.plot(prodata, metricData, precursors[j],L,U, method=plot.method, y.title2, 2)
    })
    
    do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
      layout(autosize = F, width = 1400, height = nlevels(prodata$Precursor)*200)
  }
  
  else {
    metricData <- getMetricData(prodata, precursorSelection, L, U, metric = normalize.metric, normalization = normalization.type)
    
    plot1 <- do.plot(prodata, metricData, precursorSelection,L,U, method=plot.method,  y.title1, 1)
    plot2 <- do.plot(prodata, metricData, precursorSelection,L,U, method=plot.method,  y.title2, 2)
    
    subplot(plot1,plot2)
  }
}
#################################################################################################
CUSUM.plot <- function(prodata, metricData, precursor, L, U,  ytitle, type) {
  plot.data <- CUSUM.data.prepare(prodata, metricData, precursor, L, U, type)
  
  #ymax=ifelse(max(plot.data$CUSUM)>=CUSUM.outrange.thld,(max(plot.data$CUSUM)),CUSUM.outrange.thld)
  #ymin=ifelse(min(plot.data$CUSUM)<=-CUSUM.outrange.thld,(min(plot.data$CUSUM)),-CUSUM.outrange.thld)
  x <- list(
    title =  paste("QCno - ", precursor),
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
    add_trace(x=c(0, max(plot.data$QCno)),y = c(CUSUM.outrange.thld,CUSUM.outrange.thld), marker=list(color="red" , size=4 , opacity=0.5), name = "UCL",showlegend = FALSE) %>%
    add_trace(x=c(0, max(plot.data$QCno)),y = c(-CUSUM.outrange.thld,-CUSUM.outrange.thld), marker=list(color="red" , size=4 , opacity=0.5), name = "LCL",showlegend = FALSE) %>%
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
    add_trace(x = plot.data[CUSUM.poz <= -CUSUM.outrange.thld, ]$QCno,
              y = plot.data[CUSUM.poz <= -CUSUM.outrange.thld, ]$CUSUM.poz,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.poz >= CUSUM.outrange.thld, ]$QCno,
              y = plot.data[CUSUM.poz >= CUSUM.outrange.thld, ]$CUSUM.poz,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    )%>%
    add_trace(x = plot.data[CUSUM.neg <= -CUSUM.outrange.thld, ]$QCno,
              y = plot.data[CUSUM.neg <= -CUSUM.outrange.thld, ]$CUSUM.neg,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.neg >= CUSUM.outrange.thld, ]$QCno,
              y = plot.data[CUSUM.neg >= CUSUM.outrange.thld, ]$CUSUM.neg,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    )
  
  return(p)
}

#########################################################################################################################
CP.plot <- function(prodata, metricData, precursor, ytitle, type) {
  precursor.data <- prodata[prodata$Precursor==precursor,]
  ## Create variables 
  plot.data <- CP.data.prepare(prodata, metricData, type)
  y.max=max(plot.data$Et) # y axis upper limit
  y.min=0 # y axis lower limit
  
  x <- list(
    title = paste("QCno - ", precursor)
  )
  y <- list(
    title = ytitle
  )
  
  plot_ly(plot.data, x = QCno, y = Et
          ,type = "scatter"
          ,line = list(shape = "linear")
          ,showlegend = FALSE,name=""
          , text=precursor.data$Annotations
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
XmR.plot <- function(prodata, metricData, precursor, L, U, ytitle, type) {
  precursor.data <- prodata[prodata$Precursor==precursor,]
  plot.data <- XmR.data.prepare(prodata, metricData, L, U, type)
  #print(plot.data)
  
  #y.max=ifelse(max(plot.data$t)>=UCL,(max(plot.data$t)),UCL)
  #y.min=ifelse(min(plot.data$t)<=LCL,(min(plot.data$t)),LCL)
  
  x <- list(
    title = paste("QCno - ", precursor)
  )
  y <- list(
    title = ytitle
  )
  plot_ly(plot.data, x = QCno, y = t, type = "scatter",
          name = "",  line = list(shape = "linear"),
          marker=list(color="dodgerblue" , size=4 , opacity=0.5)
          ,showlegend = FALSE
          , text=precursor.data$Annotations
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
XmR.Summary.plot <- function(prodata,L,U) {
  data.rt.1   <- XmR.Summary.prepare(prodata, metric = "Retention Time", L, U,type = 1) %>% 
    mutate(group = "Individual value") %>% mutate(metric = "Retention Time")
  data.rt.2   <- XmR.Summary.prepare(prodata, metric = "Retention Time", L, U,type = 2)%>% 
    mutate(group = "Moving Range") %>% mutate(metric = "Retention Time")
  data.rt.2$pr.y <- -(data.rt.2$pr.y)
  
  
  data.pa.1   <- XmR.Summary.prepare(prodata, metric = "Peak Assymetry", L, U,type = 1)%>% 
    mutate(group = "Individual value") %>% mutate(metric = "Peak Assymetry")
  data.pa.2   <- XmR.Summary.prepare(prodata, metric = "Peak Assymetry", L, U,type = 2)%>% 
    mutate(group = "Moving Range") %>% mutate(metric = "Peak Assymetry")
  data.pa.2$pr.y <- -(data.pa.2$pr.y)
  
  
  data.fwhm.1 <- XmR.Summary.prepare(prodata, metric = "FWHM", L, U,type = 1)%>% 
    mutate(group = "Individual value") %>% mutate(metric = "FWHM")
  data.fwhm.2 <- XmR.Summary.prepare(prodata, metric = "FWHM", L, U,type = 2)%>% 
    mutate(group = "Moving Range") %>% mutate(metric = "FWHM")
   data.fwhm.2$pr.y <- -(data.fwhm.2$pr.y)
  
  
  
  data.ta.1   <- XmR.Summary.prepare(prodata, metric = "Total Area", L, U,type = 1)%>% 
    mutate(group = "Individual value") %>% mutate(metric = "Total Area")
  data.ta.2   <- XmR.Summary.prepare(prodata, metric = "Total Area", L, U,type = 2)%>% 
    mutate(group = "Moving Range") %>% mutate(metric = "Total Area")
   data.ta.2$pr.y <- -(data.ta.2$pr.y)
  
  
  dat <- rbind(data.rt.1,data.rt.2,data.pa.1,data.pa.2,data.fwhm.1,data.fwhm.2,data.ta.1,data.ta.2)
  theme_set(theme_gray(base_size = 10))
  gg <- ggplot(dat)
  gg <- gg + geom_hline(yintercept=0, alpha=0.5)
  gg <- gg + geom_point(aes(x=dat$QCno, y=dat$pr.y,colour = group, group = group))
  gg <- gg + geom_line(aes(x=dat$QCno, y=dat$pr.y, colour = group, group = group), size=0.3)
  gg <- gg + facet_wrap(~metric,nrow = 2)
  gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.1,1.1),breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
  gg <- gg + labs(x = "QC Numbers", y = "Percentage of peptides with signal")
  gg <- gg + ggtitle("XmR Chart")
  gg <- gg + theme(plot.title = element_text(size=20, face="bold",margin = margin(10, 0, 10, 0)),
                   axis.text.x=element_text(size=12, vjust=0.5),
                   axis.text.y=element_text(size=12, hjust=0.5),
                   axis.title.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   legend.text = element_text(size = 12),
                   legend.title=element_blank()
  )
  gg
  
}
###############################################################################################
CUSUM.Summary.plot <- function(prodata, L, U) {
  h <- 5
  data.rt.1   <- CUSUM.Summary.prepare(prodata, metric = "Retention Time", L, U,type = 1)
  data.rt.2   <- CUSUM.Summary.prepare(prodata, metric = "Retention Time", L, U,type = 2)
  data.rt.2$pr.y <- -(data.rt.2$pr.y)
  
  
  
  data.pa.1   <- CUSUM.Summary.prepare(prodata, metric = "Peak Assymetry", L, U,type = 1)
  data.pa.2   <- CUSUM.Summary.prepare(prodata, metric = "Peak Assymetry", L, U,type = 2)
  data.pa.2$pr.y <- -(data.pa.2$pr.y)
 
  
  
  data.fwhm.1 <- CUSUM.Summary.prepare(prodata, metric = "FWHM", L, U,type = 1)
  data.fwhm.2 <- CUSUM.Summary.prepare(prodata, metric = "FWHM", L, U,type = 2)
  data.fwhm.2$pr.y <- -(data.fwhm.2$pr.y)
  
  
  
  data.ta.1   <- CUSUM.Summary.prepare(prodata, metric = "Total Area", L, U,type = 1)
  data.ta.2   <- CUSUM.Summary.prepare(prodata, metric = "Total Area", L, U,type = 2)
  data.ta.2$pr.y <- -(data.ta.2$pr.y)
  
  
   dat <- rbind(data.rt.1,data.rt.2,data.pa.1,data.pa.2,data.fwhm.1,data.fwhm.2,data.ta.1,data.ta.2)
   #write.csv(file="dataHAHA.csv",dat)
   gg <- ggplot(dat)
   gg <- gg + geom_hline(yintercept=0, alpha=0.5)
   gg <- gg + geom_point(aes(x=dat$QCno, y=dat$pr.y,colour = group, group = group))
   gg <- gg + geom_line(aes(x=dat$QCno, y=dat$pr.y, colour = group, group = group), size=0.3)
   #gg <- gg + geom_smooth(aes(x=dat$QCno, y=dat$pr.y, colour = group, group = group))
   gg <- gg + facet_wrap(~metric,nrow = 2)
   gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.1,1.1),
                                 breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
   gg <- gg + ggtitle("CUSUM Chart")
   
   gg <- gg + labs(x = "QC Numbers", y = "Percentage of peptides with signal")
   gg <- gg + theme(plot.title = element_text(size=20, face="bold",margin = margin(10, 0, 10, 0)),
                    axis.text.x=element_text(size=12, vjust=0.5),
                    axis.text.y=element_text(size=12, hjust=0.5),
                    axis.title.y=element_text(size=15),
                    axis.title.x=element_text(size=15),
                    legend.text = element_text(size = 12),
                    legend.title=element_blank()
              )
   gg
  
}
 

####################################################################
 XmR.Radar.Plot <- function(prodata,L,U,metric) {
   precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))

   df <- rbind(XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 1, group = "Individual Value"),
               XmR.Radar.Plot.prepare(prodata,L,U,metric = metric, type = 2, group = "Moving Range"))

   coords <- by(df, df[,"group"], function(r){
     x <- getPolarCoord(r[,2])
     x <- cbind(x$x, x$y)
     x <- data.frame(rbind(r, r[1,]), x = x[,1], y = x[,2])
     return(x)
   })
   coords <- rbind(coords[[1]], coords[[2]])
   df <- data.frame(coords)

   # Plot
   smooth <- 1
   bgcolor <- "white"

   p <- plot_ly(data = df,
                x = x, y = y, mode = "lines",
                group = group,
                fill = "toself",
                line = list(smoothing = smooth, shape = "spline"),
                hoverinfo = "text") %>%

     add_trace(data = df,
               x = x, y = y, mode = "markers",
               marker = list(color = "white",
                             size = 10,
                             line = list(width = 2)),
               hoverinfo = "none",
               showlegend = F) %>%

   layout(xaxis = list(title = "", showgrid = F, zeroline = F, showticklabels = F,
                       domain = c(0.02, 0.48)),
          yaxis = list(title = "", showgrid = F, zeroline = F, showticklabels = F,
                       domain = c(0, 0.92)),
          font = list(family = "serif", size = 15),
          legend = list(x = 0.55, y = 0.9, bgcolor = "transparent"),
          plot_bgcolor = bgcolor,
          paper_bgcolor = bgcolor)
   # Add grids
   grid <- rbind(getPolarCoord(rep(5, 50), matrix = T, na = T),
                 getPolarCoord(rep(10, 80), matrix = T, na = T),
                 getPolarCoord(rep(15, 150), matrix = T, na = T),
                 getPolarCoord(rep(20, 170), matrix = T, na = T),
                 getPolarCoord(rep(25, 200), matrix = T, na = T),
                 getPolarCoord(rep(30, 220), matrix = T, na = T),
                 getPolarCoord(rep(35, 250), matrix = T, na = T))


   grid <- as.data.frame(grid)
   p <- add_trace(p, data = grid,
                  x = x, y = y, mode = "lines",
                  line = list(color = "#57788e", dash = "4px", width = 1),
                  showlegend = F,
                  hoverinfo = "none")

   inner <- getPolarCoord(rep(5.01, length(precursors)))
   outer <- getPolarCoord(rep(35.02, length(precursors)))

   x = t(cbind(inner$x, outer$x))
   y = t(cbind(inner$y, outer$y))

   x <- as.numeric(apply(x, 2, function(vec){
     return(c(vec, NA))
   }))

   y <- as.numeric(apply(y, 2, function(vec){
     return(c(vec, NA))
   }))

   linegrid <- data.frame(x = x, y = y)

   p <- add_trace(p, data = linegrid,
                  x = x, y = y, mode = "lines",
                  line = list(color = "#57788e", dash = "4px", width = 1),
                  showlegend = F,
                  hoverinfo = "none")

   # Add text
   labels <- precursors
   p <- add_trace(p, data = getPolarCoord(rep(35.04, 6)),
                  x = x, y = y, mode = "text", text = labels,
                  showlegend = F,
                  hoverinfo = "none",
                  textfont = list(family = "serif", color = "#808080"))
   
   # Add titles, description etc
   p <- layout(p,
               annotations = list(
                 list(xref = "paper", yref = "paper",
                      xanchor = "left", yanchor = "top",
                      x = 0.03, y = 1,
                      showarrow = F,
                      text = paste(metric,"-XmR"),
                      #text = "hello",
                      font = list(family = "serif",
                                  size = 25,
                                  color = "#4080bf")),

                 list(xref = "paper", yref = "paper",
                      xanchor = "left", yanchor = "top",
                      x = 0.03, y = 0.95,
                      showarrow = F,
                      text = '',
                      font = list(family = "serif",
                                  size = 16,
                                  color = "#679bcb")),

                 list(xref = "paper", yref = "paper",
                      xanchor = "left", yanchor = "top",
                      x = 0.60, y = 0.20,
                      showarrow = F,
                      align = "left",
                      text = "",
                      font = list(family = "arial",
                                  size = 12))
               ),

               shapes = list(
                 list(
                   xref = "paper", yref = "paper",
                   x0 = 0, x1 = 0.95,
                   y0 = 0, y1 = 1,
                   type = "rect",
                   layer = "above",
                   fillcolor = "rgba(191, 191, 191, 0.1)",
                   line = list(color = "transparent"))
               ))
p

 }
 ##########################################################################
 XmR.Radar.Plot.combine <- function(prodata,L,U) {
   subplot(
     XmR.Radar.Plot(prodata,L,U,metric = "Retention Time"),
     XmR.Radar.Plot(prodata,L,U,metric = "Peak Assymetry"),
     XmR.Radar.Plot(prodata,L,U,metric = "FWHM"),
     XmR.Radar.Plot(prodata,L,U,metric = "Total Area"),
     nrows = 2)
 }

#################################################################################################################
 CUSUM.Radar.Plot <- function(prodata,L,U, metric) {
   precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
   
   df <- rbind(
     CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 1, group = "Individual Value CUSUM+", CUSUM.type = "poz"),
     CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 1, group = "Individual Value CUSUM-", CUSUM.type = "neg"),
     CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 2, group = "Moving Range CUSUM+", CUSUM.type = "poz"),
     CUSUM.Radar.Plot.prepare(prodata,L,U, metric = metric, type = 2, group = "Moving Range CUSUM-", CUSUM.type = "neg")
             )
   write.csv(file="CUSUMRadar.csv",df)
   coords <- by(df, df[,"group"], function(r){
     x <- getPolarCoord(r[,2])
     x <- cbind(x$x, x$y)
     x <- data.frame(rbind(r, r[1,]), x = x[,1], y = x[,2])
     return(x)
   })
   coords <- rbind(coords[[1]], coords[[2]], coords[[3]],coords[[4]])
   df <- data.frame(coords)
   
   # Plot
   smooth <- 1
   bgcolor <- "white"
   
   p <- plot_ly(data = df,
                x = x, y = y, mode = "lines",
                group = group,
                fill = "toself",
                line = list(smoothing = smooth, shape = "spline"),
                hoverinfo = "text") %>%
     
     add_trace(data = df,
               x = x, y = y, mode = "markers",
               marker = list(color = "white",
                             size = 10,
                             line = list(width = 2)),
               hoverinfo = "none",
               showlegend = F) %>%
     
     layout(xaxis = list(title = "", showgrid = F, zeroline = F, showticklabels = F,
                         domain = c(0.02, 0.48)),
            yaxis = list(title = "", showgrid = F, zeroline = F, showticklabels = F,
                         domain = c(0, 0.92)),
            font = list(family = "serif", size = 15),
            legend = list(x = 0.55, y = 0.9, bgcolor = "transparent"),
            plot_bgcolor = bgcolor,
            paper_bgcolor = bgcolor)
   
   # Add grids
   grid <- rbind(getPolarCoord(rep(5, 50), matrix = T, na = T),
                 getPolarCoord(rep(10, 80), matrix = T, na = T),
                 getPolarCoord(rep(15, 150), matrix = T, na = T),
                 getPolarCoord(rep(20, 170), matrix = T, na = T),
                 getPolarCoord(rep(25, 200), matrix = T, na = T),
                 getPolarCoord(rep(30, 220), matrix = T, na = T),
                 getPolarCoord(rep(35, 250), matrix = T, na = T),
                 getPolarCoord(rep(40, 270), matrix = T, na = T),
                 getPolarCoord(rep(45, 290), matrix = T, na = T))
   
   
   grid <- as.data.frame(grid)
   p <- add_trace(p, data = grid,
                  x = x, y = y, mode = "lines",
                  line = list(color = "#57788e", dash = "4px", width = 1),
                  showlegend = F,
                  hoverinfo = "none")
   
   inner <- getPolarCoord(rep(5.01, length(precursors)))
   outer <- getPolarCoord(rep(45.02,length(precursors)))
   
   x = t(cbind(inner$x, outer$x))
   y = t(cbind(inner$y, outer$y))
   
   x <- as.numeric(apply(x, 2, function(vec){
     return(c(vec, NA))
   }))
   
   y <- as.numeric(apply(y, 2, function(vec){
     return(c(vec, NA))
   }))
   
   linegrid <- data.frame(x = x, y = y)
   
   p <- add_trace(p, data = linegrid,
                  x = x, y = y, mode = "lines",
                  line = list(color = "#57788e", dash = "4px", width = 1),
                  showlegend = F,
                  hoverinfo = "none")
   
   # Add text
   labels <- precursors
   p <- add_trace(p, data = getPolarCoord(rep(35.04, 6)),
                  x = x, y = y, mode = "text", text = labels,
                  showlegend = F,
                  hoverinfo = "none",
                  textfont = list(family = "serif", color = "#808080"))
   
   # Add titles, description etc
   p <- layout(p,
               annotations = list(
                 list(xref = "paper", yref = "paper",
                      xanchor = "left", yanchor = "top",
                      x = 0.03, y = 1,
                      showarrow = F,
                      text =  paste(metric,"-CUSUM") ,
                      #text = "hello",
                      font = list(family = "serif",
                                  size = 25,
                                  color = "#4080bf")),
                 
                 list(xref = "paper", yref = "paper",
                      xanchor = "left", yanchor = "top",
                      x = 0.03, y = 0.95,
                      showarrow = F,
                      text = '',
                      font = list(family = "serif",
                                  size = 16,
                                  color = "#679bcb")),
                 
                 list(xref = "paper", yref = "paper",
                      xanchor = "left", yanchor = "top",
                      x = 0.60, y = 0.20,
                      showarrow = F,
                      align = "left",
                      text = "",
                      font = list(family = "arial",
                                  size = 12))
               ),
               
               shapes = list(
                 list(
                   xref = "paper", yref = "paper",
                   x0 = 0, x1 = 0.95,
                   y0 = 0, y1 = 1,
                   type = "rect",
                   layer = "above",
                   fillcolor = "rgba(191, 191, 191, 0.1)",
                   line = list(color = "transparent"))
               ))
   p
 } 
################################################################################################################# 
 CUSUM.Radar.Plot.combine <- function(prodata,L,U) {
  subplot(
    CUSUM.Radar.Plot(prodata,L,U,metric = "Retention Time"),
    CUSUM.Radar.Plot(prodata,L,U,metric = "Peak Assymetry"),
    CUSUM.Radar.Plot(prodata,L,U,metric = "FWHM"),
    CUSUM.Radar.Plot(prodata,L,U,metric = "Total Area"),
    nrows = 2
  )
 }
#################################################################################################################
#################################################################################################################

do.plot <- function(prodata, z, precursor, L, U, method,  y.title, type) {
  if(method=="CUSUM") {
    CUSUM.plot(prodata, z, precursor, L, U,  y.title, type)
  } else if(method=="CP") {
    CP.plot(prodata, z, precursor, y.title, type)
  } else if(method=="XmR") {
    XmR.plot(prodata, z, precursor, L, U, y.title, type)
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
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  for (j in 1:nlevels(prodata$Precursor)) {
    z <- getMetricData(prodata, precursors[j], L, U, metric, normalization)
    multidata[1:length(z),j]<-z
  }
  colnames(multidata) <- precursors
  multidata=data.frame(multidata)
  pairs(multidata, upper.panel = panel.cor, col = "blue")
}
#########################################################################################################################
metrics_box.plot <- function(prodata) {
  prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$BestRetentionTime) # to plot boxplots y axis (Retention Time) in decreasing order
  RT <- plot_ly(prodata, y = BestRetentionTime, color = PrecursorRT, type = "box") %>% 
    layout(yaxis = list(title = "Retention Time"),showlegend = FALSE)
  
  prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$MaxEndTime - prodata$MinStartTime) # to plot boxplots in increasing order
  PA <- plot_ly(prodata, y = (MaxEndTime-MinStartTime), color = PrecursorPA, type = "box") %>%
  layout(yaxis = list(title = "Peak Assymetry"),showlegend = FALSE)
  
  prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$TotalArea) # to plot boxplots in decreasing order
  TPA <- plot_ly(prodata, y = TotalArea, color = PrecursorTA, type = "box") %>% 
    layout(yaxis = list(title = "Total Peak Area"),showlegend = FALSE)
  
  prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$MaxFWHM) 
  FWHM <- plot_ly(prodata, y = MaxFWHM, color = PrecursorFWHM, type = "box") %>% 
    layout(yaxis = list(title = "FWHM"),showlegend = FALSE)
  
  return(subplot(RT, PA, TPA, FWHM, nrows = 4) %>%
           layout(autosize = F, width = 700, height = 1000))
}