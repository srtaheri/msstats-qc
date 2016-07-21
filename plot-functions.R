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
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  plots <- list()
  
  if(precursorSelection == "all peptides") {
    results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
      metricData <- getMetricData(prodata, precursors[j], L, U, metric = normalize.metric, normalization = normalization.type)
      plots[[2*j-1]] <<- do.plot(prodata, metricData, precursors[j],L,U, method=plot.method, y.title1, type = 1)
      plots[[2*j]] <<- do.plot(prodata, metricData, precursors[j],L,U, method=plot.method, y.title2, type = 2)
    })
    
    do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>% 
      layout(autosize = F, width = 1400, height = nlevels(prodata$Precursor)*200)
  }
  
  else {
    metricData <- getMetricData(prodata, precursorSelection, L, U, metric = normalize.metric, normalization = normalization.type)
    
    plot1 <- do.plot(prodata, metricData, precursorSelection,L,U, method=plot.method,  y.title1, type = 1)
    plot2 <- do.plot(prodata, metricData, precursorSelection,L,U, method=plot.method,  y.title2, type = 2)
    
    subplot(plot1,plot2)
  }
}
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
XmR.Summary.plot <- function(prodata,data.metrics, L, U) {
  
  dat <- XmR.Summary.DataFrame(prodata,data.metrics, L, U)
  tho.hat.df <- get_CP_tho.hat(prodata, L, U, data.metrics)

  gg <- ggplot(dat)
  gg <- gg + geom_hline(yintercept=0, alpha=0.5)
  gg <- gg + geom_smooth(method="loess",aes(x=dat$QCno, y=dat$pr.y,colour = group, group = group))
  gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y), color = "red")
  gg <- gg + scale_color_manual(breaks = c("Individual Value XmR+",
                                           "Individual Value XmR-",
                                           "Moving Range XmR+",
                                           "Moving Range XmR-"),
                                values = c("Individual Value XmR+" = "#E69F00",
                                           "Individual Value XmR-" = "#56B4E9",
                                           "Moving Range XmR+" = "#009E73",
                                           "Moving Range XmR-" = "#D55E00"))
  gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
  gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.2,1.2),breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
  gg <- gg + labs(x = "QC Numbers", y = "Percentage of peptides with signal")
  gg <- gg + ggtitle("XmR Chart")
  theme_set(theme_gray(base_size = 15)) # this will change the size of all the texts in all ggplot functions
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
CUSUM.Summary.plot <- function(prodata, data.metrics, L, U) {
   h <- 5
   dat <- CUSUM.Summary.DataFrame(prodata, data.metrics, L, U)
   tho.hat.df <- get_CP_tho.hat(prodata, L, U, data.metrics)
   
   gg <- ggplot(dat)
   gg <- gg + geom_hline(yintercept=0, alpha=0.5)
   gg <- gg + stat_smooth(method="loess", aes(x=dat$QCno, y=dat$pr.y, colour = group, group = group))
   gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y), color = "red")
   gg <- gg + scale_color_manual(breaks = c("Individual Value CUSUM+",
                                            "Individual Value CUSUM-",
                                            "Moving Range CUSUM+",
                                            "Moving Range CUSUM-"),
                                 values = c("Individual Value CUSUM+" = "#E69F00",
                                            "Individual Value CUSUM-" = "#56B4E9",
                                            "Moving Range CUSUM+" = "#009E73",
                                            "Moving Range CUSUM-" = "#D55E00"))
   gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
   gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.2,1.2),
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
XmR.Radar.Plot <- function(prodata, data.metrics, L,U) {

  dat <- XmR.Radar.Plot.DataFrame(prodata, data.metrics, L,U)
  #write.csv(file="dataRadar.csv",dat)

  ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                  group = group, colour = group, fill=group)) +
    coord_polar() +
    geom_point() +
    scale_fill_manual(breaks = c("Individual Value XmR+",
                                  "Individual Value XmR-",
                                  "Moving Range XmR+",
                                  "Moving Range XmR-"),
                       values = c("Individual Value XmR+" = "#E69F00",
                                  "Individual Value XmR-" = "#56B4E9",
                                  "Moving Range XmR+" = "#009E73",
                                  "Moving Range XmR-" = "#D55E00")) +
    scale_color_manual(breaks = c("Individual Value XmR+",
                                 "Individual Value XmR-",
                                 "Moving Range XmR+",
                                 "Moving Range XmR-"),
                      values = c("Individual Value XmR+" = "#E69F00",
                                 "Individual Value XmR-" = "#56B4E9",
                                 "Moving Range XmR+" = "#009E73",
                                 "Moving Range XmR-" = "#D55E00")) +
    facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
    geom_polygon(alpha=0.6)+
    ggtitle("Radar plot \n XmR Chart") +
    xlab("") +
    ylab("Number of Signals") +
    theme(
          axis.title.y=element_text(size=15),
          axis.text.y=element_text(size=12, hjust=0.5),
          plot.title = element_text(size=20, face="bold",margin = margin(10, 0, 10, 0)),
          legend.title=element_blank(),
          legend.text = element_text(size = 12)
  )
  #geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "curve"), data = df) 
}

#################################################################################################################
CUSUM.Radar.Plot <- function(prodata, data.metrics, L,U) {
  dat <- CUSUM.Radar.Plot.DataFrame(prodata, data.metrics, L,U)
  
  ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                  group = group, colour = group, fill = group)) +
    coord_polar() +
    geom_point() +
    scale_fill_manual(breaks = c("Individual Value CUSUM+",
                                 "Individual Value CUSUM-",
                                 "Moving Range CUSUM+",
                                 "Moving Range CUSUM-"),
                      values = c("Individual Value CUSUM+" = "#E69F00",
                                 "Individual Value CUSUM-" = "#56B4E9",
                                 "Moving Range CUSUM+" = "#009E73",
                                 "Moving Range CUSUM-" = "#D55E00")) +
    scale_color_manual(breaks = c("Individual Value XmR+",
                                  "Individual Value XmR-",
                                  "Moving Range XmR+",
                                  "Moving Range XmR-"),
                       values = c("Individual Value XmR+" = "#E69F00",
                                  "Individual Value XmR-" = "#56B4E9",
                                  "Moving Range XmR+" = "#009E73",
                                  "Moving Range XmR-" = "#D55E00")) +
    facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
    #geom_path(linejoin = "mitre", lineend = "butt") +
    geom_polygon(alpha=0.6)+
    ggtitle("Radar plot \n CUSUM Chart") +
    xlab("") +
    ylab("Number of Signals") +

    theme(
      axis.title.y=element_text(size=15),
      axis.text.y=element_text(size=12, hjust=0.5),
      plot.title = element_text(size=20, face="bold",margin = margin(10, 0, 10, 0)),
      legend.title=element_blank(),
      legend.text = element_text(size = 12)
    )
}
#################################################################################################################

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
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  for (j in 1:nlevels(prodata$Precursor)) {
    z <- getMetricData(prodata, precursors[j], L, U, metric, normalization)
    if(is.null(z))
      return(NULL)
    multidata[1:length(z),j]<-z
  }
  colnames(multidata) <- substring(precursors, first = 1, last = 3)
  multidata=data.frame(multidata)
  pairs(multidata, upper.panel = panel.cor, col = "blue")
}
#########################################################################################################################
metrics_box.plot <- function(prodata, data.metrics) {
  plots <- list()
  for(i in 1:length(data.metrics)) {
    metric <- data.metrics[i]
    precursor.data <- substring(reorder(prodata$Precursor,prodata[,metric]), first = 1, last = 3)
    plots[[i]] <- plot_ly(prodata, y = prodata[,metric], color = precursor.data, type = "box") %>% 
      layout(yaxis = list(title = metric),showlegend = FALSE)
  }

  p <- do.call(subplot,c(plots,nrows=length(plots))) %>% 
    layout(autosize = F, width = 700, height = 1000)
  return(p)
}