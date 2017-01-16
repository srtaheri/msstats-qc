source("QCMetrics.R")
#source("http://pcwww.liv.ac.uk/~william/Geodemographic%20Classifiability/func%20CreateRadialPlot.r")
source('ggradar.R')
source('helper-functions.R')
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(tidyr)
library(viridis)
library(extrafont)

CUSUM.outrange.thld <- 5

#################################################################################################################
#INPUT : "prodata" is the data user uploads.
#        "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#        "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#        "metric" is one of these metrics: COL.BEST.RET,COL.FWHM, COL.TOTAL.AREA,COL.PEAK.ASS or a metric that user defines in his data set
#        "plot.method" is one of XmR, CUSUM or CP methods
#        "normalization" is either TRUE or FALSE
#        "y.title1" and "y.title2" are titles of left and right plot which are Individual Value and Moving Range
# DESCRIPTION : Draws together the "Individual Value" and "Moving Range plots" (left and right plot) for each metric and method.
render.QC.chart <- function(prodata, precursorSelection, L, U, metric, plot.method, normalization, y.title1, y.title2){
  validate(
    need(!is.null(prodata), "Please upload your data")
  )
  precursors <- levels(reorder(prodata$Precursor,prodata[,COL.BEST.RET]))
  plots <- list()

  if(precursorSelection == "all peptides") {
    results <- lapply(c(1:nlevels(prodata$Precursor)), function(j) {
      metricData <- getMetricData(prodata, precursors[j], L, U, metric = metric, normalization = normalization)
      plots[[2*j-1]] <<- do.plot(prodata, metricData, precursors[j],L,U, plot.method, y.title1, type = 1)
      plots[[2*j]] <<- do.plot(prodata, metricData, precursors[j],L,U, plot.method, y.title2, type = 2)
    })

    do.call(subplot,c(plots,nrows=nlevels(prodata$Precursor))) %>%
      layout(autosize = F, width = 1400, height = nlevels(prodata$Precursor)*200)
  }

  else {
    metricData <- getMetricData(prodata, precursorSelection, L, U, metric = metric, normalization)

    plot1 <- do.plot(prodata, metricData, precursorSelection,L,U, plot.method,  y.title1, type = 1)
    plot2 <- do.plot(prodata, metricData, precursorSelection,L,U, plot.method,  y.title2, type = 2)

    subplot(plot1,plot2)
  }
}
#################################################################################################################
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#          "plot.method" is one of XmR, CUSUM or CP methods
#          "y.title" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION : draw one plot (which is either Individual Value or Moving Range based on the type user chooses) for each metric and method
do.plot <- function(prodata, metricData, precursorSelection, L, U, plot.method,  y.title, type) {
  if(plot.method=="CUSUM") {
    CUSUM.plot(prodata, metricData, precursorSelection, L, U,  y.title, type)
  } else if(plot.method=="CP") {
    CP.plot(prodata, metricData, precursorSelection, y.title, type)
  } else if(plot.method=="XmR") {
    XmR.plot(prodata, metricData, precursorSelection, L, U, y.title, type)
  }
}
#################################################################################################
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#          "ytitle" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION: draws one CUSUM plot based on type for each given metric
CUSUM.plot <- function(prodata, metricData, precursorSelection, L, U,  ytitle, type) {
  plot.data <- CUSUM.data.prepare(prodata, metricData, precursorSelection, L, U, type)

  #ymax=ifelse(max(plot.data$CUSUM)>=CUSUM.outrange.thld,(max(plot.data$CUSUM)),CUSUM.outrange.thld)
  #ymin=ifelse(min(plot.data$CUSUM)<=-CUSUM.outrange.thld,(min(plot.data$CUSUM)),-CUSUM.outrange.thld)
  x <- list(
    title =  paste("QCno - ", precursorSelection),
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
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "ytitle" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION: draws one CP plot based on type for each given metric
CP.plot <- function(prodata, metricData, precursorSelection, ytitle, type) {
  precursor.data <- prodata[prodata$Precursor==precursorSelection,]
  ## Create variables
  plot.data <- CP.data.prepare(prodata, metricData, type)
  y.max=max(plot.data$Et) # y axis upper limit
  y.min=0 # y axis lower limit

  x <- list(
    title = paste("QCno - ", precursorSelection)
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
# INPUTS : "prodata" is the data user uploads.
#          "metricData" is the column of the data related to the metric we want. Forexample if we want retention time, it gives retention time column
#          "precursorSelection" is the precursor that user selects in Data Import tab. it can be either one precursor(peptide) or it can be "all peptides"
#          "L" and "U" are lower and upper bound of guide set that user choose in Data Import tab.
#          "ytitle" is the title of the plot which is either Individual Value or Moving Range
#          "type" is either 1 or 2. one is "Individual Value" plot and other "Moving Range" plot
#DESCRIPTION: draws one XmR plot based on type for each given metric
XmR.plot <- function(prodata, metricData, precursorSelection, L, U, ytitle, type) {
  precursor.data <- prodata[prodata$Precursor==precursorSelection,]
  plot.data <- XmR.data.prepare(prodata, metricData, L, U, type)

  #y.max=ifelse(max(plot.data$t)>=UCL,(max(plot.data$t)),UCL)
  #y.min=ifelse(min(plot.data$t)<=LCL,(min(plot.data$t)),LCL)

  x <- list(
    title = paste("QCno - ", precursorSelection)
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
  gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y, colour = "Change point"))
  gg <- gg + scale_color_manual(breaks = c("Metric mean increase",
                                           "Metric mean decrease",
                                           "Metric dispersion increase",
                                           "Metric dispersion decrease",
                                           "Change point"),
                                values = c("Metric mean increase" = "#E69F00",
                                           "Metric mean decrease" = "#56B4E9",
                                           "Metric dispersion increase" = "#009E73",
                                           "Metric dispersion decrease" = "#D55E00",
                                           "Change point" = "red"),
                                guide='legend')
  gg <- gg + guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,1,0),shape=c(NA,NA,NA,NA,16))))
  gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
  gg <- gg + annotate("text", x = 15, y = 1.3, label = "Mean")
  gg <- gg + annotate("text", x = 25, y = -1.3, label = "Dispersion")
  gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.4,1.4),breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
  gg <- gg + labs(x = "QC No", y = "% of out of control \nprecursors")
  gg <- gg + ggtitle("Overall Summary \nXmR")
  theme_set(theme_gray(base_size = 15)) # this will change the size of all the texts in all ggplot functions
  gg <- gg + theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
                   axis.text.x=element_text(size=12, vjust=0.5),
                   axis.text.y=element_text(size=12, hjust=0.5),
                   axis.title.y=element_text(size=12),
                   axis.title.x=element_text(size=12),
                   legend.text = element_text(size = 12),
                   legend.title=element_blank(),
                   plot.margin = unit(c(1,3,1,1), "lines")
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
   gg <- gg + geom_point(data = tho.hat.df, aes(x = tho.hat.df$tho.hat, y = tho.hat.df$y, colour = "Change point"))
   gg <- gg + scale_color_manual(breaks = c("Metric mean increase",
                                            "Metric mean decrease",
                                            "Metric dispersion increase",
                                            "Metric dispersion decrease",
                                            "Change point"),
                                 values = c("Metric mean increase" = "#E69F00",
                                            "Metric mean decrease" = "#56B4E9",
                                            "Metric dispersion increase" = "#009E73",
                                            "Metric dispersion decrease" = "#D55E00",
                                            "Change point" = "red"),
                                 guide='legend')
   gg <- gg + guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,1,0),shape=c(NA,NA,NA,NA,16))))
   gg <- gg + facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4))
   gg <- gg + annotate("text", x = 15, y = 1.3, label = "Mean")
   gg <- gg + annotate("text", x = 25, y = -1.3, label = "Dispersion")
   gg <- gg + scale_y_continuous(expand=c(0,0), limits = c(-1.4,1.4),
                                 breaks = c(1,0.5,0,-0.5,-1) ,labels = c(1,0.5,0,"0.5","1"))
   gg <- gg + ggtitle("Overall Summary \nCUSUM")

   gg <- gg + labs(x = "QC No", y = "% of out of control \nprecursors")
   gg <- gg + theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
                    axis.text.x=element_text(size=12, vjust=0.5),
                    axis.text.y=element_text(size=12, hjust=0.5),
                    axis.title.y=element_text(size=12),
                    axis.title.x=element_text(size=12),
                    legend.text = element_text(size = 12),
                    legend.title=element_blank(),
                    plot.margin = unit(c(1,3,1,1), "lines")
                    )

   gt <- ggplot_gtable(ggplot_build(gg))
   gt$layout$clip[gt$layout$name == "panel"] <- "off"
   grid.draw(gt)
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
    scale_fill_manual(breaks = c("Metric mean increase",
                                 "Metric mean decrease",
                                 "Metric dispersion increase",
                                 "Metric dispersion decrease"),
                      values = c("Metric mean increase" = "#E69F00",
                                 "Metric mean decrease" = "#56B4E9",
                                 "Metric dispersion increase" = "#009E73",
                                 "Metric dispersion decrease" = "#D55E00")) +
    scale_color_manual(breaks = c("Metric mean increase",
                                  "Metric mean decrease",
                                  "Metric dispersion increase",
                                  "Metric dispersion decrease"),
                       values = c("Metric mean increase" = "#E69F00",
                                  "Metric mean decrease" = "#56B4E9",
                                  "Metric dispersion increase" = "#009E73",
                                  "Metric dispersion decrease" = "#D55E00")) +
    facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
    geom_polygon(alpha=0.5)+
    ggtitle("Precursor Level Summary \nXmR") +
    xlab("") +
    ylab("# of out of control \nQC samples") +
    theme(
          axis.text.x = element_text(face="bold",size = rel(0.7)),
          axis.title.y=element_text(size=12),
          axis.text.y=element_text(size=12, hjust=0.5),
          plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
          legend.title=element_blank(),
          legend.text = element_text(size = 12),
          panel.grid.major = element_line(colour = "firebrick3",linetype = "dotted"),
          plot.margin = unit(c(1,3,1,1), "lines")
  )

}

#################################################################################################################
CUSUM.Radar.Plot <- function(prodata, data.metrics, L,U) {
  dat <- CUSUM.Radar.Plot.DataFrame(prodata, data.metrics, L,U)

  ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby),
                  group = group, colour = group, fill = group)) +
    coord_polar() +
    geom_point() +
    scale_fill_manual(breaks = c("Metric mean increase",
                                 "Metric mean decrease",
                                 "Metric dispersion increase",
                                 "Metric dispersion decrease"),
                      values = c("Metric mean increase" = "#E69F00",
                                 "Metric mean decrease" = "#56B4E9",
                                 "Metric dispersion increase" = "#009E73",
                                 "Metric dispersion decrease" = "#D55E00")) +
    scale_color_manual(breaks = c("Metric mean increase",
                                  "Metric mean decrease",
                                  "Metric dispersion increase",
                                  "Metric dispersion decrease"),
                       values = c("Metric mean increase" = "#E69F00",
                                  "Metric mean decrease" = "#56B4E9",
                                  "Metric dispersion increase" = "#009E73",
                                  "Metric dispersion decrease" = "#D55E00")) +
    facet_wrap(~metric,nrow = ceiling(length(data.metrics)/4)) +
    geom_polygon(alpha=0.5)+
    ggtitle("Precursor Level Summary \nCUSUM") +
    xlab("") +
    ylab("# of out of control \nQC samples") +

    theme(
      axis.text.x = element_text(face="bold",size = rel(0.7)),
      axis.title.y=element_text(size=12),
      axis.text.y=element_text(size=12, hjust=0.5),
      plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0)),
      legend.title=element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(colour = "firebrick3",linetype = "dotted"),
      plot.margin = unit(c(1,3,1,1), "lines")
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
#####################################################################################################
metrics_heat.map <- function(prodata,precursorSelection,data.metrics, method,peptideThresholdGood,peptideThresholdWarn, L, U, type) {

  data <- heatmap.DataFrame(prodata,precursorSelection, data.metrics,method,peptideThresholdGood,peptideThresholdWarn, L, U, type)

  p <- ggplot(data,aes(time,metric, group = bin, fill = bin))
  # p <- p + scale_fill_gradient2(low="#F0E442", high="#000000", mid="#D55E00",
  #                               midpoint=0.5,
  #                               #, limit=c(0.2,0.8)
  #                               name="Correlation\n(Pearson)",guide = "legend"
  #                               #, na.value = "red"
  #                               )

  # p <- p + scale_color_manual(values=c("Good" = "green","Bad" = "red","Warning" = "orange"),
  #                             breaks=c("Good","Bad","Warning"),
  #                             guide='legend')
  p <- p + geom_tile(colour="white",size=.1)
  p <- p + coord_equal()
  p <- p + theme_minimal(base_size = 10, base_family = "Trebuchet MS")
  p <- p + removeGrid()
  p <- p + rotateTextX()
  if(type == 1) {
     p <- p + ggtitle("XmR heat map - Mean",subtitle = "# Events per metric per date and time")
   }
   else {
     p <- p + ggtitle("XmR heat map - Dispersion",subtitle = "# Events per metric per date and time")
   }

  p <- p + labs(x=NULL, y=NULL)
  #p <- p + theme(plot.title=element_text(hjust=0))
  #p <- p + theme(axis.ticks=element_blank())
  p <- p +  theme(axis.text=element_text(size=12))
  #p <- p +  theme(legend.title=element_text(size=16))
  #p <- p +  theme(legend.text=element_text(size=12))


   p
}
