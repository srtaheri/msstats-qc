source("QCMetrics.R")
library(dplyr)
library(ggplot2)

CUSUM_plot <- function(prodata, z, j, L, U, Main.title, ytitle, type) {
  h <- 5 
  precursor.level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
  plot.data <- CUSUM.data.prepare(prodata, precursor.level, z, L, U, type)
  
  #ymax=ifelse(max(plot.data$CUSUM)>=h,(max(plot.data$CUSUM)),h)
  #ymin=ifelse(min(plot.data$CUSUM)<=-h,(min(plot.data$CUSUM)),-h)
  
  Main=Main.title
  
  x <- list(
    title =  paste("QCno - ", precursor.level),
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
    add_trace(x=c(0, max(plot.data$QCno)),y = c(h,h), marker=list(color="red" , size=4 , opacity=0.5), name = "UCL",showlegend = FALSE) %>%
    add_trace(x=c(0, max(plot.data$QCno)),y = c(-h,-h), marker=list(color="red" , size=4 , opacity=0.5), name = "LCL",showlegend = FALSE) %>%
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
    add_trace(x = plot.data[CUSUM.poz <= -h, ]$QCno,
              y = plot.data[CUSUM.poz <= -h, ]$CUSUM.poz,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.poz >= h, ]$QCno,
              y = plot.data[CUSUM.poz >= h, ]$CUSUM.poz,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    )%>%
    add_trace(x = plot.data[CUSUM.neg <= -h, ]$QCno,
              y = plot.data[CUSUM.neg <= -h, ]$CUSUM.neg,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    ) %>%
    add_trace(x = plot.data[CUSUM.neg >= h, ]$QCno,
              y = plot.data[CUSUM.neg >= h, ]$CUSUM.neg,
              mode = "markers",
              marker=list(color="red" , size=5 , opacity=0.5),
              showlegend = FALSE,name=""
    )
  
  return(p)
}

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
CP_plot <- function(prodata,z,j,Main.title,type, ytitle) {
  
  precursor_level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
  prodata_grouped_by_precursor <- prodata[prodata$Precursor==precursor_level,]
  ## Create variables 
  plot.data <- CP.data.prepare(prodata,z,j, type)
  y.max=max(plot.data$Et) # y axis upper limit
  y.min=0 # y axis lower limit
  
  x <- list(
    title = paste("QCno - ", levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j])
  )
  y <- list(
    title = ytitle
  )
  
  plot_ly(plot.data, x = QCno, y = Et
          ,type = "scatter"
          ,line = list(shape = "linear")
          ,showlegend = FALSE,name=""
          , text=prodata_grouped_by_precursor$Annotations
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
#########################################################################################################################
#########################################################################################################################
XmR_plot <- function(prodata,z,j,L,U,Main.title, type, ytitle) {
  
  precursor_level <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j]
  prodata_grouped_by_precursor <- prodata[prodata$Precursor==precursor_level,]
  plot.data <- XmR.data.prepare(prodata, z, j,L,U, type)
  print(plot.data)
  
  #y.max=ifelse(max(plot.data$t)>=UCL,(max(plot.data$t)),UCL)
  #y.min=ifelse(min(plot.data$t)<=LCL,(min(plot.data$t)),LCL)
  
  x <- list(
    title = paste("QCno - ", levels(reorder(prodata$Precursor,prodata$BestRetentionTime))[j])
  )
  y <- list(
    title = ytitle
  )
  plot_ly(plot.data, x = QCno, y = t, type = "scatter",
          name = "",  line = list(shape = "linear"),
          marker=list(color="dodgerblue" , size=4 , opacity=0.5)
          ,showlegend = FALSE
          , text=prodata_grouped_by_precursor$Annotations
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
  
  A1 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "Retention Time", type = 1)
  A2 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "Retention Time", type = 2)
  B1 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "Peak Assymetry", type = 1)
  B2 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "Peak Assymetry", type = 2)
  C1 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "FWHM", type = 1)
  C2 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "FWHM", type = 2)
  D1 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "Total Area", type = 1)
  D2 <- Compute.QCno.OutOfRange.XmR(prodata,L,U, metric = "Total Area", type = 2)
  my_data <- data.frame(y = c(A1,A2,B1,B2,C1,C2,D1,D2),
                        x = c(rep("Retention Time",length(A1)+length(A2)), 
                              rep("Peak Assymetry",length(B1)+length(B2)), 
                              rep("FWHM",length(C1)+length(C2)), 
                              rep("Total Area",length(D1)+length(D2))),
                        m = c(rep("individual value",length(A1)), rep("moving range",length(A2)),
                              rep("individual value",length(B1)), rep("moving range",length(B2)),
                              rep("individual value",length(C1)), rep("moving range",length(C2)),
                              rep("individual value",length(D1)), rep("moving range",length(D2))))
 
  pdat <- my_data %>%
    group_by(x, m) %>%
    do(data.frame(loc = density(.$y)$x,
                  dens = density(.$y)$y))
  pdat$dens <- pdat$dens * 7
  pdat$dens <- ifelse(pdat$m == 'individual value', pdat$dens * -1, pdat$dens)
  pdat$dens <- ifelse(pdat$x == 'Peak Assymetry', pdat$dens + 1, pdat$dens)
  pdat$dens <- ifelse(pdat$x == 'FWHM', pdat$dens + 2, pdat$dens)
  pdat$dens <- ifelse(pdat$x == 'Total Area', pdat$dens + 3, pdat$dens)
  
  
  ggplot(pdat, aes(dens, loc, fill = m, group = interaction(m, x))) +
     geom_polygon() +
    scale_x_continuous(breaks = 0:3, labels = c('Retention Time', 'Peak Assymetry','FWHM','Total Area')) +
    ylab('QCno') +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    ggtitle("Percentage of peptides with signal - XmR Chart")
}
#################################################################################################################
XmR.Radar.Plot <- function(prodata,L,U) {
  precursors <- levels(reorder(prodata$Precursor,prodata$BestRetentionTime))
  dat1 = data.frame(peptides = precursors,
                   OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Retention Time",type = 1),
                   group = rep("individual value",length(precursors)),
                   orderby = seq(1:length(precursors))
                   )
  dat2 = data.frame(peptides = precursors,
                    OutRangeQCno = Compute.QCno.OutOfRangePeptide.XmR(prodata,L,U,metric = "Retention Time",type = 2),
                    group = rep("moving range",length(precursors)),
                    orderby = seq(1:length(precursors))
                    )
  dat <- rbind(dat1,dat2)
  print(dat)
  ggplot(dat, aes(y = OutRangeQCno, x = reorder(peptides,orderby), group = group, colour = group)) +
    coord_polar() +
    geom_point() +
    geom_path() +
    labs(x = NULL)
  
}
#########################################################################################################################
#  CUSUM.summary.plot <- function(prodata, L, U) {
#   h <- 5
#   
# 
# 
# }
#################################################################################################################

#################################################################################################################
#################################################################################################################
#################################################################################################################

do.plot <- function(prodata, z, j, L, U, method, main.title, y.title, type) {
  if(method=="CUSUM") {
    CUSUM_plot(prodata, z, j, L, U, main.title, y.title, type)
  } else if(method=="CP") {
    CP_plot(prodata, z, j, main.title, type, y.title)
  } else if(method=="XmR") {
    XmR_plot(prodata, z, j, L, U,main.title, type, y.title)
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
  for (j in 1:nlevels(prodata$Precursor)) {
    z <- prepare_column(prodata, j, L, U, metric, normalization)
    multidata[1:length(z),j]<-z
  }
  colnames(multidata) <- levels(prodata$Precursor)
  multidata=data.frame(multidata)
  pairs(multidata, upper.panel = panel.cor, col = "blue")
}
#########################################################################################################################
metrics_box.plot <- function(prodata) {
  prodata$PrecursorRT <- reorder(prodata$Precursor,prodata$BestRetentionTime) # to plot boxplots y axis (Retention Time) in decreasing order
  RT <- plot_ly(prodata, y = BestRetentionTime, color = PrecursorRT, type = "box") %>% 
    layout(yaxis = list(title = "Retention Time"),showlegend = FALSE)
  
##########HEAD
  prodata$PrecursorPA <- reorder(prodata$Precursor,prodata$MaxEndTime - prodata$MinStartTime) # to plot boxplots in increasing order
  PA <- plot_ly(prodata, y = (MaxEndTime-MinStartTime), color = PrecursorPA, type = "box") %>%
  layout(yaxis = list(title = "Peak Assymetry"),showlegend = FALSE)

##########origin/master
  
  prodata$PrecursorTA <- reorder(prodata$Precursor,prodata$TotalArea) # to plot boxplots in decreasing order
  TPA <- plot_ly(prodata, y = TotalArea, color = PrecursorTA, type = "box") %>% 
    layout(yaxis = list(title = "Total Peak Area"),showlegend = FALSE)
  
  prodata$PrecursorFWHM <- reorder(prodata$Precursor,prodata$MaxFWHM) 
  FWHM <- plot_ly(prodata, y = MaxFWHM, color = PrecursorFWHM, type = "box") %>% 
    layout(yaxis = list(title = "FWHM"),showlegend = FALSE)
  
  return(subplot(RT, PA, TPA, FWHM, nrows = 4) %>%
           layout(autosize = F, width = 700, height = 1000))
}
################################################### vioplot2 #####################################
vioplot2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE, side="both") 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q2 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  radj <- ifelse(side == "right", 0, 1)
  ladj <- ifelse(side == "left", 0, 1)
  if (!(is.null(h))) 
    args <- c(args, h = h)
  med.dens <- rep(NA, n)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q2[i] <- quantile(data, 0.5)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    med.dat <- do.call("sm.density", 
                       c(list(data, xlim=est.xlim,
                              eval.points=med[i], display = "none")))
    med.dens[i] <- med.dat$estimate
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    med.dens[i] <- med.dens[i] * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(x = c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              y = c(base[[i]], rev(base[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - radj*boxwidth/2, 
             q1[i], 
             at[i] + ladj*boxwidth/2, 
             q3[i], col = rectCol)
        # median line segment
        lines(x = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              y = rep(med[i],3))
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), 
              c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - radj*boxwidth/2, q3[i], at[i] + 
               ladj*boxwidth/2, col = rectCol)
        lines(y = c(at[i] - radj*med.dens[i], 
                    at[i], 
                    at[i] + ladj*med.dens[i]),
              x = rep(med[i],3))
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}