normalize <- function(prodata, j, L, U, metric) {
  precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
  x <- 0
  if(metric == "Retention Time"){
    x = precursdata$Best.RT # raw data for retention time
  } else if(metric == "Peak Assymetry") {
    x = precursdata$Max.End.Time-precursdata$Min.Start.Time # raw data for peak assymetry
  } else if(metric == "FWHM") {
    x=precursdata$Max.FWHM
  } else if(metric == "Total Area") {
    x=precursdata$Total.Area # raw data for total area
  } else {
    print("Error")
  }
  
  mu=mean(x[L:U]) # in-control process mean
  sd=sd(x[L:U]) # in-control process variance
  z=scale(x[1:length(x)],mu,sd) # transformation for N(0,1) )
  return(z)
}