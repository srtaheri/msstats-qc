prepare_column <- function(prodata, j, L, U, metric, normalization) {
  prodata$Precursor <- reorder(prodata$Precursor,prodata$BestRetentionTime)
  precursdata<-prodata[prodata$Precursor==levels(prodata$Precursor)[j],]
  z <- 0
  if(metric == "Retention Time"){
    z = precursdata$BestRetentionTime # raw data for retention time
    #paste("bestsss")
  } else if(metric == "Peak Assymetry") {
    z = 2*precursdata$MinStartTime/(precursdata$MaxEndTime+precursdata$MinStartTime) # raw data for peak assymetry
  } else if(metric == "FWHM") {
    z = precursdata$MaxFWHM
  } else if(metric == "Total Area") {
    z = precursdata$TotalArea # raw data for total area
  } else {
    print("Error")
  }
  if(normalization == TRUE) {
    mu=mean(z[L:U]) # in-control process mean
    sd=sd(z[L:U]) # in-control process variance
    z=scale(z[1:length(z)],mu,sd) # transformation for N(0,1) )
    return(z)
  } else if(normalization == FALSE){
    return(z)
  }
  
}
#########################################################################################################
find_metrics <- function(prodata) {
  all_metrics_availabe <- c("BestRetentionTime",
                            "MaxFWHM",
                            "TotalArea","metric1","metric2","meric3","metric4","metric5")
 
  a <- all_metrics_availabe[which(all_metrics_availabe %in% colnames(prodata)==T)]
  
  return(a)
  
}
