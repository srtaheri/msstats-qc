# here we put a selection of most column names that users use. The first element of each vector should be the best name that
# we suggest users to use and  which our code is based on. for example "BestRetentionTime" and "Max FWHM" which are the first element
# of each vector in the list, are our suggestion so we wrote them in the fisrt place.
best_colnames <- list(
  c("BestRetentionTime" ,"Best.RT","best retention time", "retention time","rt","best ret time","intensity"),
  c("MaxFWHM","fwhm","max.fwhm"),
  c("TotalArea","total area","TA","T.Area"),
  c("MinStartTime","min start time"),
  c("MaxEndTime", "max end time"),
  c("Precursor","PeptideSequence"),
  c("Annotations","anotations")
)
#### camelCaseSplit function ##############################################################################################
camelCaseSplit <- function(x) {
  # This function get a camelCase word and splits it.
  # Ex : camelCaseSplit("myComputerIsHere") ---> my Computer Is Here
  return(gsub("([a-z])([A-Z])", "\\1 \\L\\2", x, perl = TRUE))
}
#### punc_remove function #################################################################################################
punc_remove <- function(x){
  # This function removes any existing punctuation in your sentence or word and transfer it to space.
  # Ex1: punc_remove(Best.RT) --> Best RT     #Ex2: punc_remove(Best_RT) --> Best RT 
  return(gsub("[[:punct:]///' ]", " ", x))
}
#### clearString function ###############################################################################################
clearString <- function(x){
  # This function, gets a word or setence, Splits it (if it is a camelCase), removes any existing punctuations, and transfer
  # all Upper Case letters to lower case letters.
  # Ex: clearString("myName_isSara.Taheri") --> my name is sara taheri
  return(tolower(punc_remove(camelCaseSplit(x))))
}
#### guessColumnName function ###########################################################################################

# This function receives the data and check the column names of data and changes the column names if it is not the
# same names as our suggested sample data to fit our suggested sample data.guessColumnName <- function(x){
guessColumnName <- function(x){  
 
a <- clearString(x)
  
  max_index <- 0
  max <- -1
  for(i in 1:length(best_colnames)){
    col <- best_colnames[[i]]
    for(j in 1:length(col)){
      sim <- levenshteinSim(a,clearString(col[j]))
      if(sim > max){
        max <- sim
        max_index <- i
      }
    }
  }
  if (max > 0.6) {
    return(best_colnames[[max_index]][1])
  } 
  else {
    return(x)
  }
}

### Input_checking function #########################################################################################
input_checking <- function(data){

  ## save process output in each step #### creating a log file ########### from Meena's code
  allfiles <- list.files()
  
  num <- 0
  filenaming <- "msstatsqc"
  finalfile <- "msstatsqc.log"
  
  while(is.element(finalfile,allfiles)) {
    num <- num+1
    finalfile <- paste(paste(filenaming,num,sep="-"),".log",sep="")
  }
  
  session <- sessionInfo()
  sink("sessionInfo.txt")
  print(session)
  sink()
  
  processout <- as.matrix(read.table("sessionInfo.txt", header=T, sep="\t"))
  write.table(processout, file=finalfile, row.names=FALSE)
  
  processout <- rbind(processout, as.matrix(c(" "," ","MSstatsqc - dataProcess function"," "),ncol=1))
  ## preparing data
  data[data==""] <- NA
  colnames(data) <- unlist(lapply(colnames(data), function(x)guessColumnName(x)))
  
  data$MaxFWHM <- as.numeric(gsub(",","",data$MaxFWHM))
  data$TotalArea <- as.numeric(gsub(",","",data$TotalArea))
  data$BestRetentionTime <- as.numeric(gsub(",","",data$BestRetentionTime))
  data$MaxEndTime <- as.numeric(gsub(",","",data$MaxEndTime))
  data$MinStartTime <- as.numeric(gsub(",","",data$MinStartTime))
  # ## conditions
  required_column_names <- c("Precursor","BestRetentionTime","MaxFWHM","TotalArea","MinStartTime"
                             ,"MaxEndTime")
  provided_column_names <- colnames(data)
  if(all(required_column_names %in% provided_column_names)) {
    processout <- rbind(processout, c("The column names : provided - okay"))
    write.table(processout, file = finalfile, row.names = FALSE)
  } else if(!all(required_column_names %in% provided_column_names)) {
    missedInput <- which(!(requiredInputUpper %in% providedInputUpper))
    processout <- rbind(processout, c(paste("ERROR : The required input : ",
                                            paste(required_column_names[missedInput], collapse = ", "),
                                            " are not provided in input - stop")))
    #paste0("ERROR : The required input :", required_column_names[missedInput], ",are not provided in input - stop")
    write.table(processout, file = finalfile, row.names = FALSE)
    stop("Please check the required input. The required input needs (Precursor, BestRetentionTime, MaxFWHM, TotalArea, MinStartTime)")
    
  }
  
  return(data)
}
