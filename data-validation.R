# here we put a selection of most column names that users use. The first element of each vector should be the best name that
# we suggest users to use and  which our code is based on. for example "Best.RT" and "Max FWHM" which are the first element
# of each vector in the list, are our suggestion so we wrote them in the fisrt place.
best_colnames <- list(
  c("Best.RT","best retention time", "retention time","rt","best ret time","intensity"),
  c("Max.FWHM","fwhm"),
  c("Total.Area","total area","TA","T.Area"),
  c("Min.Start.Time","min start time"),
  c("Max.End.Time", "max end time")
  #c("Precursor")
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
guessColumnName <- function(x){
  # best_colnames <- best_colnames()
  # This function receives the data and check the column names of data and changes the column names if it is not the
  # same names as our suggested sample data to fit our suggested sample data.
  a <- clearString(x)
  
  max_index <- 0
  max <- -1
  for(i in 1:length(best_colnames)){
    col <- best_colnames[[i]]
    for(j in 1:length(col)){
      sim <- levenshteinSim(a,col[j])
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
  #######################################
  data[data==""] <- NA
  colnames(data) <- unlist(lapply(colnames(data), function(x)guessColumnName(x)))
  
  data$Max.FWHM <- as.numeric(gsub(",","",data$Max.FWHM))
  data$Total.Area <- as.numeric(gsub(",","",data$Total.Area))
  data$Best.RT <- as.numeric(gsub(",","",data$Best.RT))
  data$Max.End.Time <- as.numeric(gsub(",","",data$Max.End.Time))
  data$Min.Start.Time <- as.numeric(gsub(",","",data$Min.Start.Time))
  
  required_column_names <- c("Precursor","Best.RT","Max.FWHM","Total.Area","Min.Start.Time"
                         ,"Max.End.Time")
  provided_column_names <- colnames(data)
  if(all(required_column_names %in% provided_column_names)) {
    processout <- rbind(processout, c("The column names : provided - okay"))
    write.table(processout, file = finalfile, row.names = FALSE)
  }
  return(data)
}