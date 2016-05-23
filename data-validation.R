# here we put a selection of most column names that users use. The first element of each vector should be the best name that
# we suggest users to use and  which our code is based on. for example "Best.RT" and "Max FWHM" which are the first element
# of each vector in the list, are our suggestion so we wrote them in the fisrt place.
best_colnames <- list(
  c("Best.RT","best retention time", "retention time","rt","best ret time","intensity"),
  c("Max.FWHM","fwhm")
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
  data[data==""] <- NA
  colnames(data) <- unlist(lapply(colnames(data), function(x)guessColumnName(x)))
  data$Max.FWHM <- as.numeric(gsub(",","",data$Max.FWHM))
  data$Total.Area <- as.numeric(gsub(",","",data$Total.Area))
  #print(is.logical(colnames(data)[c(1,2)] == c("Precursor","Sample.File")) )
  #print(colnames(sample_data()))
  #print(colnames(data))
  return(data)
}