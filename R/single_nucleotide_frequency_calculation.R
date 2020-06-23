#' Load a data frame
#'
#' Calculate single nucleotide change frequence for each patient based on single nucleotide change
#' @param  dataframe TO output dataframe.
#' @return Dataframe includes single nucleotide change frequency for each patient.
#' @export



single_nucleotide_change_frequency <- function(dataframe){
  sample_set2 <- data.frame(REF = character())
  number_col <- ncol(dataframe)
  dataframe$REF <- ifelse(nchar(dataframe$REF) == 1, dataframe$REF, "M")
  dataframe$ALT <- ifelse(nchar(dataframe$ALT) == 1, dataframe$ALT, "M")
  for (i in 1:length(unique(dataframe$SampleID))){
    sample_set <- dataframe[dataframe$SampleID == unique(dataframe$SampleID)[i],]
    table <- as.data.frame(table(sample_set$ALT, sample_set$REF))
    table <- table[table$Freq != 0,]
    for (x in 1:nrow(table)){
      sample_set[[paste0(table[x,2], "_", table[x,1])]] <- table[x,3]
    }
    sample_set2 <- base::merge.data.frame(x= sample_set, y = sample_set2, all = TRUE)
  }
  number_col2 <- ncol(sample_set2)
  sample_set3 <- sample_set2[,(number_col+1):number_col2]
  sample_set3[is.na(sample_set3)] <- 0
  sample_set4 <- cbind(sample_set2[,1:number_col], sample_set3)
  return(sample_set4)
}

