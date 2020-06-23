#' Load a data frame
#'
#' Calculate triplet based on variant's chromosome sequence and position.
#' @param  dataframe TO output dataframe to be loaded
#' @return A triplet change for variant.
#' @export

triplet <- function(dataframe){
  for(i in 1:nrow(dataframe)){
    if (nchar(dataframe$REF[i]) ==1){
      dataframe$FR_REF[i] <- subseq(Hsapiens[[paste0("chr", dataframe$Sequence[i])]], dataframe$POS[i] - 1, dataframe$POS[i] +1) %>% as.character
    } else {
      dataframe$FR_REF[i] <- paste0(Hsapiens[[paste0("chr", dataframe$Sequence[i])]][dataframe$POS[i]-1], "M", Hsapiens[[paste0("chr", dataframe$Sequence[i])]][dataframe$POS[i]+nchar(dataframe$REF[i])]) %>% as.character
    }
  }
  dataframe$ALT <- ifelse(nchar(dataframe$ALT) == 1, dataframe$ALT, "M")
  dataframe$FR_ALT <- dataframe$FR_REF
  substr(dataframe$FR_ALT, 2, 2) <- dataframe$ALT
  return(dataframe)
}



#' Load a data frame
#'
#' Calculate triplet change frequence for each patient based on triplet change
#' @param  dataframe TO output dataframe after calculating triplet change.
#' @return Triplet frequency change for each patient.
#' @export

triplet_frequency <- function(dataframe){
  sample_set2 <- data.frame(REF = character())
  number_col <- ncol(dataframe)
  for (i in 1:length(unique(dataframe$SampleID))){
    sample_set <- dataframe[dataframe$SampleID == unique(dataframe$SampleID)[i],]
    table <- as.data.frame(table(sample_set$FR_ALT, sample_set$FR_REF))
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
