#' Load a data frame
#'
#' feature manipulation
#' @param  dataframe TO output dataframe to be loaded
#' @return A dataframe after feature manipulation
#' @export


feature_manipulation <- function(dataframe){
  dataframe <- dataframe %>%
    select(-c("ESP6500.Total", "GERP", "SIFT.Score", "MutationTaster.Score", "Variant.ID","POS",
              "Genomic.Variant", "NCBI.Gene.ID", "Transcript.ID", "Preferred.Transcript",
              "Transcript.Variant", "Protein.ID", "Exon.Number", "COSMIC.Mutation.ID",
              "COSMIC.Amino.Acid.Change", "COSMIC.Transcript.Change", "X1KG.Total", "CADD.Phred",
              "dbSNP.ID", "dbSNP.Build", "Dataset","Codon.Change", "MutationTaster.Pred"))
  dataframe$ExAC.Total <- replace_na(dataframe$ExAC.Total, replace = mean(dataframe$ExAC.Total, na.rm = TRUE))


  dataframe$Category <- ifelse(dataframe$Category == "FP", 1, 0) %>% as.numeric()
  dataframe$Sequence <- revalue(dataframe$Sequence, c(X = 23, Y =24)) %>% as.numeric()


  dataframe$REF <- ifelse(nchar(dataframe$REF) == 1, dataframe$REF, "M")
  dataframe$ALT <- ifelse(nchar(dataframe$ALT) == 1, dataframe$ALT, "M")
  dataframe$REF_ALT <- paste0(dataframe$REF, "_", dataframe$ALT)
  dataframe$FR_REF_ALT <- paste0(dataframe$FR_REF, "_", dataframe$FR_ALT)
  for(unique_value in unique(dataframe$REF_ALT)){
    dataframe[paste("REF_ALT", unique_value, sep = ".")] <- ifelse(dataframe$REF_ALT == unique_value, 1, 0) %>% as.numeric()
  }
  for(unique_value in unique(dataframe$FR_REF_ALT)){
    dataframe[paste("FR_REF_ALT", unique_value, sep = ".")] <- ifelse(dataframe$FR_REF_ALT == unique_value, 1, 0) %>% as.numeric()
  }

  for(unique_value in unique(dataframe$Variant.Type)){
    dataframe[paste("Variant.Type", unique_value, sep = ".")] <- ifelse(dataframe$Variant.Type == unique_value, 1, 0) %>% as.numeric()
  }

  for(unique_value in unique(dataframe$Transcript.Biotype)){
    dataframe[paste("Transcript.Biotype", unique_value, sep = ".")] <- ifelse(dataframe$Transcript.Biotype == unique_value, 1, 0) %>% as.numeric()
  }


  Var_freq <- dataframe%>% dplyr::count(Variant.Effect) %>% as.data.frame
  Var_freq <- Var_freq[order(Var_freq, decreasing = TRUE), ]
  oth <- Var_freq[Var_freq$n < 1000, ]$Variant.Effect
  dataframe$Variant.Effect <- ifelse(dataframe$Variant.Effect %in% oth , "others",dataframe$Variant.Effect)
  for(unique_value in unique(dataframe$Variant.Effect)){
    dataframe[paste("Variant.Effect", unique_value, sep = ".")] <- ifelse(dataframe$Variant.Effect == unique_value, 1, 0) %>% as.numeric()
  }

  for(unique_value in unique(dataframe$Functional.Class)){
    dataframe[paste("Functional.Class", unique_value, sep = ".")] <- ifelse(dataframe$Functional.Class == unique_value, 1, 0) %>% as.numeric()
  }

  dataframe$Variant.Effect.Impact[dataframe$Variant.Effect.Impact == "HIGH"] <-3
  dataframe$Variant.Effect.Impact[dataframe$Variant.Effect.Impact == "MODIFIER"]<-0
  dataframe$Variant.Effect.Impact[dataframe$Variant.Effect.Impact == "MODERATE"]<-2
  dataframe$Variant.Effect.Impact[dataframe$Variant.Effect.Impact == "LOW"]<-1
  dataframe$Variant.Effect.Impact <- as.numeric(dataframe$Variant.Effect.Impact)

  dataframe$Seen.as.Somatic <- ifelse(dataframe$Seen.as.Somatic == "Y", 1, 0) %>% as.numeric()

  dataframe$Seen.as.Germline <- ifelse(dataframe$Seen.as.Germline == "Y", 1, 0) %>% as.numeric()

  dataframe$Cancer.Gene.Census <- ifelse(dataframe$Cancer.Gene.Census == "Y", 1, 0) %>% as.numeric()


  dataframe <- dataframe %>%
    select(-c("REF", "ALT", "Variant.Type", "Transcript.Biotype", "Variant.Effect", "Functional.Class",
              "REF_ALT", "FR_REF", "FR_ALT", "FR_REF_ALT"))


  dataframe$Category <- as.factor(dataframe$Category)

  base_features <- colnames(dataframe)

  other_features <- base_features[!base_features %in% all_possible_features]
  dataframe <- dataframe %>% select(-all_of(other_features))
  return(dataframe)
}




