#' Load directory where TO_output and path of TN_output are in.
#'
#' Generate a dataframe including different classes (TP, FP, FN).
#' @param  directory The directory includes TO_output and TN_output folders. Also, the folder name of TO_output should be
#' "TO_output". The folder name of TN_output should be "TN_output".
#' @return A dataframe including different classes(TP, FP, FN), datatype, dataset and SampleID.
#' @export
#'
#' @examples
#' original_data(sprintf("%s/avelumab-001-personalis-TNvTO_xome", Sys.getenv("GITCLONES")))


original_data <- function(directory){
  setwd(directory)
  mydir = "TO_output"
  TO_files = list.files(path=mydir, pattern="*.tsv", full.names=TRUE)

  foo<-sapply(strsplit(TO_files, "_"), "[", 3)
  data_type<-sapply(strsplit(TO_files, "_"), "[", 8)


  TO_data = lapply(TO_files, read.table, sep = "\t", header=T, stringsAsFactors = F)


  for(i in seq_along(TO_data)){
    TO_data[[i]]$SampleID <- rep(foo[i],nrow(TO_data[[i]]))
  }


  for(i in seq_along(TO_data)){
    TO_data[[i]]$data_type <- rep(data_type[i],nrow(TO_data[[i]]))
  }


  for (i in seq_along(TO_data)){

    TO_data[[i]]$Sequence <- as.character(TO_data[[i]]$Sequence)
    TO_data[[i]]$ALT<-as.character(TO_data[[i]]$ALT)
    TO_data[[i]]$ALT[TO_data[[i]]$ALT == "TRUE"] <- "T"
    TO_data[[i]]$REF<-as.character(TO_data[[i]]$REF)
    TO_data[[i]]$REF[TO_data[[i]]$REF == "TRUE"] <- "T"
  }


  mydir2<- "TN_output"

  TN_files<-list.files(path=mydir2, pattern="*.tsv", full.names=TRUE)

  TN_sampleID<-sapply(strsplit(TN_files, "_"), "[", 3)
  TN_data_type <- sapply(strsplit(TN_files, "_"), "[", 8)

  TN_Data<-lapply(TN_files, read.table, sep = "\t", header=T, stringsAsFactors = F)

  for(i in seq_along(TN_Data)){
    TN_Data[[i]]$SampleID <- rep(TN_sampleID[i],nrow(TN_Data[[i]]))
  }
  for(i in seq_along(TN_Data)){
    TN_Data[[i]]$data_type <- rep(TN_data_type [i],nrow(TN_Data[[i]]))
  }


  for (i in seq_along(TN_Data)){

    TN_Data[[i]]$Sequence <- as.character(TN_Data[[i]]$Sequence)
    TN_Data[[i]]$ALT<-as.character(TN_Data[[i]]$ALT)
    TN_Data[[i]]$ALT[TN_Data[[i]]$ALT == "TRUE"] <- "T"
    TN_Data[[i]]$REF<-as.character(TN_Data[[i]]$REF)
    TN_Data[[i]]$REF[TN_Data[[i]]$REF == "TRUE"] <- "T"
  }


  TN_Data<-data.table::rbindlist(TN_Data)
  TN_Data<-as_tibble(TN_Data)
  TN_Data <- na_if(TN_Data, "")
  TN_Data<-unique(TN_Data)

  TO_data<-data.table::rbindlist(TO_data)%>%as_tibble()
  TO_data <- na_if(TO_data, "")
  TO_data<-unique(TO_data)


  common_parentID<-intersect(unique(TN_Data$SampleID),unique(TO_data$SampleID))



  temp1<-TN_Data%>%filter(.$SampleID %in% common_parentID)
  temp1$Dataset=rep_len("T/N",length(temp1$Variant.ID))
  temp2<-TO_data%>%filter(.$SampleID %in% common_parentID)
  temp2$Dataset=rep_len("TO",length(temp2$Variant.ID))

  combined_TN_TO<-rbind(temp1,temp2)



  TP <- lapply(temp1$SampleID %>% unique(), function(x){
    a<-temp1%>%filter(.$SampleID == x)
    a <- a[!duplicated(a$Variant.ID),]
    b<-temp2%>%filter(.$SampleID == x)
    b <- b[!duplicated(b$Variant.ID),]
    e <- a %>% select(-c(Variant.ID, data_type, Dataset))
    f <- b %>% select(-c(Variant.ID, data_type, Dataset))
    c <- merge(e,f)
    a <- merge(a, c)
    b <- merge(b, c)
    d<-rbind(a,b)
    return(d)
  }
  )

  TP<-data.table::rbindlist(TP)
  TP$Category<-rep_len("TP",length(TP$Variant.ID))


  FP <- lapply(temp1$SampleID %>% unique(), function(x){
    a<-temp1%>%filter(.$SampleID == x)
    a <- a[!duplicated(a$Variant.ID),]
    b<-temp2%>%filter(.$SampleID == x)
    b <- b[!duplicated(b$Variant.ID),]
    e <- a %>% select(-c(Variant.ID, data_type, Dataset))
    f <- b %>% select(-c(Variant.ID, data_type, Dataset))
    c <- merge(e,f)
    g <- merge(b, c)
    d <- anti_join(b, g)
    return(d)
  }
  )

  FP<-data.table::rbindlist(FP)
  FP$Category<-rep_len("FP",length(FP$Variant.ID))



  FN <- lapply(temp1$SampleID %>% unique(), function(x){
    a<-temp1%>%filter(.$SampleID == x)
    a <- a[!duplicated(a$Variant.ID),]
    b<-temp2%>%filter(.$SampleID == x)
    b <- b[!duplicated(b$Variant.ID),]
    e <- a %>% select(-c(Variant.ID, data_type, Dataset))
    f <- b %>% select(-c(Variant.ID, data_type, Dataset))
    c <- merge(e,f)
    g <- merge(a, c)
    d <- anti_join(a, g)
    return(d)
  }
  )

  FN<-data.table::rbindlist(FN)
  FN$Category<-rep_len("FN",length(FN$Variant.ID))

  TO_TN<-rbind(TP,FP)

  TO_TN<-rbind(TO_TN,FN)
  return(TO_TN)

}



#' Load directory where TO_output is in. (use this function is only TO_output is available)
#'
#' Generate a dataframe including different classes.
#' @param  directory The directory includes TO_output folder. Also, the folder name of TO_output should be
#' "TO_output".
#' @return A dataframe including different classes, datatype, dataset and SampleID.
#' @export
#'
#' @examples
#' original_data(sprintf("%s/avelumab-001-personalis-TNvTO_xome", Sys.getenv("GITCLONES")))


TO_generate <- function(directory){
  setwd(directory)
  mydir = "TO_output"
  TO_files = list.files(path=mydir, pattern="*.tsv", full.names=TRUE)

  foo<-sapply(strsplit(TO_files, "_"), "[", 3)
  data_type<-sapply(strsplit(TO_files, "_"), "[", 8)


  TO_data = lapply(TO_files, read.table, sep = "\t", header=T, stringsAsFactors = F)


  for(i in seq_along(TO_data)){
    TO_data[[i]]$SampleID <- rep(foo[i],nrow(TO_data[[i]]))
  }


  for(i in seq_along(TO_data)){
    TO_data[[i]]$data_type <- rep(data_type[i],nrow(TO_data[[i]]))
  }


  for (i in seq_along(TO_data)){

    TO_data[[i]]$Sequence <- as.character(TO_data[[i]]$Sequence)
    TO_data[[i]]$ALT<-as.character(TO_data[[i]]$ALT)
    TO_data[[i]]$ALT[TO_data[[i]]$ALT == "TRUE"] <- "T"
    TO_data[[i]]$REF<-as.character(TO_data[[i]]$REF)
    TO_data[[i]]$REF[TO_data[[i]]$REF == "TRUE"] <- "T"
  }


  TO_data<-data.table::rbindlist(TO_data)%>%as_tibble()
  TO_data <- na_if(TO_data, "")
  TO_data<-unique(TO_data)

  TO_data2 <- lapply(TO_data$SampleID %>% unique, function(x){
    a <- TO_data %>% filter(.$SampleID == x)
    a <- a[!duplicated(a$Variant.ID),]
    return(a)
  }
  )

  TO_data3 <- data.table::rbindlist(TO_data2)
  TO_data3$Dataset <- "TO"
  TO_data3$Category <- "Unknown"
  return(TO_data3)
}
