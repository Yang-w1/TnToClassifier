#' Load directory where patient vcf files are in
#'
#' Calculate cosine similarities with COSMIC muational signatures version 2.
#' @param  directory The directory which includes patients' vcf files.
#' @return 30 consine similarities for each patient
#' @export



cosine_sim <- function(directory){
  setwd(directory)
  set.seed(5)
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  vcf_files <- list.files(pattern = ".vcf", full.names = TRUE)
  foo <- sapply(strsplit(vcf_files, "_"), "[", 2)
  ff <- data.frame()
  gg <- as.data.frame(matrix(ncol=0, nrow=96))

  for (i in 1:length(foo)){
    vcfs <- read_vcfs_as_granges(vcf_files[i], foo[i], ref_genome)
    type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
    ff <- rbind(ff, type_occurrences)
    mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
    gg <- cbind(gg, mut_mat)
  }
  sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
  cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
  new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
  cancer_signatures = cancer_signatures[as.vector(new_order),]
  row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
  cos_sim_samples_signatures = cos_sim_matrix(gg, cancer_signatures)
  cos_sim_samples_signatures <- as.data.frame(cos_sim_samples_signatures)
  return(cos_sim_samples_signatures)
}
