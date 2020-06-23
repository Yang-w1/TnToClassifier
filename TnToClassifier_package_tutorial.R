# Introduction of TnToClassifier package.
# Installation:You can use devtools. Then use install_local function.
# For example: install_local("/home/x217659/clones/avelumab-001-personalis-TNvTO_xome/ml_Yang/FalseVariantClassifier", force = TRUE)
# Since we use h2o package to make prediction. Please make sure the right h2o package version installed.  (version 3.28.0.4). You can install it from https://h2o-release.s3.amazonaws.com/h2o/rel-yu/4/index.html




library(TnToClassifier)    ## step 1
### If we only have TO data and want to make prediction. Just use TO_generate function, the input should be the path "TO_output" folder in. tsv files should be in folder named "TO_output".
# For example, if the path of "TO_output" folder is "/home/x217659/clones/avelumab-001-personalis-TNvTO_xome/TO_output", input should be "/home/x217659/clones/avelumab-001-personalis-TNvTO_xome".
# TO_generate function can generate lots of information like "SamppleID" and data_type. As for tsv files name, please make sure sample id is between first "_" and second "_". data_type is between
# sixth "_" and seventh "_". For example:
# "DNA_AE890967_somatic_dna_small_variant_cancer_research_report_preferred.tsv"
TO_TN_data <- TO_generate(sprintf("%s/avelumab-001-personalis-TNvTO_xome", Sys.getenv("GITCLONES")))   ## step 2



# original_data: Import the path where TO_output and TN_output folders exist.
# TO_output folder should be named "TO_output". TN_output folder should be named "TN_output".
# Generate TO and TN total data including extra information like SampleID, data_type(lowpopula.., cancer...), dataset(TO or TN), Category(TP,FP, FN).

# TO_TN_data <- original_data(sprintf("%s/avelumab-001-personalis-TNvTO_xome", Sys.getenv("GITCLONES")))


# In this case example, we just need to extract part of the TO_TN_data for the guide example to save time.
TO_data <- TO_TN_data %>% filter(Dataset == "TO")
set.seed(123)
random_set <- sample(1:nrow(TO_data), size = 10000, replace = FALSE)
sample_data <- TO_TN_data[random_set,]
sample_data <- sample_data[order(sample_data$SampleID),]

# triplet:  generate the triplet change based on chromosome sequence and position.
sample_data <- triplet(sample_data)     ## step 3


# single_nucleotide_change_frequency: calculate the frequency of single nucleotide change for each patient.
sample_data <- single_nucleotide_change_frequency(sample_data)   ## step 4

# triplet_frequency: calculate triplet change frequency for each patient.
sample_data <- triplet_frequency(sample_data)   ## step 5


# cosine_sim: Calculate cosine similarities with COSMIC muational signatures version 2.
# Input is the path where the vcf files are in. It does not in my local machine. But it needs to be provides. In vcf files name, sample id should be between first "_" and second "_"
# for example: "DNA_AG354769_somatic_dna.vcf"
# suppose: cos_sim <- cosine_sim(path)       ## step 6
# cos_sim$SampleID <- rownames(cos_sim)
# sample_data <- merge(sample_data, cos_sim, all.x = TRUE, by = "SampleID")


# If we have cancer type information for each patient. If not, we do not need to care about it.
# suppose: arm <- arm.inforamtion
# sample_data <- merge(sample_data, arm, all.x = TRUE, by = "SampleID)


# feature_manipulation: includes handle missing value, remove useless columns, convert categorical features to numeric...

sample_data <- feature_manipulation(sample_data)    ## step 7

sample_data1 <- sample_data %>% select(-c(SampleID, Protein.Variant,Gene.Symbol ,data_type))  ## step 8


# make_prediction: use GBM model making prediction for variants. Positive class is FP variant.Inputs are dataframe and path of package
prediction <- make_prediction(sample_data1, "/home/x217659/clones/avelumab-001-personalis-TNvTO_xome/ml_Yang/FalseVariantClassifier")    ## step 9
sample_data2 <- cbind(sample_data, prediction$predict)
colnames(sample_data2)[ncol(sample_data2)] <- "prediction"





