#' Load a data frame
#'
#' Make predictions for TO_output using GBM model.
#' @param  dataframe TO output dataframe to be loaded
#' @param  package_path The path of TnToClassifier package you have downloaded.
#' @return A prediction dataframe.
#' @export

make_prediction <- function(dataframe, package_path){
  h2o.init(max_mem_size = "8g", nthreads = -1)
  newdata <- as.h2o(dataframe)
  set.seed(123)
  model <- h2o.loadModel(paste0(package_path, "/data/","Grid_GBM_base_data1_sid_af67_1_model_R_1590467726503_19783_model_2"))
  pred <- h2o.predict(object = model, newdata = newdata)

  pred1 <- as.data.frame(pred)
}


