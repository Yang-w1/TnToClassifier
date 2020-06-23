#' Load a data frame
#'
#' Make predictions for TO_output using GBM model.
#' @param  dataframe TO output dataframe to be loaded
#' @return A prediction dataframe.
#' @export

make_prediction <- function(dataframe){
  h2o.init(max_mem_size = "8g", nthreads = -1)
  newdata <- as.h2o(dataframe)
  set.seed(123)
  pred <- h2o.predict(object = gbm_model, newdata = newdata)

  pred1 <- as.data.frame(pred)
}


