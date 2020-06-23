load(sprintf("%s/avelumab-001-personalis-TNvTO_xome/ml_Yang/mutational_signatures/TO/all_possible_features.Rdata", Sys.getenv("GITCLONES")))

use_data(all_possible_features)


gbm_model <- h2o.loadModel(paste0(sprintf("%s/avelumab-001-personalis-TNvTO_xome/ml_Yang/scripts/final_work/new_R_package/",Sys.getenv("GITCLONES")),
                                  "Grid_GBM_base_data1_sid_af67_1_model_R_1590467726503_19783_model_2"))

use_data(gbm_model)
