# wrapper for the tune functions
tune_ml <- function(Y, X, ml_method, ml_par){
  if(is.null(X)){
    fit_par <- list()
  } else {
    if(ml_method == "gam"){
      fit_par <- tune_gam(Y, X, ml_par)
    } else if(ml_method == "xgboost"){
      fit_par <- tune_xgboost(Y, X, ml_par)
    } else if(ml_method == "randomForest"){
      fit_par <- tune_randomForest(Y, X, ml_par)
    }
    else {
      stop("ML method not implemented.")
    }
  }
  return(fit_par)
}
