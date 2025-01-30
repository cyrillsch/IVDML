#### helper functions for the nuisance function estimates ####


## gam ##

## no tuning
tune_gam <- function(Y, X, ml_par){
  return(ml_par)
}

# input: Y, X and list fit_par with vector ind_linear with indices such that X_j is modeled as linear
# output: fitted gam object of Y vs. X, where
fit_gam <- function(Y, X, fit_par){
  p <- NCOL(X)
  dat <- data.frame(Y,X)
  varnames <- character(length = p)
  for(j in 1:p){
    varnames[j] <- paste("X", j, sep = "")
  }
  colnames(dat) <- c("Y", varnames)
  unique_lengths <- apply(as.matrix(X), 2, function(x){length(unique(x))})
  max_k_gam <- pmin(30, unique_lengths - 1)
  if(all(max_k_gam == 0)){
    form <- "Y ~ 1"
  } else {
    form <- "Y ~ "
    for(j in 1:p){
      if((max_k_gam[j] %in% c(1, 2)) | (j %in% fit_par$ind_linear)){
        form <- paste(form, " + X", j, sep = "")
      }
      else if(max_k_gam[j] > 2){
        form <- paste(form, " + s(X", j, ", k = ", max_k_gam[j], ")", sep = "")
      }
    }
  }
  mod <- mgcv::gam(formula = formula(form), data = dat)
  return(mod)
}

# predict a fitted gam object at observations Xnew
predict_gam <- function(Xnew, mod){
  p <- NCOL(Xnew)
  varnames <- character(length = p)
  for(j in 1:p){
    varnames[j] <- paste("X", j, sep = "")
  }
  Xnew <- data.frame(Xnew)
  colnames(Xnew) <- varnames
  pred <- predict(mod, newdata = Xnew)
  return(pred)
}


## xgboost ##

# input: Y, X, list ml_par of hyperparameters or empty
# output: list containing hyperparameters optimized using cross_validation
tune_xgboost <- function(Y, X, ml_par){
  if(!exists("max_nrounds", ml_par)){
    ml_par$max_nrounds <- 500
  }
  if(!exists("eta", ml_par)){
    ml_par$eta <- c(0.1, 0.2, 0.3, 0.5)
  }
  if(!exists("max_depth", ml_par)){
    ml_par$max_depth <-  c(1, 2, 3, 4, 5, 6, 7)
  }
  if(!exists("early_stopping_rounds", ml_par)){
    ml_par$early_stopping_rounds <- 10
  }
  if(!exists("k_cv", ml_par)){
    ml_par$k_cv <- 10
  }
  par_grid <- with(ml_par, expand.grid(eta = eta, max_depth = max_depth))
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = Y, nthread = 1)
  get_min_mse_nrounds <- function(pars){
    fit_cv <- xgboost::xgb.cv(params = list(nthread = 1, eta = pars[1], max_depth = pars[2]),
                     data = dtrain, nrounds = ml_par$max_nrounds, nfold = ml_par$k_cv,
                     early_stopping_rounds = ml_par$early_stopping_rounds, verbose = FALSE)
    min_mse <- min(fit_cv$evaluation_log$test_rmse_mean)
    indopt <- which(fit_cv$evaluation_log$test_rmse_mean == min_mse)
    if(length(indopt) > 1){
      indopt <- indopt[1]
    }
    return(c(min_mse, indopt))
  }
  par_mses <- apply(par_grid, 1, get_min_mse_nrounds)
  ind_min <- which(par_mses[1,] == min(par_mses[1,]))
  if(length(ind_min) > 1){
    ind_min <- ind_min[1]
  }
  par_min <- par_grid[ind_min,]
  min_nrounds <- par_mses[2, ind_min]
  if(min_nrounds == ml_par$max_nrounds){
    warning("CV chooses nrounds equal to max_nrounds. Consider increasing max_nrounds.")
  }
  return(list(eta = par_min$eta, max_depth = par_min$max_depth, nrounds = min_nrounds))
}

# input: Y, X and list fit_par with values for eta, max_depth and nrounds
# output: fitted xgboost object of Y vs. X
fit_xgboost <- function(Y, X, fit_par){
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = Y, nthread = 1)
  mod <-  xgboost::xgboost(params = list(nthread = 1, eta =  fit_par$eta, max_depth = fit_par$max_depth),
                  data = dtrain, nrounds = fit_par$nrounds, verbose = FALSE)
  return(mod)
}

# predict a fitted xgboost object at observations Xnew
predict_xgboost <- function(Xnew, mod){
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(Xnew), nthread = 1)
  pred <- predict(mod, dtest)
  return(pred)
}

# input: Y, X, list ml_par of hyperparameters or empty
# output: list containing hyperparameters optimized using out of bag error.
tune_randomForest <- function(Y, X, ml_par){
  if(!exists("num.trees", ml_par)){
    ml_par$num.trees  <- 250
  }
  # number of mtry values to try, default 15
  if(!exists("num_mtry", ml_par)){
    ml_par$num_mtry <- 15
  }
  if(!exists("mtry", ml_par)){
    # default: try some selection of num_mtry integers in [min(pX/3, sqrt(pX)), 2*pX/3]
    pX <- NCOL(X)
    ml_par$mtry <- unique(round(seq(max(1, floor(min(sqrt(pX),pX/3))), ceiling(2*pX/3), length.out = ml_par$num_mtry)))
  }
  # by default, we do not restrict max_depth, but regularize min_node_size
  if(!exists("max.depth", ml_par)){
    ml_par$max.depth <-  0
  }
  # number of min.node.size values to try, default: 15
  if(!exists("num_min.node.size", ml_par)){
    ml_par$num_min.node.size <- 15
  }
  # default: try some selection  of num_min.node.size integers from a logarithmic grid up to nX/4
  if(!exists("min.node.size", ml_par)){
    nX <- NROW(X)
    ml_par$min.node.size <- unique(round(5 * exp(seq(0, log(nX/20), length.out = ml_par$num_min.node.size))))
  }
  par_grid <- with(ml_par, expand.grid(num.trees = num.trees, mtry = mtry, max.depth = max.depth, min.node.size = min.node.size))
  p <- NCOL(X)
  dat <- data.frame(Y,X)
  varnames <- character(length = p)
  for(j in 1:p){
    varnames[j] <- paste("X", j, sep = "")
  }
  colnames(dat) <- c("Y", varnames)
  form <- "Y ~ X1"
  if(p > 1){
    for(j in 2:p){
      form <- paste(form, " + X", j, sep = "")
    }
  }
  oob_forest <- function(pars){
    rf <- ranger::ranger(formula = formula(form), data = dat, num.trees = pars[1], mtry = pars[2], max.depth = pars[3], min.node.size = pars[4], write.forest = FALSE)
    return(rf$prediction.error)
  }
  oob_errors <- apply(par_grid, 1, oob_forest)
  ind_min <- which(oob_errors == min(oob_errors))
  if(length(ind_min) > 1){
    ind_min <- ind_min[1]
  }
  par_min <- par_grid[ind_min,]
  return(par_min)
}


# input: Y, X and list fit_par with values for num.trees, mtry, max.depth and min.node.size
# output: fitted ranger object of Y vs. X
fit_randomForest <- function(Y, X, fit_par){
  p <- NCOL(X)
  dat <- data.frame(Y,X)
  varnames <- character(length = p)
  for(j in 1:p){
    varnames[j] <- paste("X", j, sep = "")
  }
  colnames(dat) <- c("Y", varnames)
  form <- "Y ~ X1"
  if(p > 1){
    for(j in 2:p){
      form <- paste(form, " + X", j, sep = "")
    }
  }
  mod <- ranger::ranger(formula = formula(form), data = dat, num.trees = fit_par$num.trees, mtry = fit_par$mtry, max.depth = fit_par$max.depth, min.node.size = fit_par$min.node.size)
  return(mod)
}

# predict a fitted ranger object at observations Xnew
predict_randomForest <- function(Xnew, mod){
  p <- NCOL(Xnew)
  varnames <- character(length = p)
  for(j in 1:p){
    varnames[j] <- paste("X", j, sep = "")
  }
  Xnew <- data.frame(Xnew)
  colnames(Xnew) <- varnames
  pred <- predict(mod, data = Xnew)$predictions
  return(pred)
}


## wrapper functions ##

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

# wrapper for the fit functions
fit_ml <- function(Y, X, ml_method, fit_par){
  if(is.null(X)){
    mod <- list(meanY = mean(Y))
  } else {
    if(ml_method == "gam"){
      mod <- fit_gam(Y, X, fit_par)
    } else if(ml_method == "xgboost"){
      mod <- fit_xgboost(Y, X, fit_par)
    } else if(ml_method == "randomForest"){
      mod <- fit_randomForest(Y, X, fit_par)
    }
    else {
      stop("ML method not implemented.")
    }
  }
  return(mod)
}

# wrapper for the predict functions
predict_ml <- function(Xnew, mod, ml_method){
  if(is.null(Xnew)){
    pred <- mod$meanY
  } else {
    if(ml_method == "gam"){
      pred <- predict_gam(Xnew, mod)
    } else if(ml_method == "xgboost"){
      pred <- predict_xgboost(Xnew, mod)
    } else if(ml_method == "randomForest"){
      pred <- predict_randomForest(Xnew, mod)
    }
    else {
      stop("ML method not implemented.")
    }
  }
  return(pred)
}


