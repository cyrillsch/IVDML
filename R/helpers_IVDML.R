


tuning_helper <- function(Y, D, Z, X, ml_method, ml_par, iv_method){
  sublist_names <- c("ml_par_D_XZ", "ml_par_D_X", "ml_par_f_X", "ml_par_Y_X", "ml_par_Z_X")
  if(all(names(ml_par) %in% sublist_names)){
    # if at least one of the individual ml parameter lists is defined or ml_par is empty
    if(!(exists("ml_par_D_XZ", ml_par))){
      ml_par$ml_par_D_XZ <- list()
    }
    if(!(exists("ml_par_D_X", ml_par))){
      ml_par$ml_par_D_X <- list()
    }
    if(!(exists("ml_par_Y_X", ml_par))){
      ml_par$ml_par_Y_X <- list()
    }
    if(!(exists("ml_par_Z_X", ml_par))){
      ml_par$ml_par_Z_X <- list()
    }
    if(!(exists("ml_par_f_X", ml_par))){
      ml_par$ml_par_f_X <- list()
    }
  }
  else if(any(names(ml_par) %in% sublist_names)){
    stop(paste("Either specify the ml parameters globally or via seperate lists, but not mixed."))
  }
  else {
    # no individual regression parameters are defined, hence they should be all the same
    ml_par <- list("ml_par_D_XZ" = ml_par, "ml_par_D_X" = ml_par, "ml_par_f_X" = ml_par, "ml_par_Y_X" = ml_par, "ml_par_Z_X" = ml_par)
  }
  # check if all the names in the ml_par lists are actually defined
  gam_par_names <- c("ind_lin_X", "ind_lin_Z")
  xgb_par_names <- c("max_nrounds", "k_cv", "early_stopping_rounds", "eta", "max_depth")
  rf_par_names <- c("num.trees", "num_mtry", "mtry", "max.depth", "num_min.node.size", "min.node.size")
  check_name <- function(sublist_name, par_names){
    ml_par_sublist <- ml_par[[sublist_name]]
    if(length(ml_par_sublist) > 0){
      for(i in 1:length(ml_par_sublist)){
        if(!(names(ml_par_sublist)[i] %in% par_names)){
          warning(paste(names(ml_par_sublist)[i], " is not a permissible parameter for your chosen ml_method."))
        }
      }
    }
  }
  switch(ml_method,
         "gam" = sapply(sublist_names, check_name, par_names = gam_par_names),
         "xgboost" = sapply(sublist_names, check_name, par_names = xgb_par_names),
         "randomForest" = sapply(sublist_names, check_name, par_names = rf_par_names))

  # if ML method is gam, for every ml_par_D_XZ, ..., ml_par_Z_X, we only want one argument ind_linear that gives the
  # indices of the corresponding predictors (X, Z or (X,Z)) that should be modeled as linear
  if(ml_method == "gam"){
    if(exists("ind_lin_X", ml_par$ml_par_D_X)){
      ml_par$ml_par_D_X <- list(ind_linear= ml_par$ml_par_D_X$ind_lin_X)
    }
    if(exists("ind_lin_X", ml_par$ml_par_Y_X)){
      ml_par$ml_par_Y_X <- list(ind_linear = ml_par$ml_par_Y_X$ind_lin_X)
    }
    if(exists("ind_lin_X", ml_par$ml_par_f_X)){
      ml_par$ml_par_f_X <- list(ind_linear = ml_par$ml_par_f_X$ind_lin_X)
    }
    if(exists("ind_lin_X", ml_par$ml_par_Z_X)){
      ml_par$ml_par_Z_X <- list(ind_linear = ml_par$ml_par_Z_X$ind_lin_X)
    }
    if(exists("ind_lin_X", ml_par$ml_par_D_XZ)){
      ind_lin <- ml_par$ml_par_D_XZ$ind_lin_X
    } else {
      ind_lin <- NULL
    }
    if(exists("ind_lin_Z", ml_par$ml_par_D_XZ)){
      ind_lin <- c(ind_lin, ml_par$ml_par_D_XZ$ind_lin_Z + NCOL(X))
    }
    ml_par$ml_par_D_XZ <- list(ind_linear = ind_lin)
  }
  tuned_ml_par <- list()
  # Y vs. X
  tuned_ml_par$ml_par_Y_X <- tune_ml(Y = Y, X = X, ml_method = ml_method, ml_par = ml_par$ml_par_Y_X)
  # D vs. X
  tuned_ml_par$ml_par_D_X <- tune_ml(Y = D, X = X, ml_method = ml_method, ml_par = ml_par$ml_par_D_X)

  if("linearIV" %in% iv_method){
    # Z vs. X
    tuned_ml_par$ml_par_Z_X <- lapply(1:NCOL(Z), function(j){tune_ml(Y = Z[, j], X = X, ml_method = ml_method, ml_par = ml_par$ml_par_Z_X)})
  }
  if(("mlIV" %in% iv_method)|("mlIV_direct" %in% iv_method)){
    # D vs. (X, Z)
    tuned_ml_par$ml_par_D_XZ <- tune_ml(Y = D, X = cbind(X, Z), ml_method = ml_method, ml_par = ml_par$ml_par_D_XZ)
  }
  if("mlIV" %in% iv_method){
    # f(X, Z) vs. X
    if(ml_method == "gam"){
      # we do not need to actually fit the model if ml_method == gam
      tuned_ml_par$ml_par_f_X <- tuned_ml_par$ml_par_D_X
    } else {
      fit0_f <- fit_ml(Y = D, X = cbind(X, Z), ml_method = ml_method, fit_par = tuned_ml_par$ml_par_D_XZ)
      fitted0_f <- predict_ml(Xnew = cbind(X, Z), mod = fit0_f, ml_method = ml_method)
      tuned_ml_par$ml_par_f_X <- tune_ml(Y = fitted0_f, X = X, ml_method = ml_method, ml_par = ml_par$ml_par_f_X)
    }
  }
  return(tuned_ml_par)
}

residuals_helper <- function(Y, D, Z, X, ml_method, tuned_ml_par, iv_method, K_dml){
  N <- length(Y)
  n <- ceiling(N/K_dml)
  dZ <- NCOL(Z)
  ind <- rep(1:K_dml, length.out = N)
  ind <- sample(ind, replace = FALSE)
  # Prepare matrices for cross-fitting residuals
  # Y vs. X
  res_Y_X <- matrix(0, nrow = n, ncol = K_dml)
  # D vs. X
  res_D_X <- matrix(0, nrow = n, ncol = K_dml)
  if("linearIV" %in% iv_method){
    # Z vs. X
    res_Z_X <- matrix(0, nrow = n, ncol = K_dml)
  }
  if("mlIV" %in% iv_method){
    # f(X, Z) vs. X
    res_f_X <- matrix(0, nrow = n, ncol = K_dml)
  }
  if("mlIV_direct" %in% iv_method){
    res_f_X_direct <- matrix(0, nrow = n, ncol = K_dml)
  }
  ## estimating the cross-fitted residuals ##
  for(k in 1:K_dml){
    test <- which(ind == k)
    train <- which(ind != k)
    Xtrain <- if(is.null(X)) NULL else as.matrix(X[train,])
    Xtest <- if(is.null(X)) NULL else as.matrix(X[test, ])
    Ztrain <- as.matrix(Z[train,])
    Ztest <- as.matrix(Z[test, ])
    ltest <- length(test)
    # Y vs. X
    fit_Y_X <- fit_ml(Y = Y[train], X = Xtrain, ml_method = ml_method, fit_par = tuned_ml_par$ml_par_Y_X)
    res_Y_X[1:ltest, k] <- Y[test] - predict_ml(Xnew = Xtest, mod = fit_Y_X, ml_method = ml_method)
    # D vs. X
    fit_D_X <- fit_ml(Y = D[train], X = Xtrain, ml_method = ml_method, fit_par = tuned_ml_par$ml_par_D_X)
    pred_phi_test <- predict_ml(Xnew = Xtest, mod = fit_D_X, ml_method = ml_method)
    res_D_X[1:ltest, k] <- D[test] - pred_phi_test
    if("linearIV" %in% iv_method){
      # Z vs. X
      if(dZ == 1){
        fit_Z_X <- fit_ml(Y = Ztrain, X = Xtrain, ml_method = ml_method, fit_par = tuned_ml_par$ml_par_Z_X[[1]])
        res_Z_X[1:ltest, k] <- as.numeric(Ztest) - predict_ml(Xnew = Xtest, mod = fit_Z_X, ml_method = ml_method)
      } else {
        # we can use the projection of res_D_X onto res_Z_X as res_Z_X
        res_Z_X_k <- matrix(NA, nrow = length(test), ncol = dZ)
        for(j in 1:dZ){
          fit_Z_X_j <- fit_ml(Y = Ztrain[, j], X = Xtrain, ml_method = ml_method, fit_par = tuned_ml_par$ml_par_Z_X[[j]])
          # res_Z_X_mult[1:ltest, k, j] <- as.numeric(Ztest[, j]) - predict_ml(Xnew = Xtest, mod = fit_Z_X_j, ml_method = ml_method)
          res_Z_X_k[ , j] <- as.numeric(Ztest[, j]) - predict_ml(Xnew = Xtest, mod = fit_Z_X_j, ml_method = ml_method)
        }
        res_Z_X[1:ltest, k] <- fitted(lm(res_D_X[1:ltest, k] ~ -1 + res_Z_X_k))
      }
    }
    if(("mlIV" %in% iv_method)|("mlIV_direct" %in% iv_method)){
      # D vs. (X, Z)
      fit_D_XZ <- fit_ml(Y = D[train], X = cbind(Xtrain, Ztrain), ml_method = ml_method, fit_par = tuned_ml_par$ml_par_D_XZ)
      pred_f_test <- predict_ml(Xnew = cbind(Xtest, Ztest), mod = fit_D_XZ, ml_method = ml_method)
    }
    if("mlIV" %in% iv_method){
      # f(X, Z) vs. X
      fitted_f <- predict_ml(Xnew = cbind(Xtrain, Ztrain), mod = fit_D_XZ, ml_method = ml_method)
      fit_f_X <- fit_ml(Y = fitted_f, X = Xtrain, ml_method = ml_method, fit_par = tuned_ml_par$ml_par_f_X)
      res_f_X[1:ltest, k] <- pred_f_test - predict_ml(Xnew = Xtest, mod = fit_f_X, ml_method = ml_method)
    }
    if("mlIV_direct" %in% iv_method){
      res_f_X_direct[1:ltest, k] <- pred_f_test - pred_phi_test
    }
  }
  res_list <- list(res_Y_X = res_Y_X, res_D_X = res_D_X)
  if("linearIV" %in% iv_method){
    res_list$res_Z_X <- res_Z_X
  }
  if("mlIV" %in% iv_method){
    res_list$res_f_X <- res_f_X
  }
  if("mlIV_direct" %in% iv_method){
    res_list$res_f_X_direct <- res_f_X_direct
  }
  res_list$cross_fitting_ind = ind
  return(res_list)
}

