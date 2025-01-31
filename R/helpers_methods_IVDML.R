
# checks if the iv_method is available for the fitted model
check_iv_method <- function(model, iv_method){
  if(!(iv_method %in% model$iv_method)){
    stop(paste("IV method ", iv_method, " is not available.", sep = ""))
  }
}

# checks if input is consistent, i.e. if some parameters for the heterogeneous
# are supplied, all the parameters should be supplied, otherwise only the
# homogeneous treatment effect can be estimated.
check_input_consistency <- function(model, iv_method, a, A, kernel_name, bandwidth, N){
  if(is.null(a)){
    if(!(is.null(A) && is.null(kernel_name) && is.null(bandwidth))){
      warning("The result is for the homogeneous treatment effect. The heterogeneous treatment effect is only possible if !is.null(a).")
    }
  } else {
    if(is.null(kernel_name) || is.null(bandwidth)){
      stop("Kernel name and bandwidth must be specified for the heterogeneous treatment effect.")
    }
    if(is.null(A)){
      if(!is.null(model$A)){
        A <- model$A
      } else {
        stop("A must be specified either in the model or as argument of the function.")
      }
    }
    stopifnot(length(A) == N)
  }
}

# helper function to obtain the kernel
# all kernel's are standardized to have \integral x^2 K(x) dx = 1
get_kernel_function <- function(kernel_name){
  if(kernel_name == "boxcar"){
    Kx0 <- function(x){0.5 * (abs(x) <= 1)}
    A <- 1/sqrt(3)
  } else if(kernel_name == "gaussian"){
    Kx0 <- function(x){exp(-x^2/2)/sqrt(2 * pi)}
    A <- 1
  } else if(kernel_name == "epanechnikov"){
    Kx0 <- function(x){3/4*(1-x^2)*(abs(x)<=1)}
    A <- 1/sqrt(5)
  } else if(kernel_name == "tricube"){
    Kx0 <- function(x){70/81*(1-abs(x)^3)^3*(abs(x)<=1)}
    A <- sqrt(35/243)
  } else {
    stop("Kernel not implemented.")
  }
  Kx <- function(x){A * Kx0(A * x)}
  return(Kx)
}

# helper function to construct a matrix of A values
# according to the cross_fitting indices
construct_Amat <- function(ind, A){
  K_dml <- max(ind)
  N <- length(ind)
  n <- ceiling(N/K_dml)
  Amat <- matrix(0, nrow = n, ncol = K_dml)
  for(k in 1:K_dml){
    test <- which(ind == k)
    Amat[1:length(test), k] <- A[test]
  }
  return(Amat)
}



# helper function to get "effective" instrument residuals (possibly reweighted by kernel)
get_instr_res <- function(result_split, iv_method, a, A, kernel_name, bandwidth){
  instr_res <- switch(iv_method,
                      "linearIV" = result_split$res_Z_X,
                      "mlIV" = result_split$res_f_X,
                      "mlIV_direct" = result_split$res_f_X_direct)
  if(!(is.null(a))){
    Amat <- construct_Amat(result_split$cross_fitting_ind, A)
    Kx <- get_kernel_function(kernel_name)
    KA <- Kx((Amat - a) / bandwidth)
    instr_res <- KA * instr_res
  }
  return(instr_res)
}

# helper function to estimate the coefficient from the residuals
estimate_coef <- function(Y_res, D_res, instr_res){
  return(sum(Y_res * instr_res)/sum(D_res * instr_res))
}

# helper function to estimate the standard error of coefficient estimate from the residuals
estimate_se <- function(Y_res, D_res, instr_res, betahat){
  return(sqrt(sum((Y_res - betahat * D_res)^2 * instr_res^2)/sum(D_res * instr_res)^2))
}

# wrapper function to get the coefficient estimate from the results of a single split
coef_single_split <- function(result_split, iv_method, a, A, kernel_name, bandwidth){
  Y_res <- result_split$res_Y_X
  D_res <- result_split$res_D_X
  instr_res <- get_instr_res(result_split, iv_method, a, A, kernel_name, bandwidth)
  return(estimate_coef(Y_res, D_res, instr_res))
}

# wrapper function to get the standard error estimate from the results of a single split
se_single_split <- function(result_split, iv_method, a, A, kernel_name, bandwidth){
  Y_res <- result_split$res_Y_X
  D_res <- result_split$res_D_X
  instr_res <- get_instr_res(result_split, iv_method, a, A, kernel_name, bandwidth)
  betahat <- estimate_coef(Y_res, D_res, instr_res)
  return(estimate_se(Y_res, D_res, instr_res, betahat))
}

# helper function to get IV strength from the results of a single split
iv_strength_single_split <- function(result_split, iv_method, N){
  D_res <- result_split$res_D_X
  instr_res <- switch(iv_method,
                      "linearIV" = result_split$res_Z_X,
                      "mlIV" = result_split$res_f_X,
                      "mlIV_direct" = result_split$res_f_X_direct)
  return((N-1) * sum(instr_res * D_res)^2/(sum(instr_res^2) * sum(D_res^2) - sum(D_res * instr_res)^2))
}

# median IV strength over the S_split cross fitting splits
iv_strength <- function(model, iv_method){
  check_iv_method(model, iv_method)
  N <- length(model$results_splits[[1]]$cross_fitting_ind)
  return(median(unlist(lapply(model$results_splits, iv_strength_single_split, iv_method = iv_method, N = N))))
}

# helper function to get the robust test statistic with standard error from the results of a single split
robust_test_statistic_single_split <- function(result_split, candidate_value, iv_method, a, A, kernel_name, bandwidth, N_eff){
  Y_res <- result_split$res_Y_X
  D_res <- result_split$res_D_X
  instr_res <- get_instr_res(result_split, iv_method, a, A, kernel_name, bandwidth)
  Qhat <- sum((Y_res - candidate_value * D_res) * instr_res)/N_eff
  if(is.null(bandwidth)){
    SEQhat <- sqrt(sum((Y_res - candidate_value * D_res)^2 * instr_res^2)/N_eff - Qhat^2)
  } else {
    SEQhat <- sqrt(sum((Y_res - candidate_value * D_res)^2 * instr_res^2)/N_eff - bandwidth * Qhat^2)
  }
  return(c(Qhat = Qhat, SEQhat = SEQhat))
}

# helper function to get the robust p-value from the results of a single split
robust_p_value_single_split <- function(result_split, candidate_value, iv_method, a, A, kernel_name, bandwidth, N_eff){
  test_stat <- robust_test_statistic_single_split(result_split, candidate_value, iv_method, a, A, kernel_name, bandwidth, N_eff)
  return(2 * (1-pnorm(sqrt(N_eff)*abs(test_stat[1])/test_stat[2])))
}

# helper function to get the aggregated robust test statistic
robust_test_statistic_aggregated <- function(model, candidate_value, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL){
  check_iv_method(model, iv_method)
  N <- length(model$results_splits[[1]]$cross_fitting_ind)
  check_input_consistency(model, iv_method, a, A, kernel_name, bandwidth, N)
  N_eff <- ifelse(is.null(a), N, N * bandwidth)
  test_stats <- do.call(rbind, lapply(model$results_splits, robust_test_statistic_single_split, candidate_value = candidate_value, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth, N_eff = N_eff))
  Q_agg <- median(test_stats[, 1])
  vars <- test_stats[, 2]^2
  SE_agg <- sqrt(median(vars + (test_stats[, 1] - Q_agg)^2))
  return(c(Q_agg = Q_agg, SE_agg_N = SE_agg/sqrt(N_eff)))
}


