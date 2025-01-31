#' Extract Treatment Effect Estimate from an IVDML Model
#'
#' This function computes the estimated (potentially heterogeneous) treatment effect from a fitted `IVDML` model (output of [fit_IVDML()]).
#'
#' @param model An object of class `IVDML`, produced by the [fit_IVDML()] function.
#' @param iv_method Character. The instrumental variable estimation method to use. Must be one of the methods specified in the fitted model.
#' @param a Numeric (optional). A specific value of `A` at which to evaluate the heterogeneous treatment effect. If `NULL`, the function returns the homogeneous treatment effect.
#' @param A Numeric vector (optional). The variable with respect to which treatment effect heterogeneity is considered. If `NULL`, the function assumes the `A` used in model fitting.
#' @param kernel_name Character (optional). The name of the kernel function to use for smoothing (if a heterogeneous treatment effect is estimated). Needs to be one of "boxcar", "gaussian", "epanechnikov" or "tricube".
#' @param bandwidth Numeric (optional). The bandwidth for the kernel smoothing (if a heterogeneous treatment effect is estimated).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return If `a` is not specified, the estimated homogeneous treatment effect is returned. If `a` is specified, the heterogeneous treatment effect \eqn{\beta(a)} at \eqn{A = a} is returned.
#'
#' @examples
#' set.seed(1)
#' Z <- rnorm(100)
#' X <- Z + rnorm(100)
#' H <- rnorm(100)
#' D <- Z^2 + sin(X) + H + rnorm(100)
#' A <- X
#' Y <- tanh(A) * D + cos(X) - H + rnorm(100)
#' fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = "gam")
#' coef(fit, iv_method = "mlIV")
#' coef(fit, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2)

#'
#' @export
coef.IVDML <- function(model, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, ...){
  check_iv_method(model, iv_method)
  N <- length(model$results_splits[[1]]$cross_fitting_ind)
  if(is.null(A)){
    A <- model$A
  }
  check_input_consistency(model, iv_method, a, A, kernel_name, bandwidth, N)
  return(median(unlist(lapply(model$results_splits, coef_single_split, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth))))
}



se <- function(model, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL){
  check_iv_method(model, iv_method)
  N <- length(model$results_splits[[1]]$cross_fitting_ind)
  check_input_consistency(model, iv_method, a, A, kernel_name, bandwidth, N)
  N_eff <- ifelse(is.null(a), N, N * bandwidth)
  coefs <- unlist(lapply(model$results_splits, coef_single_split, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth))
  med_coef <- median(coefs)
  vars <- N_eff * unlist(lapply(model$results_splits, se_single_split, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth))^2
  return(sqrt(median(vars + (coefs - med_coef)^2)/N_eff))
}





standard_confint <- function(model, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, level = 0.95){
  # the se and coef functions function already checks the input, so we do not have to do it again
  coef_hat <- coef(model, iv_method, a, A, kernel_name, bandwidth)
  se_hat <- se(model, iv_method, a, A, kernel_name, bandwidth)
  z <- -qnorm((1-level)/2)
  lower <- coef_hat - z * se_hat
  upper <- coef_hat + z * se_hat
  if(is.null(a)){
    heterogeneous_parameters <-  NULL
  } else {
    heterogeneous_parameters <- list(a = a, kernel_name = kernel_name, bandwidth = bandwidth)
  }
  return(list(CI = c("lower" = lower, "upper" = upper), level = level, heterogeneous_parameters = heterogeneous_parameters))
}





# aggregate the robust p_values from the splits using the method by Meinshausen, Meier, Bühlmann (2009) or DML type
robust_p_value_aggregated <- function(model, candidate_value, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, agg_method = "DML_agg", gamma = 0.5){
  check_iv_method(model, iv_method)
  N <- length(model$results_splits[[1]]$cross_fitting_ind)
  check_input_consistency(model, iv_method, a, A, kernel_name, bandwidth, N)
  N_eff <- ifelse(is.null(a), N, N * bandwidth)
  if(agg_method == "MMB_agg"){
    p_values <- unlist(lapply(model$results_splits, robust_p_value_single_split, candidate_value = candidate_value, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth, N_eff = N_eff))
    if(length(p_values) == 1){
      return(unname(p_values))
    }
    else{
      return(min(1, quantile(p_values, probs = gamma, names = FALSE)/gamma))
    }
  }
  else if(agg_method == "DML_agg"){
    test_stat <- robust_test_statistic_aggregated(model, candidate_value, iv_method, a, A, kernel_name, bandwidth)
    return(unname(2 * (1-pnorm(abs(test_stat[1])/test_stat[2]))))
  }
  else{
    stop("agg_method must be either DML_agg or MMB_agg.")
  }
}




# agg_method can be DML_agg for the method inspired by DML (TODO: Need to check the aggregation scale for the variance)
# or MMB_agg for the method by Meinshausen, Meier, Bühlmann (2009)
# gamma is only relevant for MMB_agg
robust_confint <- function(model, iv_method, level = 0.95, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, CI_range = NULL, agg_method = "DML_agg", gamma = 0.5){
  # input is checked in coef_estimate
  coef_estimate <- coef(model, iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth)
  if(is.null(CI_range)){
    warning("No CI_range was specified. It is set to a range that is 4 times as large as the standard CI and centered at the point estimate.")
    se <- se(model, iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth)
    z <- -qnorm((1-level)/2)
    CI_range = c(coef_estimate - 4 * z *se, coef_estimate + 4 * z * se)
  }
  if(!((coef_estimate <= CI_range[2])&&(coef_estimate >= CI_range[1]))){
    warning("The point estimate does not lie in the specified CI_range.")
  }
  alpha <- 1 - level
  if(agg_method == "MMB_agg"){
    agg_test <- function(candidate_value){
      return(alpha - robust_p_value_aggregated(model, candidate_value, iv_method, a, A, kernel_name, bandwidth, agg_method = "MMB_agg", gamma))
    }
  }
  else if (agg_method == "DML_agg"){
    agg_test <- function(candidate_value){
      Tstat <- robust_test_statistic_aggregated(model, candidate_value, iv_method, a, A, kernel_name, bandwidth)
      z <- -qnorm((1-level)/2)
      return(abs(Tstat[1]) - z * Tstat[2])
    }
  } else {
    stop("agg_method must be either DML_agg or MMB_agg.")
  }
  ci_left <- NA
  ci_right <- NA
  try(
    {
      tleft <- agg_test(CI_range[1])
      tright <- agg_test(CI_range[2])
      tmin <- optimize(agg_test, interval = CI_range)
      bmin <- tmin$minimum
      tmin <- tmin$objective
      tmax <- optimize(agg_test, interval = CI_range, maximum = TRUE)
      bmax <- tmax$maximum
      tmax <- tmax$objective
      # it might happen that optimize does not find the global minimum if
      # the function is constant in some regions
      # Double check with coef_estimate
      if(agg_test(coef_estimate)<= tmin){
        tmin <- agg_test(coef_estimate)
        bmin <- coef_estimate
      }
      if(tmax <= 0){
        ci_left <- CI_range[1]
        ci_right <- CI_range[2]
        message_text <- "The method returns the whole interval CI_range."
      } else if(tmin > 0){
        ci_left <- CI_range[1]
        ci_right <- CI_range[2]
        message_text <- "The method returns an empty confidence set for method. This might happen if the optimization routine fails. Try using a smaller interval for CI_range."
      } else if((tleft > 0) & (tright > 0)){# hence, the maximum is bigger than zero and the minimum less than zero
        ci_left <- uniroot(agg_test, lower = CI_range[1], upper = bmin)$root
        ci_right <- uniroot(agg_test, lower = bmin, upper = CI_range[2])$root
        message_text <- "The interval is contained in CI_range."
      } else if((tleft > 0) & (tright < 0)){
        ci_left <- uniroot(agg_test, lower = CI_range[1], upper = CI_range[2])$root
        ci_right <- CI_range[2]
        message_text <- "The right endpoint of the interval coincides with the right endpoint of CI_range."
      } else if((tleft < 0) & (tright > 0)){
        ci_left <- CI_range[1]
        ci_right <- uniroot(agg_test, lower = CI_range[1], upper = CI_range[2])$root
        message_text <- "The left endpoint of the interval coincides with the left endpoint of CI_range."
      } else{ # The confidence set is not connected
        ci_left <- CI_range[1]
        ci_right <- CI_range[2]
        r1 <- uniroot(agg_test, lower = CI_range[1], upper = bmax)$root
        r2 <- uniroot(agg_test, lower = bmax, upper = CI_range[2])$root
        message_text <- paste("The method returns a non-connected confidence set. Precisely, it is given by [", CI_range[1], ", ", r1,"] U [", r2, ", ", CI_range[2], "].", sep = "")
      }
    }, silent = TRUE
  )
  if(is.na(ci_left)||is.na(ci_right)){
    ci_left <- CI_range[1]
    ci_right <- CI_range[2]
    message_text <- "The method returns a non-connected confidence set. The precise form cannot be determined."
  }
  if(is.null(a)){
    heterogeneous_parameters <-  NULL
  } else {
    heterogeneous_parameters <- list(a = a, kernel_name = kernel_name, bandwidth = bandwidth)
  }
  return(list(CI = c("lower" = ci_left, "upper" = ci_right), level = level, message = message_text, heterogeneous_parameters = heterogeneous_parameters))
}

# normal reference rule
bandwidth_normal <- function(A){
  Q <- quantile(A, 0.75) - quantile(A, 0.25)
  return(1.06 * min(sd(A), Q/1.34) /length(A)^0.2)
}


