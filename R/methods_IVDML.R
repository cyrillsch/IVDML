#' Extract Treatment Effect Estimate from an IVDML Object
#'
#' This function computes the estimated (potentially heterogeneous) treatment effect from a fitted `IVDML` object (output of [fit_IVDML()]).
#'
#' @param object An object of class `IVDML`, produced by the [fit_IVDML()] function.
#' @param iv_method Character. The instrumental variable estimation method to use. Must be one of the methods specified in the fitted object.
#' @param a Numeric (optional). A specific value of `A` at which to evaluate the heterogeneous treatment effect. If `NULL`, the function returns the homogeneous treatment effect.
#' @param A Numeric vector (optional). The variable with respect to which treatment effect heterogeneity is considered. If `NULL`, the function assumes the `A` used in model fitting.
#' @param kernel_name Character (optional). The name of the kernel function to use for smoothing (if a heterogeneous treatment effect is estimated). Needs to be one of "boxcar", "gaussian", "epanechnikov" or "tricube".
#' @param bandwidth Numeric (optional). The bandwidth for the kernel smoothing (if a heterogeneous treatment effect is estimated).
#' @param ... Further arguments passed to or from other methods.
#'
#' @returns If `a` is not specified, the estimated homogeneous treatment effect is returned. If `a` is specified, the heterogeneous treatment effect \eqn{\beta(a)} at \eqn{A = a} is returned.
#'
#' @examples
#' set.seed(1)
#' Z <- rnorm(100)
#' X <- Z + rnorm(100)
#' H <- rnorm(100)
#' D <- Z^2 + sin(X) + H + rnorm(100)
#' A <- X
#' Y <- tanh(A) * D + cos(X) - H + rnorm(100)
#' fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
#' coef(fit, iv_method = "mlIV")
#' coef(fit, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2)

#'
#' @export
coef.IVDML <- function(object, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, ...){
  check_iv_method(object, iv_method)
  N <- length(object$results_splits[[1]]$cross_fitting_ind)
  if(is.null(A)){
    A <- object$A
  }
  check_input_consistency(object, iv_method, a, A, kernel_name, bandwidth, N)
  return(median(unlist(lapply(object$results_splits, coef_single_split, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth))))
}


#' Compute Standard Error for the Treatment Effect Estimate in an IVDML Object
#'
#' This function calculates the standard error of the estimated (potentially heterogeneous) treatment effect from a fitted `IVDML` object (output of [fit_IVDML()]).
#'
#' @param object An object of class `IVDML`, produced by the [fit_IVDML()] function.
#' @param iv_method Character. The instrumental variable estimation method to use. Must be one of the methods specified in the fitted object.
#' @param a Numeric (optional). A specific value of `A` at which to evaluate the standard error of the heterogeneous treatment effect. If `NULL`, the function returns the standard error of the homogeneous treatment effect.
#' @param A Numeric vector (optional). The variable with respect to which treatment effect heterogeneity is considered. If `NULL`, the function assumes the `A` used in model fitting.
#' @param kernel_name Character (optional). The name of the kernel function to use for smoothing (if a heterogeneous treatment effect is estimated). Must be one of "boxcar", "gaussian", "epanechnikov", or "tricube".
#' @param bandwidth Numeric (optional). The bandwidth for the kernel smoothing (if a heterogeneous treatment effect is estimated).
#'
#' @returns A numeric value representing the estimated standard error of the treatment effect estimate. If `a` is not specified, the function returns the standard error of the homogeneous treatment effect. If `a` is specified, it returns the standard error of the heterogeneous treatment effect estimate at \eqn{A = a}.
#'
#'
#' @examples
#' set.seed(1)
#' Z <- rnorm(100)
#' X <- Z + rnorm(100)
#' H <- rnorm(100)
#' D <- Z^2 + sin(X) + H + rnorm(100)
#' A <- X
#' Y <- tanh(A) * D + cos(X) - H + rnorm(100)
#' fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
#' se(fit, iv_method = "mlIV")
#' se(fit, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2)
#'
#' @export
se <- function(object, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL){
  check_iv_method(object, iv_method)
  N <- length(object$results_splits[[1]]$cross_fitting_ind)
  if(is.null(A)){
    A <- object$A
  }
  check_input_consistency(object, iv_method, a, A, kernel_name, bandwidth, N)
  N_eff <- ifelse(is.null(a), N, N * bandwidth)
  coefs <- unlist(lapply(object$results_splits, coef_single_split, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth))
  med_coef <- median(coefs)
  vars <- N_eff * unlist(lapply(object$results_splits, se_single_split, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth))^2
  return(sqrt(median(vars + (coefs - med_coef)^2)/N_eff))
}


#' Compute Standard Confidence Interval for the Treatment Effect Estimate in an IVDML Object
#'
#' This function calculates a standard confidence interval for the estimated (potentially heterogeneous) treatment effect from a fitted `IVDML` object (output of [fit_IVDML()]). The confidence interval is computed using the normal approximation method using the standard error computed by [se()] and the treatment effect estimate from [coef()].
#'
#' @param object An object of class `IVDML`, produced by the [fit_IVDML()] function.
#' @param iv_method Character. The instrumental variable estimation method to use. Must be one of the methods specified in the fitted object.
#' @param a Numeric (optional). A specific value of `A` at which to compute the confidence interval for the heterogeneous treatment effect. If `NULL`, the function returns the confidence interval for the homogeneous treatment effect.
#' @param A Numeric vector (optional). The variable with respect to which treatment effect heterogeneity is considered. If `NULL`, the function assumes the `A` used in object fitting.
#' @param kernel_name Character (optional). The name of the kernel function to use for smoothing (if a heterogeneous treatment effect is estimated). Must be one of "boxcar", "gaussian", "epanechnikov", or "tricube".
#' @param bandwidth Numeric (optional). The bandwidth for the kernel smoothing (if a heterogeneous treatment effect is estimated).
#' @param level Numeric (default: 0.95). The confidence level for the interval (e.g., 0.95 for a 95% confidence interval).
#'
#' @returns description A list containing:
#'   - `CI`: A numeric vector of length 2 with the lower and upper confidence interval bounds.
#'   - `level`: The confidence level used.
#'   - `heterogeneous_parameters`: A list with values of `a`, `kernel_name`, and `bandwidth` (if applicable), or `NULL` if a homogeneous treatment effect is estimated.
#'
#' @examples
#' set.seed(1)
#' Z <- rnorm(100)
#' X <- Z + rnorm(100)
#' H <- rnorm(100)
#' D <- Z^2 + sin(X) + H + rnorm(100)
#' A <- X
#' Y <- tanh(A) * D + cos(X) - H + rnorm(100)
#' fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
#' standard_confint(fit, iv_method = "mlIV")
#' standard_confint(fit, iv_method = "mlIV", a = 0, A = A,
#'                  kernel_name = "boxcar", bandwidth = 0.2, level = 0.95)
#'
#' @export
standard_confint <- function(object, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, level = 0.95){
  if(is.null(A)){
    A <- object$A
  }
  # the se and coef functions function already check the input, so we do not have to do it again
  coef_hat <- coef.IVDML(object, iv_method, a, A, kernel_name, bandwidth)
  se_hat <- se(object, iv_method, a, A, kernel_name, bandwidth)
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


#' Compute Aggregated Robust p-Value for Treatment Effect in an IVDML Object
#'
#' This function calculates an aggregated robust (with respect to weak IV) p-value for testing a candidate treatment effect value in a fitted `IVDML` object (output of [fit_IVDML()]), using either the the standard Double Machine Learning aggregation method ("DML_agg") or the method by Meinshausen, Meier, and Bühlmann (2009) ("MMB_agg") to aggregate the p-values corresponding to the `S_split` cross-fitting sample splits (where `S_split` was an argument of the [fit_IVDML()] function).
#'
#' @param object An object of class `IVDML`, produced by the [fit_IVDML()] function.
#' @param candidate_value Numeric. The candidate treatment effect value to test.
#' @param iv_method Character. The instrumental variable estimation method to use. Must be one of the methods specified in the fitted object.
#' @param a Numeric (optional). A specific value of `A` at which to compute the p-value for the heterogeneous treatment effect. If `NULL`, the function returns the p-value for the homogeneous treatment effect.
#' @param A Numeric vector (optional). The variable with respect to which treatment effect heterogeneity is considered. If `NULL`, the function assumes the `A` used in model fitting.
#' @param kernel_name Character (optional). The name of the kernel function to use for smoothing (if a heterogeneous treatment effect is estimated). Must be one of "boxcar", "gaussian", "epanechnikov", or "tricube".
#' @param bandwidth Numeric (optional). The bandwidth for the kernel smoothing (if a heterogeneous treatment effect is estimated).
#' @param agg_method Character (default: "DML_agg"). The aggregation method for computing the p-value. Options are:
#'   - `"DML_agg"`: Uses the Double Machine Learning (DML) aggregation approach.
#'   - `"MMB_agg"`: Uses the quantile-based aggregation method of Meinshausen, Meier, and Bühlmann (2009).
#' @param gamma Numeric (default: 0.5). Quantile level for the `"MMB_agg"` method. Ignored if `agg_method = "DML_agg"`.
#'
#' @returns The aggregated robust p-value for testing the candidate treatment effect.
#'
#' @references
#' Meinshausen, N., Meier, L., & Bühlmann, P. (2009). *P-values for high-dimensional regression*. Journal of the American Statistical Association, 104(488), 1671–1681.
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
#' robust_p_value_aggregated(fit, candidate_value = 0, iv_method = "mlIV")
#' robust_p_value_aggregated(fit, candidate_value = 0, iv_method = "mlIV",
#'                           a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2)
#'
#' @export
robust_p_value_aggregated <- function(object, candidate_value, iv_method, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, agg_method = "DML_agg", gamma = 0.5){
  check_iv_method(object, iv_method)
  N <- length(object$results_splits[[1]]$cross_fitting_ind)
  check_input_consistency(object, iv_method, a, A, kernel_name, bandwidth, N)
  N_eff <- ifelse(is.null(a), N, N * bandwidth)
  if(agg_method == "MMB_agg"){
    p_values <- unlist(lapply(object$results_splits, robust_p_value_single_split, candidate_value = candidate_value, iv_method = iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth, N_eff = N_eff))
    if(length(p_values) == 1){
      return(unname(p_values))
    }
    else{
      return(min(1, quantile(p_values, probs = gamma, names = FALSE)/gamma))
    }
  }
  else if(agg_method == "DML_agg"){
    test_stat <- robust_test_statistic_aggregated(object, candidate_value, iv_method, a, A, kernel_name, bandwidth)
    return(unname(2 * (1-pnorm(abs(test_stat[1])/test_stat[2]))))
  }
  else{
    stop("agg_method must be either DML_agg or MMB_agg.")
  }
}




#' Compute Robust Confidence Interval for Treatment Effect in an IVDML Object
#'
#' This function computes a robust (with respect to weak IV) confidence interval/confidence set for the estimated treatment effect in a fitted `IVDML` object (output of [fit_IVDML()]). The confidence interval/confidence set is constructed by inverting the robust test from the [robust_p_value_aggregated()] function, which either uses the Double Machine Learning aggregation method (`"DML_agg"`) or the quantile-based method of Meinshausen, Meier, and Bühlmann (2009) (`"MMB_agg"`) to aggregate the p-values corresponding to the `S_split` cross-fitting sample splits (where `S_split` was an argument of the [fit_IVDML()] function).
#'
#' @param object An object of class `IVDML`, produced by the [fit_IVDML()] function.
#' @param iv_method Character. The instrumental variable estimation method to use. Must be one of the methods specified in the fitted object.
#' @param level Numeric (default: 0.95). The confidence level for the confidence interval.
#' @param a Numeric (optional). A specific value of `A` at which to compute the confidence interval for the heterogeneous treatment effect. If `NULL`, the function returns the confidence interval for the homogeneous treatment effect.
#' @param A Numeric vector (optional). The variable with respect to which treatment effect heterogeneity is considered. If `NULL`, the function assumes the `A` used in model fitting.
#' @param kernel_name Character (optional). The name of the kernel function to use for smoothing (if a heterogeneous treatment effect is estimated). Must be one of `"boxcar"`, `"gaussian"`, `"epanechnikov"`, or `"tricube"`.
#' @param bandwidth Numeric (optional). The bandwidth for the kernel smoothing (if a heterogeneous treatment effect is estimated).
#' @param CI_range Numeric vector of length 2 (optional). The search range for the confidence interval. If `NULL`, the function sets `CI_range` to be four times as large as the standard confidence interval centered at the point estimate of the treatment effect.
#' @param agg_method Character (default: `"DML_agg"`). The aggregation method for computing the confidence interval. Options are:
#'   - `"DML_agg"`: Uses the Double Machine Learning (DML) aggregation approach.
#'   - `"MMB_agg"`: Uses the quantile-based aggregation method of Meinshausen, Meier, and Bühlmann (2009).
#' @param gamma Numeric (default: 0.5). Quantile level for the `"MMB_agg"` method. Ignored if `agg_method = "DML_agg"`.
#'
#' @returns A list with the following elements:
#'   - `CI`: A named numeric vector with the lower and upper bounds of the confidence interval.
#'   - `level`: The confidence level used.
#'   - `message`: A message describing the nature of the confidence set (e.g., whether it spans the full range, is non-connected, or is empty due to optimization failure).
#'   - `heterogeneous_parameters`: A list of parameters (`a`, `kernel_name`, `bandwidth`) if a heterogeneous treatment effect is considered; otherwise, `NULL`.
#'
#' @references
#' Meinshausen, N., Meier, L., & Bühlmann, P. (2009). *P-values for high-dimensional regression*. Journal of the American Statistical Association, 104(488), 1671–1681.
#'
#' @examples
#' set.seed(1)
#' Z <- rnorm(100)
#' X <- Z + rnorm(100)
#' H <- rnorm(100)
#' D <- Z^2 + sin(X) + H + rnorm(100)
#' A <- X
#' Y <- tanh(A) * D + cos(X) - H + rnorm(100)
#' fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
#' robust_confint(fit, iv_method = "mlIV", CI_range = c(-10, 10))
#' robust_confint(fit, iv_method = "mlIV", a = 0, A = A,
#'                kernel_name = "boxcar", bandwidth = 0.2, CI_range = c(-10, 10))
#'
#' @export
robust_confint <- function(object, iv_method, level = 0.95, a = NULL, A = NULL, kernel_name = NULL, bandwidth = NULL, CI_range = NULL, agg_method = "DML_agg", gamma = 0.5){
  # input is checked in coef_estimate
  coef_estimate <- coef.IVDML(object, iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth)
  if(is.null(CI_range)){
    warning("No CI_range was specified. It is set to a range that is 4 times as large as the standard CI and centered at the point estimate.")
    se <- se(object, iv_method, a = a, A = A, kernel_name = kernel_name, bandwidth = bandwidth)
    z <- -qnorm((1-level)/2)
    CI_range = c(coef_estimate - 4 * z *se, coef_estimate + 4 * z * se)
  }
  if(!((coef_estimate <= CI_range[2])&&(coef_estimate >= CI_range[1]))){
    warning("The point estimate does not lie in the specified CI_range.")
  }
  alpha <- 1 - level
  if(agg_method == "MMB_agg"){
    agg_test <- function(candidate_value){
      return(alpha - robust_p_value_aggregated(object, candidate_value, iv_method, a, A, kernel_name, bandwidth, agg_method = "MMB_agg", gamma))
    }
  }
  else if (agg_method == "DML_agg"){
    agg_test <- function(candidate_value){
      Tstat <- robust_test_statistic_aggregated(object, candidate_value, iv_method, a, A, kernel_name, bandwidth)
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

#' Compute Bandwidth Using the Normal Reference Rule
#'
#' This function calculates the bandwidth for kernel smoothing using the Normal Reference Rule. The rule is based on Silverman's rule of thumb, which selects the bandwidth as a function of the standard deviation and interquartile range (IQR) of the data. The bandwidth is computed as: \eqn{h = 1.06 \times \min(\mathrm{sd}(A), \mathrm{IQR}(A) / 1.34) / N^{0.2}}, where \eqn{\mathrm{sd}(A)} is the standard deviation of `A`, \eqn{\mathrm{IQR}(A)} is the interquartile range and `N` is the length of `A`.
#'
#' @param A Numeric vector. The data for which the bandwidth is to be computed.
#'
#' @returns A numeric value representing the computed bandwidth.
#'
#' @references
#' Silverman, B. W. (1986). *Density Estimation for Statistics and Data Analysis*.  Chapman & Hall/CRC monographs on statistics and applied probability.  Chapman & Hall.
#'
#' @examples
#' set.seed(1)
#' A <- rnorm(100)
#' bandwidth_normal(A)
#'
#' @export
bandwidth_normal <- function(A){
  Q <- quantile(A, 0.75) - quantile(A, 0.25)
  return(1.06 * min(sd(A), Q/1.34) /length(A)^0.2)
}


#' Print IVDML
#'
#' Print information for an IVDML object.
#' @param x Fitted object of class `IVDML`.
#' @param ... Further arguments passed to or from other methods.
#' @returns No return value, called for side effects
#' @method print IVDML
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
#' print(fit)
#'
#' @export
print.IVDML <- function(x, ...){
  cat("Fitted IVDML object\n")
  cat("Machine learning method: ", x$ml_method, "\n")
  cat("IV methods: ", paste(x$iv_method, collapse = ", "), "\n")
  cat("Number of cross-fitting sample splits: ", length(x$results_splits))
}

