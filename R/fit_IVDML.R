#' Fitting Double Machine Learning Models with Instrumental Variables and Potentially Heterogenous Treatment Effect
#'
#' This function is used to fit a Double Machine Learning (DML) model with Instrumental Variables (IV) with the goal to perform inference on potentially heterogeneous treatment effects. The model under study is \eqn{Y = \beta(A)D + g(X) + \epsilon}, where the error \eqn{\epsilon} is potentially correlated with the treatment \eqn{D}, but there is an IV \eqn{Z} satisfying \eqn{\mathbb E[\epsilon|Z,X] = 0}. The object of interest is the treatment effect \eqn{\beta} of the treatment \eqn{D} on the response \eqn{Y}. The treatment effect \eqn{\beta} is either constant or can depend on the univariate quantity \eqn{A}, which is typically a component of the covariates \eqn{X}.
#'
#' @param Y Numeric vector. Response variable.
#' @param D Numeric vector. Treatment variable.
#' @param Z Matrix, vector, or data frame. Instrumental variables.
#' @param X Matrix, vector, or data frame. Additional covariates (default: NULL).
#' @param A Numeric vector. Variable with respect to which treatment effect heterogeneity is considered. Usually equal to a column of X and in this case it can also be specified later (default: NULL).
#' @param ml_method Character. Machine learning method to use. Options are "gam", "xgboost", and "randomForest".
#' @param ml_par List. Parameters for the machine learning method:
#'   - If `ml_method == "gam"`, can specify `ind_lin_Z` and `ind_lin_X` for components of `Z` and `X` to be modeled linearly.
#'   - If `ml_method == "xgboost"`, can specify `max_nrounds`, `k_cv`, `early_stopping_rounds`, and vectors `eta` and `max_depth`.
#'   - If `ml_method == "randomForest"`, can specify `num.trees`, `num_mtry` (number of different mtry values to try out) or a vector `mtry`, a vector `max.depth`, `num_min.node.size` (number of different min.node.size values to try out) or a vector `min.node.size`.
#'   - To specify different parameters for the different nuisance function regressions, `ml_par` should be a list of lists: `ml_par_D_XZ` (parameters for nuisance function \eqn{\mathbb E[D|Z, X]}, needed for `iv_method` "mlIV" and "mlIV_direct"), `ml_par_D_X` (parameters for nuisance function \eqn{\mathbb E[D|X]}, needed for `iv_method` "linearIV", "mlIV" and "mlIV_direct"), `ml_par_f_X` (parameters for nuisance function \eqn{\mathbb E[\widehat{\mathbb E}[D|Z, X]|X]}, needed for `iv_method` "mlIV"), `ml_par_Y_X` (parameters for nuisance function \eqn{\mathbb E[Y|X]}, needed for `iv_method` "linearIV", "mlIV" and "mlIV_direct"), `ml_par_Z_X` (parameters for nuisance function \eqn{\mathbb E[Z|X]}, needed for `iv_method` "linearIV").
#' @param A_deterministic_X Logical. Whether `A` is a deterministic function of `X` (default: TRUE).
#' @param K_dml Integer. Number of cross-fitting folds (default: 5).
#' @param iv_method Character vector. Instrumental variables estimation method. Options:
#' "linearIV", "mlIV", "mlIV_direct" (default: c("linearIV", "mlIV")). "linearIV" corresponds to using instruments linearly and "mlIV" corresponds to using machine learning instruments. "mlIV_direct" is a variant of "mlIV" that uses the same estimate of \eqn{\mathbb E[D|X]} for both the residuals \eqn{X - \mathbb E[D|X]} and \eqn{\mathbb E[D|Z, X] - \mathbb E[D|X]}, whereas "mlIV" uses a two-stage estimate of \eqn{\mathbb E[\widehat{\mathbb E}[D|Z, X]|X]} for the residuals \eqn{\mathbb E[D|Z, X] - \mathbb E[D|X]}.
#' @param S_split Integer. Number of sample splits for cross-fitting (default: 1).
#'
#' @return An object of class `IVDML`, containing:
#'   - `results_splits`: A list of S_split lists of cross-fitted residuals from the different sample splits.
#'   - `A`: The argument `A` of the function.
#'   - `ml_method`: The argument `ml_method` of the function.
#'   - `A_deterministic_X`: The argument `A_deterministic_X` of the function.
#'   - `iv_method`: The argument `iv_method` of the function.
#'   The treatment effect estimates, standard errors and confidence intervals can be calculated from the `IVDML` object using the functions TODO:ADD REFERENCE TO OTHER FUNCTIONS.
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
#' coef(fit, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2)
#'
#' @references ?Scheidegger et al. 2025, to be published ?
#'
#' @seealso Inference for a fitted `IVDML` object is done with the functions [coef.IVDML()], [se()], [standard_confint()] and [robust_confint()].
#'
#' @export
fit_IVDML <- function(Y, D, Z, X = NULL, A =  NULL, ml_method, ml_par = list(), A_deterministic_X = TRUE, K_dml = 5, iv_method = c("linearIV", "mlIV"), S_split = 1){
  stopifnot(is.numeric(Y), is.numeric(D), length(Y) == length(D))
  N <- length(Y)
  matrix_ZXA <- function(var, var_name) {
    if (!is.null(var)) {
      mat <- try(as.matrix(var), silent = TRUE)
      if (inherits(mat, "try-error")) {
        stop(paste(var_name, " cannot be converted to a matrix. Make sure that it is a vector, matrix or data frame.", sep = ""))
      }
      if (nrow(mat) != N) {
        stop(paste("The number of rows in ", var_name, " must match the length of Y and D.", sep = ""))
      }
      return(mat)
    }
    return(NULL)
  }
  Z <- matrix_ZXA(Z, "Z")
  X <- matrix_ZXA(X, "X")
  A <- matrix_ZXA(A, "A")
  if(!(ml_method %in% c("gam", "xgboost", "randomForest"))){
    stop(paste("ML method ", ml_method, " is not implemented.", sep = ""))
  }
  stopifnot(is.list(ml_par))
  stopifnot(is.logical(A_deterministic_X))
  if(!(A_deterministic_X)){
    if(is.null(A)){
      warning("If A is not a deterministic function of X, it should be specified, otherwise no heterogeneous treatment effect can be estimated.")
    } else {
      X <- cbind(X, A)
    }
  }
  if(!is.null(A)){
    if(!dim(A)[2]==1){
      warning("Heterogenous treatment effects can only be estimated with respect to univariate A.")
    }
  }
  if(!all(iv_method %in% c("linearIV", "mlIV", "mlIV_direct"))){
    stop(paste("Only IV methods linearIV, mlIV and mlIV_direct are implemented.", sep = ""))
  }
  stopifnot(is.numeric(K_dml), K_dml >= 2)
  stopifnot(is.numeric(S_split), S_split >= 1)

  tuned_ml_par <- tuning_helper(Y, D, Z, X, ml_method, ml_par, iv_method)
  one_rep <- function(i){residuals_helper(Y, D, Z, X, ml_method, tuned_ml_par, iv_method, K_dml)}
  results_splits <- lapply(1:S_split, one_rep)
  model <- list(results_splits = results_splits, A = A, ml_method = ml_method, A_deterministic_X = A_deterministic_X, iv_method = iv_method)
  class(model) <- "IVDML"

  return(model)
}
