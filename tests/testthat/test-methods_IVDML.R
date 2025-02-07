test_that("coef.IVDML works (and does not give NULL)", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  fitted_coef <- coef(fitted_object, iv_method = "mlIV")
  expect_equal(is.null(fitted_coef), FALSE)
})

test_that("coef.IVDML works with heterogeneous treatment effect", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  A <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_no_error(coef(fitted_object, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 1))
})

test_that("coef.IVDML works with heterogeneous treatment effect, when A is only specified in the fit_IVDML function", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  A <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = "gam")
  expect_no_error(coef(fitted_object, iv_method = "mlIV", a = 0, kernel_name = "boxcar", bandwidth = 1))
})

test_that("coef.IVDML raises error when A is neither specified in fit_IVDML nor in coef", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  A <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_error(coef(fitted_object, iv_method = "mlIV", a = 0, kernel_name = "boxcar", bandwidth = 1))
})

test_that("se works (and does not give NULL)", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  fitted_se <- se(fitted_object, iv_method = "mlIV")
  expect_equal(is.null(fitted_se), FALSE)
})

test_that("se works with heterogeneous treatment effect", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  A <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_no_error(se(fitted_object, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 1))
})

test_that("se works with heterogeneous treatment effect, when A is only specified in the fit_IVDML function", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  A <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = "gam")
  expect_no_error(se(fitted_object, iv_method = "mlIV", a = 0, kernel_name = "boxcar", bandwidth = 1))
})

test_that("se raises error when A is neither specified in fit_IVDML nor in coef", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  A <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_error(se(fitted_object, iv_method = "mlIV", a = 0, kernel_name = "boxcar", bandwidth = 1))
})


test_that("standard_confint works (and does not give NULL)", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  fitted_confint <- standard_confint(fitted_object, iv_method = "mlIV")
  expect_equal(is.null(fitted_confint), FALSE)
})

test_that("DML_agg and MMB_agg give the same result when S_split == 1", {
  set.seed(1)
  Z <- rnorm(100)
  X <- Z + rnorm(100)
  H <- rnorm(100)
  D <- Z^2 + sin(X) + H + rnorm(100)
  A <- X
  Y <- tanh(A) * D + cos(X) - H + rnorm(100)
  fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = "gam")
  p_DML <- robust_p_value_aggregated(fit, candidate_value = 0, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2, agg_method = "DML_agg")
  p_MMB <- robust_p_value_aggregated(fit, candidate_value = 0, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2, agg_method = "MMB_agg")
  expect_equal(p_DML, p_MMB)
})

test_that("DML_agg and MMB_agg work without error when S_split > 1", {
  set.seed(1)
  Z <- rnorm(100)
  X <- Z + rnorm(100)
  H <- rnorm(100)
  D <- Z^2 + sin(X) + H + rnorm(100)
  A <- X
  Y <- tanh(A) * D + cos(X) - H + rnorm(100)
  fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = "gam", S_split = 3)
  expect_no_error(p_DML <- robust_p_value_aggregated(fit, candidate_value = 0, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2, agg_method = "DML_agg"))
  expect_no_error(p_MMB <- robust_p_value_aggregated(fit, candidate_value = 0, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2, agg_method = "MMB_agg"))
})


test_that("robust_confint yields warning when no CI_range is specified", {
  set.seed(1)
  Z <- rnorm(100)
  X <- Z + rnorm(100)
  H <- rnorm(100)
  D <- Z^2 + sin(X) + H + rnorm(100)
  A <- X
  Y <- tanh(A) * D + cos(X) - H + rnorm(100)
  fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_warning(robust_confint(fit, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2))
})


test_that("robust_confint yields no warning when CI_range is specified", {
  set.seed(1)
  Z <- rnorm(100)
  X <- Z + rnorm(100)
  H <- rnorm(100)
  D <- Z^2 + sin(X) + H + rnorm(100)
  A <- X
  Y <- tanh(A) * D + cos(X) - H + rnorm(100)
  fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_no_warning(robust_confint(fit, iv_method = "mlIV", a = 0, A = A, kernel_name = "boxcar", bandwidth = 0.2, CI_range = c(-10, 10)))
})

test_that("print works", {
  set.seed(1)
  Z <- rnorm(100)
  X <- Z + rnorm(100)
  H <- rnorm(100)
  D <- Z^2 + sin(X) + H + rnorm(100)
  A <- X
  Y <- tanh(A) * D + cos(X) - H + rnorm(100)
  fit <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  expect_no_message(print(fit))
})


