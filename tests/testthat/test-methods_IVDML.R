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

test_that("se.IVDML works (and does not give NULL)", {
  set.seed(1)
  Y <- rnorm(10)
  D <- rnorm(10)
  Z <- rnorm(10)
  X <- rnorm(10)
  fitted_object <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam")
  fitted_se <- se(fitted_object, iv_method = "mlIV")
  expect_equal(is.null(fitted_se), FALSE)
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
