test_that("tune_ml stops if ml_method is not avaliable", {
  expect_error(tune_ml(Y = rnorm(10), X = rnorm(10), ml_method = "gaaaaaaaaaaam", ml_par = list()))
})



### TODO: where do we handle X = NULL? ###
