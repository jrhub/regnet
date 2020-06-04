test_that("cv.regnet_data_format_cont", {
  n = 27; p = 3;
  X = scale(matrix(rnorm(n*p,0,5), n, p), scale=TRUE)
  Y= 5 + X%*%c(0,3,0) + rnorm(n)
  expect_error(cv.regnet(X[,1], Y, "c", "n"), "too less variables for network penalty")
  expect_error(cv.regnet(X[1:4,], Y[1:4], "c", "n"), "sample size too small")
  expect_error(cv.regnet(X, Y[1:4], "c", "n"), "length of Y does not match")
  expect_error(cv.regnet(X, Y, "c", "n", folds = -5), "incorrect value of folds")
  expect_error(cv.regnet(X, Y, "c", "n", alpha.i = 2), "alpha.i should be between 0 and 1")
  expect_error(cv.regnet(X, Y, "c", "n", ncore = 0), "incorrect value of ncores")

  out = cv.regnet(X, Y, "c", "n", lamb.2 = 1)
  expect_equal(ncol(out$lambda), 2)
  expect_equal(colnames(out$CVM), "1")
  expect_equal(out$penalty, "network")

  out = cv.regnet(X, Y, "c", "n")
  expect_equal(colnames(out$CVM), c("0.1", "1", "10"))

  out = cv.regnet(X, Y, "c", "m")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "mcp")

  out = cv.regnet(X, Y, "c", "l")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "lasso")
  # fit = regnet(X, Y, "c", "n", out$lambda[1], lamb.2=1); fit$coeff
})

test_that("cv.regnet_data_format_surv", {
  n = 27; p = 3;
  X = scale(matrix(rnorm(n*p,0,1), n, p), scale=TRUE)
  Y0= exp(2 + X%*%c(0,3,0) + rnorm(n)); Y1 = sample(rep(c(0,1,1,1),n/2),n)
  Y = data.frame(time=(Y0+Y0*(Y1-1)*0.2), status=Y1)

  expect_error(cv.regnet(X, Y0, "s", "n"), "Y should be a two-column matrix")
  expect_error(cv.regnet(X, cbind(Y[,1], Y[,2]), "s", "n"), "columns named 'time' and 'status'")
  expect_error(cv.regnet(X, data.frame(time=Y[,1]-10, status=Y[,2]), "s", "n"), "survival times need to be positive")
  expect_error(cv.regnet(X, data.frame(time=Y[,1], status=Y[,2]*2), "s", "n"), "binary variable of 1 and 0")
  expect_error(cv.regnet(X, Y[1:4,], "s", "n"), "the number of rows of Y")

  out = cv.regnet(X, Y, "s", "n", lamb.2 = 1)
  expect_equal(ncol(out$lambda), 2)
  expect_equal(colnames(out$CVM), "1")
  expect_equal(out$penalty, "network")

  out = cv.regnet(X, Y, "s", "m")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "mcp")
  # fit = regnet(X, Y, "s", "n", out$lambda[1], lamb.2=1); fit$coeff
})

test_that("cv.regnet_data_format_logit", {
  n = 52; p = 3; Y = rep(0,n)
  X = scale(matrix(rnorm(n*p,0,5), n, p), scale=TRUE)
  Y0 = 1/(1+exp(-X%*%c(0,5,0)))
  Y = rbinom(n,1,Y0);

  expect_error(cv.regnet(X, Y0, NULL, "n"), "Y must be a binary variable")
  expect_error(cv.regnet(X, Y[1:4], "b", "n"), "length of Y does not match")
  expect_message(cv.regnet(X, Y, "b", "n", robust=TRUE), "robust methods are not available")

  out = cv.regnet(X, Y, "b", "n", lamb.2 = 1)
  expect_equal(ncol(out$lambda), 2)
  expect_equal(colnames(out$CVM), "1")
  expect_equal(out$penalty, "network")

  out = cv.regnet(X, Y, "b", "m")
  expect_equal(ncol(out$CVM), 1)
  expect_equal(out$penalty, "mcp")
  # fit = regnet(X, Y, "b", "m", out$lambda[1]); fit$coeff

})
