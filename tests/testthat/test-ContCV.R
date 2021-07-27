test_that("ContCV", {
  n = 40; p = 5; Y = rep(0,n)
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5])); #Adjacency(X)
  Y = X[,4:7]%*%c(3,2,3,-2)

  X[,1] = 0; X[c(4,8),1] = c(1, -1);
  expect_error(CV.Cont(X, Y, "network", robust = TRUE, debugging = FALSE))
  out = CV.Cont(X, Y, "network", robust = TRUE, debugging = TRUE)
  expect_equal(out$penalty, "network")
  expect_equal(ncol(out$lambda), 2)


  expect_error(CV.Cont(X, Y, "mcp", lamb.2=0, robust = TRUE, debugging = FALSE))
  out = CV.Cont(X, Y, "mcp", lamb.2=0, robust = TRUE, debugging = TRUE)
  expect_equal(out$penalty, "mcp")
  expect_null(ncol(out$lambda))


  expect_error(CV.Cont(X, Y, "lasso", lamb.2=0, robust = TRUE, debugging = FALSE))
  out = CV.Cont(X, Y, "lasso", lamb.2=0, robust = TRUE, debugging = TRUE)
  expect_equal(out$penalty, "lasso")


  X[,1] = 0;
  expect_error(CV.Cont(X, Y, "network", robust = TRUE, debugging = TRUE), "standard deviation equal zero")
  out = CV.Cont(X, Y, "mcp", lamb.2=0, robust = TRUE, debugging = TRUE)
  expect_equal(out$penalty, "mcp")
  out = CV.Cont(X, Y, "lasso", lamb.2=0, robust = TRUE, debugging = TRUE)
  expect_equal(out$penalty, "lasso")


})
