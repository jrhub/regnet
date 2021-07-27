test_that("LogitCV", {
  n = 40; p = 5; Y = rep(0,n)
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y0 = 1/(1+exp(-X[,4:7]%*%c(3,2,3,-2)))
  Y = rbinom(n,1,Y0); #hist(Y)

  out = CV.Logit(X, Y, "network", debugging=TRUE)
  out$lambda
  expect_length(unique(out$para$CVM2[out$para$inds]),1)
  expect_equal(min(out$para$CVM2[out$para$inds0]), out$para$CVM2[out$para$inds][1])
  expect_equal(dim(out$CVM), dim(out$para$CVM2))
  expect_equal(out$para$lamb.1[out$para$inds[,1]], as.numeric(out$lambda[,1]))
  expect_equal(out$para$lamb.2[out$para$inds[,2]], as.numeric(out$lambda[,2]))

  out = CV.Logit(X, Y, "mcp", lamb.2=0, debugging = TRUE)
  expect_equal(out$para$lamb.2,0)
  # par(mfrow=c(2,2)); for(i in 1:ncol(out$CVM)){plot(out$CVM[,i])}
  # expect_equal(min(out$para$CVM2[out$para$inds0]), out$para$CVM2[out$para$inds][1])
  expect_null(out$para$inds0)
  expect_equal(dim(out$CVM), dim(out$para$CVM2))
  expect_null(dim(out$lambda))

  out = CV.Logit(X, Y, "lasso", lamb.2=0, debugging = TRUE)
  expect_null(out$para$inds0)
  expect_null(dim(out$lambda))
  expect_equal(out$para$lamb.1[out$para$inds[,1]], out$lambda)

  X[,1] = 0;
  expect_error(CV.Logit(X, Y, "network"), "standard deviation equal zero")

})
