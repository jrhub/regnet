test_that("LogitCV", {
  n = 40; p = 5; Y = rep(0,n)
  x = matrix(rnorm(n*p,0,5), n, p)
  X = scale(data.frame(x,X6=x[,4]+x[,5]*0.5,X7=x[,4]*0.2-x[,5]), scale=TRUE); #Adjacency(X)
  Y0 = 1/(1+exp(-X[,4:7]%*%c(3,2,3,-2)))
  Y = rbinom(n,1,Y0); #hist(Y)

  out = cv.regnet(X, Y, NULL, "n")
  expect_output(print(out), "CV error")
})
