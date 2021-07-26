test_that("example_surv_robust", {
  skip_on_cran() # skip the test
  skip_on_travis()

  data(SurvExample)
  X = rgn.surv$X
  Y = rgn.surv$Y
  clv = c(1:5)
  fit = regnet(X, Y, "survival", "network", rgn.surv$lamb1, rgn.surv$lamb2, clv=clv, robust=TRUE)
  index = which(rgn.surv$beta[-(1:6)] != 0)  # [-(1:6)] removes the intercept and clinical variables that are not subject to selection.
  pos = which(fit$coeff[-(1:6)] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  expect_gt(tp,45)
  expect_lt(fp,5)
})


test_that("example_surv_nonrobust", {
  skip_on_cran() # skip the test
  skip_on_travis()

  data(SurvExample)
  X = rgn.surv$X
  Y = rgn.surv$Y
  clv = c(1:5)
  out = cv.regnet(X, Y, response="survival", penalty="network", clv=clv, robust=FALSE, verbo=FALSE)
  par(mfrow=c(2,2)); for(i in 1:ncol(out$CVM)){plot(out$CVM[,i])}
  fit = regnet(X, Y, "survival", "network", out$lambda[1,1], out$lambda[1,2], clv=clv, robust=FALSE)
  index = which(rgn.surv$beta[-(1:6)] != 0)  # [-(1:6)] removes the intercept and clinical variables that are not subject to selection.
  pos = which(fit$coeff[-(1:6)] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  expect_gt(tp,10)
  expect_lt(fp,100)
})


test_that("example_cont_robust", {
  skip_on_cran() # skip the test
  skip_on_travis()

  data(HeteroExample)
  X = rgn.htr$X
  Y = rgn.htr$Y
  fit = regnet(X, Y, "c", "n", 0.19, 1, r = 1.5, robust=TRUE)
  index = which(rgn.htr$beta[-1] != 0)   # [-1] removes the intercept
  pos = which(fit$coeff[-1] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  net = plot(fit)

  expect_gt(tp,25)
  expect_lt(fp,5)
  expect_length(igraph::V(net), length(pos))
})


test_that("example_cont_nonrobust", {
  skip_on_cran() # skip the test
  skip_on_travis()

  data(HeteroExample)
  X = rgn.htr$X
  Y = rgn.htr$Y
  robust= FALSE
  penalty="n"
  ncores = 2
  out = cv.regnet(X, Y, response="c", penalty=penalty, folds=5, r = 1.5, robust=robust, ncores=ncores, verbo=FALSE)
  par(mfrow=c(2,2)); for(i in 1:ncol(out$CVM)){plot(out$CVM[,i])}
  fit = regnet(X, Y, "c", penalty, out$lambda[1,1], out$lambda[1,2], r = 1.5, robust=robust)
  index = which(rgn.htr$beta[-1] != 0)   # [-1] removes the intercept
  pos = which(fit$coeff[-1] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  net = plot(fit)

  expect_gt(tp,10)
  expect_lt(fp,100)
  expect_length(igraph::V(net), length(pos))
})


test_that("example_cont_tcga", {
  skip_on_cran() # skip the test
  skip_on_travis()

  data(ContExample)
  X = rgn.tcga$X
  Y = rgn.tcga$Y
  clv = (1:2)
  fit = regnet(X, Y, "continuous", "network", rgn.tcga$lamb1, rgn.tcga$lamb2, clv =clv, alpha.i=0.5)
  expect_named(fit$coeff)

  nets = plot(fit,subnetworks = TRUE)
  expect_gte(length(nets),1)
  expect_gt(length(igraph::V(nets[[1]])$name),1)
})


test_that("example_logit", {
  skip_on_cran() # skip the test
  skip_on_travis()

  data(LogisticExample)
  X = rgn.logi$X
  Y = rgn.logi$Y
  penalty="n"
  out = cv.regnet(X, Y, response="binary", penalty=penalty, folds=5, r = 4.5)
  par(mfrow=c(2,2)); for(i in 1:ncol(out$CVM)){plot(out$CVM[,i])}
  fit = regnet(X, Y, "binary", penalty=penalty, out$lambda[1,1], out$lambda[1,2], r = 4.5)
  index = which(rgn.logi$beta[-1] != 0)   # [-1] removes the intercept
  pos = which(fit$coeff[-1] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  expect_gt(tp,10)
  expect_lt(fp,100)

  out = cv.regnet(X, Y, response="binary", penalty="m", folds=5, r = 4.5)
  par(mfrow=c(2,2)); for(i in 1:ncol(out$CVM)){plot(out$CVM[,i])}
  fit = regnet(X, Y, "binary", penalty="l", out$lambda[1], NULL, r = 4.5)
  index = which(rgn.logi$beta[-1] != 0)   # [-1] removes the intercept
  par(mfrow=c(2,2)); for(i in 1:ncol(out$CVM)){plot(out$CVM[,i])}
  pos = which(fit$coeff[-1] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  expect_gt(tp,1)
  expect_lt(fp,100)

  fit = regnet(X, Y, "binary", penalty=penalty, 0.055, 1, r = 4.5)
  index = which(rgn.logi$beta[-1] != 0)   # [-1] removes the intercept
  pos = which(fit$coeff[-1] != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
  expect_gt(tp,25)
  expect_lt(fp,5)

  net = plot(fit)
  a.tri = fit$Adj;a.tri[lower.tri(a.tri)]=0
  a.sort = sort(a.tri[which(a.tri!=0, arr.ind = TRUE)])
  expect_length(igraph::V(net), length(pos))
  expect_equal(sort(igraph::E(net)$weight),a.sort)
})
