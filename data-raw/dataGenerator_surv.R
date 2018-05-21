corName = "AR"
corr = 0.8
errName = "Cauchy25" # ("Normal", "Cauchy25", "T", "Cauchy15")
sim = 1

cat("corName:", corName, " corr:", corr, " errName:", errName, " sim:", sim, "\n\n")

library("glmnet")
library("MASS")
library("Rlab")
library("Rcpp")
library("RcppArmadillo")
# library("survAUC")
# library("regnet")

source("./data-raw/Network2.R")
source("./data-raw/correlation6.R")
#res <- try(sourceCpp("./WMRcpp/CD_WMR_HF.cpp", rebuild = T))
# sourceCpp(paste("./src/SurvCD.cpp", sep=""), rebuild = T)

n = 40; p = 60; p.c = 5+1;
n.sub = 5
#corr = 0.1

sig.c = SimCorStru.Clinic(p.c-1, 0.7);
if(corName == "AR"){
  sig = SimCorStru.AR(p, corr, n.sub);
  print("AR!")
}else if(corName == "BA.1"){
  sig = SimCorStru.BA.1(p, corr, n.sub);
}else if(corName == "BA.2"){
  sig = SimCorStru.BA.2(p, corr, n.sub);
}else if(corName == "BA.3"){
  sig = SimCorStru.BA.3(p, corr, n.sub);
}else{
  sig = SimCorStru.BLK(p, corr, n.sub);
}

r = r.m = 3
Pr = 1     # prob for bernulli
quant = 0.75


#------------------------------------------- validation / evaluation  ------------------------------------------

TruePos <- function(b, b.true){
  pos = setdiff(which(b != 0), c(1:p.c))
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
}

#--------------------------------------- simulate data set ------------------------
simulation <- function(error){
  x.g = mvrnorm(n,rep(0,p),sig)
  x.c = mvrnorm(n,rep(0,p.c-1),sig.c)
  # x = cbind(x.c, x.g)
  x = cbind(rep(1, nrow(x.c)), x.c, x.g)    # add intercept

  xb = c(x %*% b.true)
  LOG = min(xb + error)
  tt = exp(xb + error)
  tt[which(tt=="NaN"|tt=="Inf")] =  max(tt[which(tt!="NaN"&tt!="Inf")])
  cens = runif(n, min = min(tt), max=quantile(tt, quant))
  #cens = Inf

  y = tt*(tt<=cens)+cens*(tt>cens) + 10^-9
  status = 1 * (tt <= cens)
  list(Y = y, X = x, status = status)
}


sparsity <- function(n.sub, Pr){
  k = rep(0,p);  g = p/n.sub;
  index = sample(seq(1, g, 1), round(g*0.1/Pr), replace = FALSE)
  begining = (index-1)*n.sub + 1
  for(i in begining){
    k[i:(i+n.sub-1)] = rbern(n.sub, Pr)
  }
  k
}
#sum(sparsity(5))
#---------------------------------------------------------------------------------------------

b.true = runif(p, 0.2, 0.8) * sparsity(n.sub, Pr)   # 5 per cluster
# b.true = c(runif(p.c, 0.2, 0.8), b.true)
b.true = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, b.true)
index = setdiff(which(b.true != 0), c(1:p.c))
T = length(index); F = p - T
cat("b.true: ", b.true, "\n")

#--------------------------------------------------------------------------------------------
# ptm <- proc.time()
# sim = 1
# for(s in 1:sim){
  if(errName == "Normal"){
    error = rnorm(n, 0, 1);
  }else if(errName == "Cauchy25"){
    error = c(rnorm(n*0.75, 0, 1),rcauchy(n*0.25, 0, 1));
    print("25%Cauchy!")
  }else if(errName == "Cauchy15"){
    error = c(rnorm(n*0.85, 0, 1),rcauchy(n*0.15, 0, 1));
    print("15%Cauchy!")
  }else if(errName == "LogNormal"){
    error = c(rnorm(n*0.85, 0, 1),rlnorm(n*0.15, 0, 1));
    print("LogNormal!")
  }else if(errName == "T"){
    error = rt(n, 1);
  }
  error = sample(error, n, replace = FALSE);
  hist(error)
  train = simulation(error);
  # test = simulation();
  cat("censoring rate = ", 1-sum(train$status)/n, "\n")


  X0 = train$X[,-1]
  colnames(X0) = c(paste(" clin", seq = (1:(p.c-1)), sep=""), paste("gene", seq = (1:p), sep=""))
  Y0 = train$Y
  status = train$status
  # LL = log(2*10^(-5)); UL = log(0.02)
  # lam1.n = rev(exp(seq(LL,UL,length.out = 30)))
  # lamb.1=rev(exp(seq(1,45,1)/7 -9))/7 + 0.0004;
  lamb.1=NULL
  # lamb.2=c(0.001, 0.01, 0.1, 1)
  lamb.2=NULL
  r=3; alpha.i=1; folds=5;
  initiation="zero";  robust=TRUE; standardize=TRUE
  clv = c(1:5)
  Y = cbind(time = Y0, status = status)
  initiation = NULL
  penalty = "network"
  robust = FALSE # TRUE FALSE

  # X0=rgn.surv$X
  # Y=rgn.surv$Y
  # b.true=rgn.surv$beta

  ptm <- proc.time()
  # result = CV.NetSurv(X0, Y0, status, lamb.1=lamb.1, lamb.2=lamb.2, clv=clv, r=3, alpha.i=1, folds=5, init=NULL,
  #            robust=robust, standardize=TRUE, verbo = TRUE)
  result = cv.regnet(X0, Y, response=c("survival"), penalty=penalty, lamb.1=lamb.1, lamb.2=lamb.2, folds=5, r=r, clv=clv,
                     initiation=initiation, alpha.i=1, robust=robust, standardize=TRUE, verbo = TRUE)
  print(proc.time() - ptm)
  result$lambda
  lamb.1=result$lambda[1];  lamb.2=result$lambda[2]

  # bn = NetSurv(X0, Y0, status, lamb.1=lamb.1, lamb.2=lamb.2, clv=clv, r=r, alpha.i=1,init=initiation,
  #         robust=robust, standardize=TRUE)

  bn = regnet(X0, Y, response=c("survival"), penalty=penalty, lamb.1=lamb.1, lamb.2=lamb.2,  r=r, clv=clv,init=initiation, alpha.i=1,
              robust=robust, standardize=TRUE)
  TF = TruePos(bn, b.true)
  cat("TP: ", TF$tp, " | FP: ", TF$fp, " | ", sum(bn!=0), " | ", bn[1:10],"\n")

  ####################################################################################

  rgn.surv = list(X=X0, Y=Y, beta=b.true)

  ##################### L2 Network ##############################################
  penalty = "network"

  ptm <- proc.time()
  result = CV.NetLogistic(rgn.logi$X, rgn.logi$Y, r = 4.5)
  print(proc.time() - ptm)
  result$lambda

  ptm <- proc.time()
  res = cv.regnet(rgn.logi$X, rgn.logi$Y, response=c("binary"), penalty=penalty, folds=5, r = 4.5)
  print(proc.time() - ptm)
  res$lambda

  # b = regnet(rgn.logi$X, rgn.logi$Y, response=c("binary"), penalty, lamb.1=res$lambda[1], r = 4.5)
  b = regnet(rgn.logi$X, rgn.logi$Y, response=c("binary"), penalty=penalty, res$lambda[1,1], res$lambda[1,2], r = 4.5)
  index = which(rgn.logi$beta != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)


  b = McpLogistic(rgn.logi$X, rgn.logi$Y, result$lambda[1,1], result$lambda[1,2])
  # b = McpLogistic(rgn.logi$X, rgn.logi$Y, result$lambda[1])
  index = which(rgn.logi$beta != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)

  b = McpLogistic(rgn.logi$X, rgn.logi$Y, res$lambda[1,1], res$lambda[1,2])
  # b = McpLogistic(rgn.logi$X, rgn.logi$Y, res$lambda[1])
  index = which(rgn.logi$beta != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
