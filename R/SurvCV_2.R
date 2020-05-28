
CV.Surv2 <- function(X0, Y0, status, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, folds=5, clv=NULL, r=5,
                    init=NULL, alpha.i=1, robust=TRUE, standardize=TRUE, ncores, verbo = FALSE)
{
  intercept = TRUE
  if(is.null(clv)){
    clv = intercept*1
  }else{
    clv = union(1, (clv+intercept))
  }

  n = nrow(X0); p.c = length(clv); p = ncol(X0)-p.c+intercept;
  if(standardize){
    V0 = apply(X0, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V0[V0==0]=1
    X1 = scale(X0, center = TRUE, scale = V0)
  }
  if(intercept) X1 = cbind(Intercept = rep(1, n), X1)
  Y1 = Y0

  out = KMweight(X1, Y1, status, robust)
  X = out$X + 10^-9
  Y = out$Y
  init = match.arg(init, choices = c("zero","cox","elnet"))

  if(is.null(lamb.1)){
    u=abs(t(X) %*% Y)
    if(!robust && penalty!="network") u=u/sqrt(n)
    LL = log(stats::quantile(u, 0.1)); UL = log(max(u))
    lamb.1 = rev(exp(seq(LL,UL,length.out = 35)))
  }

  if(is.null(lamb.2)){
    if(robust){
      lamb.2 = c(0.001, 0.01, 0.1, 1)
    }else{
      lamb.2 = c(0.01, 0.1, 1, 10)
    }
  }
  rs <- sample(c(1:n))
  # CVM = matrix(0, length(lamb.1), length(lamb.2));
  if(init == "cox"){
    b0 = initiation_cox(out$Xo, out$Yo, out$So)
  } else if(init == "elnet"){
    b0 = initiation(X, Y, alpha.i)
  } else{
    b0 = rep(0, (p+p.c))
  }
  a = Adjacency(X[,-clv,drop=FALSE])
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------

  Xc=X[, clv, drop = FALSE]; Xg = X[, -clv, drop = FALSE];
  CVM = SurvCV(Xc, Xg, Y, folds, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, robust, method, ncores)
  CVM = CVM/n
  mcvm = min(CVM)
  inds = which(CVM == mcvm, arr.ind=TRUE)
  lambda = lambda1 = lamb.1[inds[,1]]
  lambda2 = lamb.2[inds[,2]]
  if(length(lambda)>1) message("multiple optimal values(pairs) of lambda(s) are found.")
  rownames(CVM) = signif(lamb.1, digits = 3)
  if(penalty == "network"){
    lambda = cbind(lambda1, lambda2)
    colnames(CVM) = lamb.2
  }
  outlist = list(lambda=lambda, mcvm=mcvm, CVM=CVM, penalty=penalty)
  class(outlist) = "cv.surv"
  outlist
}
