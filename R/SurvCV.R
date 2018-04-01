
CV.Surv <- function(X0, Y0, status, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, clv=NULL, folds=5, r=5,
                    init=NULL, alpha.i=1, robust=TRUE, standardize=TRUE, verbo = FALSE)
{
  intercept = TRUE
  status = as.numeric(status)
  if(is.null(clv)){
    clv = intercept*1
  }else{
    clv = union(1, (clv+intercept))
  }

  n = nrow(X0); p.c = length(clv); p = ncol(X0)-p.c+intercept;
  if(standardize) X1 = scale(X0, center = FALSE, scale = apply(X0, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
  if(intercept) X1 = cbind(Intercept = rep(1, n), X1)
  Y1 = Y0

  out = KMweight(X1, Y1, status, robust)
  X = out$X + 10^-9
  Y = out$Y
  init = match.arg(init, choices = c("zero","cox","elnet"))

  if(is.null(lamb.1)){
    u=abs(t(X) %*% Y)
    LL = log(stats::quantile(u, 0.1)); UL = log(max(u))
    lamb.1 = rev(exp(seq(LL,UL,length.out = 30)))
  }
  # if(penalty != "network") lamb.2 = 0
  if(is.null(lamb.2)){
    if(robust){
      lamb.2 = c(0.001, 0.01, 0.1, 1)
    }else{
      lamb.2 = c(0.01, 0.1, 1, 10)
    }
  }
  rs <- sample(c(1:n))
  CVM = matrix(0, length(lamb.1), length(lamb.2));
  #---------------------------------------------- Main Loop -----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    index = c(1: ceiling(n/folds)) + (f-1)*ceiling(n/folds)
    test = rs[intersect(index, seq(1,n,1))]

    x = X[-test,]; y = Y[-test];
    x.c=x[, clv, drop = FALSE]; x.g = x[, -clv, drop = FALSE];

    x2 = X[test,]; x2 = cbind(x2[,clv], x2[,-clv])
    y2 = Y[test]

    if(init == "cox"){
      b0 = initiation_cox(out$Xo[-test,], out$Yo[-test], out$So[-test])
      # cat("#b0: ", sum(b0 != 0), "TP FP: ", TruePos(b0)$tp, "|", TruePos(b0)$fp, " | b0: ", b0[1:100], "\n")
    } else if(init == "elnet"){
      b0 = initiation(x, y, alpha.i)
      # cat("#b0: ", sum(b0 != 0), "TP FP: ", TruePos(b0)$tp, "|", TruePos(b0)$fp, " | b0: ", b0[1:100], "\n")
    } else{
      b0 = rep(0, (p+p.c))
    }

    if(penalty == "network"){
      a = Adjacency(x.g)
      CVM = CVM + NetGrid(x.c, x.g, y, x2, y2, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, robust)
    }else if(penalty == "mcp"){
      CVM = CVM + MCPGrid(x.c, x.g, y, x2, y2, lamb.1, b0[clv], b0[-clv], r, p, p.c, robust)
    }else{
      CVM = CVM + LassoGrid(x.c, x.g, y, x2, y2, lamb.1, b0[clv], b0[-clv], p, p.c, robust)
    }

  }
  CVM = CVM/n
  mcvm = min(CVM)
  inds = which(CVM == mcvm, arr.ind=TRUE)
  lambda1 = lamb.1[inds[,1]]
  lambda2 = lamb.2[inds[,2]]
  lambda = lambda1
  if(penalty == "network") lambda = cbind(lambda1, lambda2)
  rownames(CVM) = signif(lamb.1, digits = 3)
  if(penalty == "network") colnames(CVM) = lamb.2
  outlist = list(lambda=lambda, mcvm=mcvm, CVM=CVM, penalty=penalty)
  class(outlist) = "cv.surv"
  outlist
}
