
CV.Logit <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, folds=5, r=5, alpha=1,
                    init=NULL, alpha.i=1, standardize=TRUE, verbo = FALSE)
{
  if(is.null(lamb.1)){
    lamb.1 = switch (penalty,
                     "network" = lambda.n,
                     "mcp" = lambda.m,
                     "lasso" = lambda.l)
  }
  if(is.null(lamb.2)) lamb.2 = c(0.1, 1, 10)
  init = match.arg(init, choices = c("zero","elnet"))

  n = nrow(X); p = ncol(X);
  X = as.matrix(X); Y = as.matrix(Y)

  b0 = rep(0, p+1)
  rs <- sample(c(1:n))
  CVM = matrix(0, length(lamb.1), length(lamb.2))
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    index = c(1: ceiling(n/folds)) + (f-1)*ceiling(n/folds)
    test = rs[intersect(index, seq(1,n,1))]

    x = X[-test,]; y = Y[-test]
    x2 = X[test,]; y2 = Y[test]
    if(standardize){
      x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
      x2 = scale(x2, scale = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
    }
    if(penalty == "network") a = Adjacency(x) else a = as.matrix(0)
    x = cbind(rep(1,n-length(test)), x)
    x2 = cbind(rep(1,length(test)), x2)

    if(init == "elnet") b0 = initiation(x, y, alpha.i)

    CVM = CVM + LogitGrid(x, y, x2, y2, lamb.1, lamb.2, b0, r, a, p, alpha, method)

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
  class(outlist) = "cv.logit"
  outlist
}
