#' @useDynLib regnet
#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for network-based logistic regression.
#'
#' This function does k-fold cross-validation for the network-based logistic regression and returns
#' a pair of lambda1 and lambda2.
#'
#' @param X a matrix of predictors.
#' @param Y a vector of the binary response.
#' @param lamb.1 a user-supplied sequence of lambda1 values, which serves as a tuning parameter to impose sparsity.
#' If it is left as NULL, a default sequence will be used.
#' @param lamb.2 a user-supplied sequence of lambda2 values, which serves as a tuning parameter to control the smoothness
#' among coefficient profiles. If it is left as NULL, a default sequence, c(0.1, 1, 10), will be used.
#' @param r the regularization parameter in MCP; default is 5.
#' @param alpha.i by default, the program uses the lasso for choosing initial values of
#' the coefficient vector. alpha.i is the Elastic-Net mixing parameter, with \eqn{0 \le alpha.i \le 1}. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 the ridge penalty. If alpha.i is assigned to be -1, the program will use zeroes
#' as initial coefficients.
#' @param folds the number of folds for cross-validation; default is 5.
#' @return a list with components:
#' \item{lambda}{the optimal pair of lambda1 and lambda2.}
#' \item{mcr}{the misclassification rate of the optimal pair of lambda1 and lambda2.}
#' \item{MCR}{a matrix of the misclassification rates for all pairs of lambdas tested.}
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. BMC Genetics, 18(1):44.
#'
#' @seealso \code{\link{NetLogistic}}
#'
#' @export

CV.NetLogistic <- function(X, Y, lamb.1=NULL, lamb.2=NULL, r=5, alpha.i=1, folds=5, verbo = FALSE){

  if(is.null(lamb.1)) lamb.1 = lambda.n
  if(is.null(lamb.2)) lamb.2 = c(0.1, 1, 10)
  n = nrow(X); p = ncol(X);
  X = as.matrix(X); Y = as.matrix(Y)

  b0 = rep(0, p+1)
  rs <- sample(c(1:n))
  tMSE = matrix(0, length(lamb.2), length(lamb.1))
  #---------------------------------------------- Main Loop -----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    index = c(1: ceiling(n/folds)) + (f-1)*ceiling(n/folds)
    test = rs[intersect(index, seq(1,n,1))]

    x = X[-test,]; y = Y[-test]
    x2 = X[test,]; y2 = Y[test]
    x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
    x2 = scale(x2, scale = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))

    a = Adjacency(x)
    x = cbind(rep(1,n-length(test)), x)
    x2 = cbind(rep(1,length(test)), x2)

    if(alpha.i != -1) b0 = initiation(x, y, alpha.i)
    n.x = nrow(x)

    for(j in 1:length(lamb.2)){
      for(i in 1:length(lamb.1)){ # Network
        # b = run.net(x, y, lamb.1[i], lamb.2[j], b0, r, a, n.x, p)
        b = RunNet(x, y, lamb.1[i], lamb.2[j], b0, r, a, n.x, p)
        tMSE[j,i] = tMSE[j,i] + validation(b, x2, y2, n)
      }
    }
    # tMSE = RunNet_Grid(x, y, x2, y2, lamb.1, lamb.2, b0, r, a, n.x, p)
  }
  mcr = min(tMSE)
  inds = which(tMSE == mcr, arr.ind=TRUE)
  lambda1 = lamb.1[inds[,2]]
  lambda2 = lamb.2[inds[,1]]
  lambda = cbind(lambda1, lambda2)
  colnames(tMSE) = signif(lamb.1, digits = 3)
  rownames(tMSE) = lamb.2
  return(list(lambda=lambda, mcr=mcr, MCR=tMSE))
}

#' Network-based logistic regression for given lambda1 and lambda2 pair.
#'
#' This function makes predictions for network-based logistic regression for a given pair of lambda1 and lambda2 values.
#' Typical usage is to have the CV.NetLogistic function compute the optimal lambdas, then provide them to the
#' NetLogistic function.
#'
#' @param X a matrix of predictors.
#' @param Y a vector of the binary response.
#' @param lamb.1 the tuning parameter (lambda1) that imposes sparsity.
#' @param lamb.2 the tuning parameter (lambda2) that controls the smoothness among coefficient profiles.
#' @param alpha.i by default, the program uses Elastic-Net for choosing initial values of
#' the coefficient vector. alpha.i is the Elastic-Net mixing parameter, with \eqn{0 \le alpha.i \le 1}. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 is the ridge penalty. If alpha.i is assigned to be -1, the program will use zeroes
#' as initial coefficients.
#' @param r the regularization parameter in MCP.
#' @param folds the number of folds for cross-validation.
#' @return the estimated coefficients vector.
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. BMC Genetics.
#'
#' @seealso \code{\link{CV.NetLogistic}}
#'
#' @examples
#' b = NetLogistic(regnet$X, regnet$Y, 0.05, 1)
#' regnet$beta
#' @export

NetLogistic <- function(X, Y, lamb.1, lamb.2, alpha.i=1, r=5, folds=5){
  n = nrow(X); p = ncol(X);
  x = as.matrix(X); y = as.matrix(Y)
  x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
  a = Adjacency(x)
  x = cbind(rep(1,n), x)
  if(alpha.i != -1) b0 = initiation(x, y, alpha.i)
  else b0 = rep(0, p+1)
  # b = run.net(x, y, lamb.1, lamb.2, b0, r, a, n, p)
  b = RunNet(x, y, lamb.1, lamb.2, b0, r, a, n, p)
}

run.net <- function(x, y, lam1, lam2, b, r, a, n, p){
  count = 0
  while(count < 20){
    b.new = Network(x, y, lam1, lam2, b, r, a, n, p)
    dif = sum(abs(b - b.new))/n
    #cat("L1 Diff: ", dif, ", for lam1: ", lam1, "\n")
    if( dif < 0.001) break
    else{
      b = b.new; count = count +1
    }
  }
  b.new
}

Network <- function(x, y, lam1, lam2, b, r, a, n, p){
  y0 = x %*% b
  pi = 1/(1+exp(-y0))
  t = (y - pi)*4
  for( k in 1: length(b)){
    b.old = b[k]
    l = t(x[,k]) %*% t/n + b[k]
    if(k == 1) b[k] = l          # intercept
    else{
      m = min(k,p)
      z = l*0.25 + lam2 * (t(a[(k-1),m:p]) %*% b[(m+1):(p+1)])
      u = 0.25 + lam2 * sum(abs(a[(k-1),m:p]))

      if(abs(z) > (r*lam1*u)) b[k] = z / u
      else b[k] = Soft(z, lam1)/ (u - 1/r)
    }
    t = t - x[,k] * (b[k] - b.old)
  }
  b
}

