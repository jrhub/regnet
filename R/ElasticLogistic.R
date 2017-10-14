#' k-folds cross-validation for Elastic-Net logistic regression.
#'
#' This function does k-fold cross-validation for the Elastic-Net logistic regression and returns
#' the optimal value of lambda.
#'
#' @param X a matrix of predictors.
#' @param Y a vector of the binary response.
#' @param lambda a user-supplied sequence of lambda values, which serves as a tuning parameter to impose sparsity.
#' If it is left as NULL, a default sequence will be used.
#' @param alpha the Elastic-Net mixing parameter, with \eqn{0 \le \alpha \le 1}. alpha=1 corresponds to the lasso penalty,
#' and alpha=0 corresponds to the ridge penalty.
#' @param alpha.i by default, the program uses the lasso penalty for choosing initial values of
#' the coefficient vector. alpha.i is the Elastic-Net mixing parameter, with \eqn{0 \le alpha.i \le 1}. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 the ridge penalty. If alpha.i is assigned to be -1, the program will use zeroes
#' as initial coefficients.
#' @param folds the number of folds for cross-validation.
#' @param verbo output progress to the console.
#' @return a list with components:
#' \item{lambda}{the optimal lambda.}
#' \item{mcr}{the misclassification rate of the optimal lambda.}
#' \item{MCR}{a matrix of the misclassification rates for all the values of lambda tested.}
#'
#' @seealso \code{\link{ElasLogistic}}
#'
#' @export
CV.ElasLogistic <- function(X, Y, lambda=NULL, alpha=0.5, alpha.i=1, folds=5, verbo = FALSE){

  if(is.null(lambda)) lambda = lambda.e
  n = nrow(X); p = ncol(X);
  X = as.matrix(X); Y = as.matrix(Y)

  b0 = rep(0, p+1)
  rs <- sample(c(1:n))
  tMSE = matrix(0, 1, length(lambda))
  #-------------------------------------------- Main Loop ----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    index = c(1: ceiling(n/folds)) + (f-1)*ceiling(n/folds)
    test = rs[intersect(index, seq(1,n,1))]

    x = X[-test,]; y = Y[-test]
    x2 = X[test,]; y2 = Y[test]
    x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
    x2 = scale(x2, scale = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))

    x = cbind(rep(1,n-length(test)), x)
    x2 = cbind(rep(1,length(test)), x2)
    if(alpha.i != -1) b0 = initiation(x, y, alpha.i)
    n.x = nrow(x)
    for(i in 1:length(lambda)){ #Elastic net
      # b = run.elastic(x, y, lambda[i], b0, alpha, n.x, p)
      b = RunElastic(x, y, lambda[i], b0, alpha, n.x, p)
      tMSE[1,i] = tMSE[1,i] + validation(b, x2, y2, n)
    }
  }

  mcr = min(tMSE)
  inds = which(tMSE == mcr, arr.ind=FALSE)
  lambda.opt = lambda[inds]
  colnames(tMSE) = signif(lambda, digits = 3)

  return(list(lambda=lambda.opt, mcr=mcr, MCR=tMSE))
}

#' Elastic-Net logistic regression for a given lambda.
#'
#' This function makes predictions for Elastic-Net logistic regression for a given value of lambda.
#' Typical usage is to have the CV.ElasLogistic function compute the optimal lambda, then provide it to
#' the ElasLogistic function.
#'
#' @param X a matrix of predictors.
#' @param Y a vector of the binary response.
#' @param lambda the tuning parameter that imposes sparsity.
#' @param alpha the Elastic-Net mixing parameter, with \eqn{0 \le \alpha \le 1}. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
#' @param alpha.i by default, the program use the lasso for choosing initial values of
#' the coefficient vector. alpha.i is the Elastic-Net mixing parameter, with \eqn{0 \le alpha.i \le 1}. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 the ridge penalty. If alpha.i is assigned as -1, the program will use zeroes
#' as initial coefficients.
#' @param folds the number of folds for cross-validation.
#' @return the estimated coefficients vector.
#'
#' @seealso \code{\link{CV.ElasLogistic}}
#'
#' @examples
#' b = ElasLogistic(regnet$X, regnet$Y, 0.04)
#' inds = which(regnet$beta != 0)
#' sel = which(b != 0)
#' tp = length(intersect(inds, sel))
#' fp = length(sel) - tp
#' list(tp=tp, fp=fp)
#' @export
ElasLogistic <- function(X, Y, lambda, alpha=0.5, alpha.i=1, folds=5){
  n = nrow(X); p = ncol(X);
  x = as.matrix(X); y = as.matrix(Y)
  x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
  x = cbind(rep(1,n), x)
  b0 = initiation(x, y, alpha.i)
  # b = run.elastic(x, y, lambda, b0, alpha, n, p)
  b = RunElastic(x, y, lambda, b0, alpha, n, p)
}

run.elastic <- function(x, y, lambda, b, alpha, n, p){
  count = 0
  while(count < 20){
    b.new = ElasticNet(x, y, lambda, b, alpha, n, p)
    dif = sum(abs(b - b.new))/n
    #cat("L1 Diff: ", dif, ", for lambda: ", lambda, "\n")
    if( dif < 0.0005) break
    else{
      b = b.new; count = count +1
    }
  }
  b.new
}

ElasticNet <- function(x, y, lambda, b, alpha, n, p){
  y0 = x %*% b
  pi = 1/(1+exp(-y0))
  r = (y - pi)*4
  for( k in 1: length(b)){
    b.old = b[k]
    z = ( t(x[,k]) %*% r /n + b[k] )*0.25
    if(k == 1) b[k] = z*4
    else b[k] = Soft(z, lambda*alpha)/(lambda*(1-alpha)+0.25)
    r = r - x[,k] * (b[k] - b.old)
  }
  b
}
