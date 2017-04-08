#' k-folds cross-validation for MCP logistic regression.
#'
#' This function dose k-fold cross-validation for the MCP logistic regression and returns
#' a value of lambda.
#'
#' @param X a matrix of predictors.
#' @param Y a vector of the binary response.
#' @param lambda a user-supplied sequence of lambda. Tuning parameter lambda imposes sparsity.
#' If it is left as NULL, a default sequence will be used.
#' @param r the regularization parameter in MCP.
#' @param alpha.i by default, the program use the lasso for choosing initial values of
#' the coefficient vector. alpha.i is the elastic-net mixing parameter, with \eqn{0 \le alpha.i \le 1}. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 the ridge penalty. If assign alpha.i to be -1, program will use zero
#' as initial coefficients.
#' @param folds the number of folds for cross-validation.
#' @return a list with components:
#' \item{lambda}{the optimal lambda.}
#' \item{mcr}{the misclassification rate of the optimal lambda.}
#' \item{MCR}{a matrix of the misclassification rates for all the values of lambda tested.}
#'
#' @seealso \code{\link{McpLogistic}}
#'
#' @export
CV.McpLogistic <- function(X, Y, lambda=NULL, r=5, alpha.i=1, folds=5){

  if(is.null(lambda)) lambda = lambda.m
  n = nrow(X); p = ncol(X);
  X = as.matrix(X); Y = as.matrix(Y)

  b0 = rep(0, p+1)
  rs <- sample(c(1:n))
  tMSE = matrix(0, 1, length(lambda))
  #------------------------------------------ Main Loop ----------------------------------------------
  for(f in 1:folds){
      cat("CrossValidation: ",f, "/", folds, "\n")
      index = c(1: ceiling(n/folds)) + (f-1)*ceiling(n/folds)
      test = rs[intersect(index, seq(1,n,1))]

      x = X[-test,]; y = Y[-test]
      x2 = X[test,]; y2 = Y[test]
      x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))        #normalize the data
      x2 = scale(x2, scale = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))        #normalize the data

      x = cbind(rep(1,n-length(test)), x)
      x2 = cbind(rep(1,length(test)), x2)

      if(alpha.i != -1) b0 = initiation(x, y, alpha.i)
      n.x = nrow(x)

      for(i in 1:length(lambda)){# MCP
        b = run.mcp(x, y, lambda[i], b0, r, n.x, p)
        tMSE[1,i] = tMSE[1,i] + validation(b, x2, y2, n)
      }
        #graph(lambda, mse.mcp, "Log MCP ",f, "/", "folds")
  }

  mcr = min(tMSE)
  inds = which(tMSE == mcr, arr.ind=FALSE)
  lambda.opt = lambda[inds]
  colnames(tMSE) = signif(lambda, digits = 3)

  return(list(lambda=lambda.opt, mcr=mcr, MCR=tMSE))
}

#' MCP logistic regression for a given lambda.
#'
#' This function makes predictions for MCP logistic for a given value of lambda1.
#' Typical usage is to have the CV.MCPLogistic function compute the optimal lambda, then provide it to
#' the McpLogistic function.
#'
#' @param X a matrix of predictors.
#' @param Y a vector of the binary response.
#' @param lambda the tunning parameter lambda imposes sparsity.
#' @param r the regularization parameter in MCP.
#' @param alpha.i by default, the program use the lasso for choosing initial values of
#' the coefficient vector. alpha.i is the elastic-net mixing parameter, with \eqn{0 \le alpha.i \le 1}. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 the ridge penalty. If assign alpha.i to be -1, program will use zero
#' as initial coefficients.
#' @param folds the number of folds for cross-validation.
#' @return the estimated coefficients vector.
#'
#' @seealso \code{\link{CV.McpLogistic}}
#'
#' @examples
#' b = McpLogistic(regnet$X, regnet$Y, 0.075)
#' regnet$beta # the ture coefficient
#' @export
McpLogistic <- function(X, Y, lambda, r=5, alpha.i=1, folds=5){
  n = nrow(X); p = ncol(X);
  x = as.matrix(X); y = as.matrix(Y)
  x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
  x = cbind(rep(1,n), x)
  if(alpha.i != -1) b0 = initiation(x, y, alpha.i)
  else b0 = rep(0, p+1)
  b = run.mcp(x, y, lambda, b0, r, n, p)
}

run.mcp <- function(x, y, lambda, b, r, n, p){
  count = 0
  while(count < 20){
    b.new = MCP(x, y, lambda, b, r, n, p)
    dif = sum(abs(b - b.new))/n
    #cat("L1 Diff: ", dif, ", for lambda: ", lambda, "\n")
    if( dif < 0.0005) break
    else{
      b = b.new; count = count +1
    }
  }
  b.new
}

MCP <- function(x, y, lambda, b, r, n, p){
  y0 = x %*% b
  pi = 1/(1+exp(-y0))
  t = (y - pi)*4
  for( k in 1: length(b)){
    b.old = b[k]
    z = ( t(x[,k]) %*% t/n + b[k] )*0.25
    if(k == 1) b[k] = z /0.25
    else if(abs(z) > (r*lambda*0.25)) b[k] = z /0.25
    else b[k] = Soft(z, lambda)/(0.25 - 1/r)
    t = t - x[,k] * (b[k] - b.old)
  }
  b
}
