#' @useDynLib regnet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for regnet
#'
#' This function does k-fold cross-validation for regnet and returns the optimal value(s) of lambda.
#'
#' @param X X matrix as in regnet.
#' @param Y response Y as in regnet.
#' @param response response type. regnet now supports two types of response: "binary" and "survival".
#' @param penalty penalty type. regnet provides three choices of penalty: Network, MCP and Lasso.
#' @param lamb.1 a user-supplied sequence of \eqn{\lambda 1} values, which serves as a tuning parameter to impose sparsity.
#' If it is left as NULL, regnet will compute its own sequence.
#' @param lamb.2 a user-supplied sequence of \eqn{\lambda 2} values for network method. \eqn{\lambda 2} control the smoothness
#' among coefficient profiles. If it is left as NULL, a default sequence will be used.
#' @param folds the number of folds for cross-validation; default is 5.
#' @param r the regularization parameter in MCP; default is 5.
#' @param clv a value or a vector, indexing variables that are not subject to penalty. clv only works for survival response for now,
#' and will be ignored for other types of responses.
#' @param initiation method for initiating the coefficient vector. For binary response, the default is elastic-net,
#' and for survival response the default is zero.
#' @param alpha.i the elastic-net mixing parameter. The program can use the elastic-net for choosing initial values of
#' the coefficient vector. alpha.i is the elastic-net mixing parameter, with \eqn{0 \le alpha.i \le 1}. \eqn{alpha.i=1} is the
#' lasso penalty, and \eqn{alpha.i=0} the ridge penalty. If user chooses a method other than elastic-net for initializing
#' coefficient, alpha.i will be ignored.
#' @param robust logical flag. Whether or not to use robust methods. Robust methods are only available for survival response.
#' @param standardize logical flag for standardizing variables in X; default is TRUE. If variables are already in the same units,
#' you might not wish to precede standardization.
#' @param verbo output progress to the console.
#' @return a list with components:
#' \item{lambda}{the optimal values(s) of \eqn{\lambda}. More than one values will be returned, if multiple lambdas have the cross-validated error =
#' min(cross-validated errors). If the network penalty was used, lambda contains optimal pair(s) of \eqn{\lambda 1} and \eqn{\lambda 2}.}
#' \item{mcvm}{the cross-validated error of the optimal \eqn{\lambda}.}
#' \item{CVM}{a matrix of the cross-validated errors of all lambdas used in the firs. The row names of CVM are the values of \eqn{\lambda 1}.
#' If the network penalty was used, the column names are the values of \eqn{\lambda 2}.}
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. BMC Genetics, 18(1):44.
#'
#' @seealso \code{\link{regnet}}
#'
#' @examples
#' data(LogisticExample)
#' X = rgn.logi$X
#' Y = rgn.logi$Y
#' out = cv.regnet(X, Y, response="binary", penalty="network", folds=5, r = 4.5)
#' out$lambda
#' b = regnet(X, Y, "binary", "network", out$lambda[1,1], out$lambda[1,2], r = 4.5)
#' index = which(rgn.logi$beta != 0)
#' pos = which(b != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export

cv.regnet <- function(X, Y, response=c("binary", "continuous", "survival"), penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL,
                      folds=5, r=NULL, clv=NULL, initiation=NULL, alpha.i=1, robust=TRUE, standardize=TRUE, verbo = FALSE)
{
  intercept = TRUE
  response = match.arg(response)
  penalty = match.arg(penalty)
  # method = paste(response, "_", penalty, sep = "")
  this.call = match.call()
  if(response == "survival"){
    if(ncol(Y) != 2) stop("y should be a two-column matrix")
    Y0 = Y[,"time"]
    status = Y[,"status"]
    if(sum(Y0<=0)>0) stop("Survival times need to be positive")
  }
  if(alpha.i>1 | alpha.i<0) stop("alpha.i should be between 0 and 1")
  folds = as.integer(folds)
  if(folds<2 | folds>ncol(X)) stop("incorrect value of folds")
  if(is.null(initiation)){
    if(response == "survival") initiation = "zero" else initiation = "elnet"
  }
  if(is.null(r)) r = 5
  if(penalty != "network") lamb.2 = 0
  alpha = alpha.i # temporary

  fit=switch (response,
    "binary" = CV.Logit(X, Y, penalty, lamb.1, lamb.2, folds, r, alpha, init=initiation, alpha.i, standardize, verbo),
    "survival" = CV.Surv(X, Y0, status, penalty, lamb.1, lamb.2, clv=clv, folds, r, init=initiation, alpha.i, robust, standardize, verbo)
    # "continuous_network" = NULL,
    # "continuous_mcp" = NULL,
    # "continuous_lasso" = NULL
  )
  fit$call = this.call
  class(fit) = "cv.regnet"
  fit
}
