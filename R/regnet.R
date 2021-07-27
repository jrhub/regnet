#' fit a regression for given lambda with network-based regularization
#'
#' Network-based penalization regression for given values of \eqn{\lambda_{1}} and \eqn{\lambda_{2}}.
#' Typical usage is to have the cv.regnet function compute the optimal lambdas, then provide them to the
#' regnet function. Users could also use MCP or Lasso.
#'
#' @keywords models
#' @param X matrix of predictors without intercept. Each row should be an observation vector. A column of 1 will be added to the X matrix
#' by the program as the intercept.
#' @param Y response variable. For response="binary", Y should be a numeric vector with zeros and ones. For response="survival", Y should be a
#' two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating a event, and '0'
#' indicating censoring.
#' @param response response type. regnet now supports three types of response: "binary", "continuous" and "survival".
#' @param penalty penalty type. regnet provides three choices for the penalty function: "network", "mcp" and "lasso".
#' @param lamb.1 the tuning parameter \eqn{\lambda_{1}} that imposes sparsity.
#' @param lamb.2 the tuning parameter \eqn{\lambda_{2}} that controls the smoothness among coefficient profiles. \eqn{\lambda_{2}} is  needed
#' for network penalty.
#' @param r the regularization parameter in MCP. For binary response, r should be larger than 4.
#' @param clv a value or a vector, indexing variables that are not subject to penalty. clv only works for continuous and survival responses
#' for now, and will be ignored for other types of responses.
#' @param initiation method for initiating the coefficient vector. The default method is elastic-net.
#' @param alpha.i the elastic-net mixing parameter. The program can use the elastic-net for choosing initial values of
#' the coefficient vector. alpha.i is the elastic-net mixing parameter, with 0 \eqn{\le} alpha.i \eqn{\le} 1. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 is the ridge penalty. If the user chooses a method other than elastic-net for initializing
#' coefficients, alpha.i will be ignored.
#' @param robust logical flag. Whether or not to use robust methods. Robust methods are only available for survival response.
#' @param debugging logical flag. If TRUE, extra information will be returned.
#'
#' @details The current version of regnet supports two types of responses: “binary”, "continuous" and “survival”.
#' \itemize{
#' \item {regnet(…, response="binary", penalty="network")} fits a network-based penalized logistic regression.
#' \item {regnet(…, response="continuous", penalty="network")} fits a network-based least square regression.
#' \item {regnet(…, response="survival", penalty="network")} fits a robust regularized AFT model using network penalty.
#' }
#' Please see the references for more details of the models. By default, regnet uses robust methods for survival response.
#' If users would like to use non-robust methods, simply set robust=FALSE. User could also use MCP or Lasso penalty.
#'
#' The coefficients are always estimated on a standardized X matrix. regnet standardizes each columns of X to have unit variance
#' (using 1/n rather than 1/(n-1) formula). If the coefficients on the original scale are needed, the user can refit a standard model
#' using the subset of variables that have non-zero coefficients.
#'
#' @return an object of class "regnet" is returned, which is a list with components:
#' \item{coeff}{a vector of estimated coefficients. Please note that, if there are variables not subject to penalty (indicated by clv),
#' the order of returned vector is c(Intercept, unpenalized coefficients of clv variables, penalized coefficients of other variables).}
#' \item{Adj}{a matrix of adjacency measures of the identified genetic variants. Identified genetic variants are those that have non-zero estimated coefficients.}
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., and Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. {\emph{BMC Genetics}, 18(1):44} \doi{10.1186/s12863-017-0495-5}
#'
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang,Y. and Wu, C. (2019). Robust network-based regularization
#' and variable selection for high dimensional genomics data in cancer prognosis.
#' {\emph{Genet. Epidemiol.}, 43:276-291} \doi{10.1002/gepi.22194}
#'
#' @seealso \code{\link{cv.regnet}}
#'
#' @examples
#' ## Survival response
#' data(SurvExample)
#' X = rgn.surv$X
#' Y = rgn.surv$Y
#' clv = c(1:5) # variables 1 to 5 are clinical variables which we choose not to penalize.
#' penalty = "network"
#' fit = regnet(X, Y, "survival", penalty, rgn.surv$lamb1, rgn.surv$lamb2, clv=clv, robust=TRUE)
#' index = which(rgn.surv$beta != 0)
#' pos = which(fit$coeff != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export

regnet <- function(X, Y, response=c("binary", "continuous", "survival"), penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL,
                   r=NULL, clv=NULL, initiation=NULL, alpha.i=1, robust=FALSE, debugging = FALSE)
{
  # intercept = TRUE
  standardize=TRUE
  response = match.arg(response)
  penalty = match.arg(penalty)
  X = as.matrix(X)
  if(penalty != "network"){
    lamb.2 = 0
  }else{
    if(ncol(X)<3) stop("too less variables for network penalty.")
  }
  if(missing(lamb.1)) stop("Both lambda1 and lambda2 need to be provided")
  if(missing(lamb.2)) stop("Lambda2 needs to be provided for network method")

  this.call = match.call()
  if(response == "survival"){
    if(ncol(Y) != 2) stop("Y should be a two-column matrix")
    if(!setequal(colnames(Y), c("time", "status"))) stop("Y should be a two-column matrix with columns named 'time' and 'status'")
    Y0 = Y[,"time"]
    status = as.numeric(Y[,"status"])
    if(any(Y0<=0)) stop("survival times need to be positive")
    if(length(Y0) != nrow(X))  stop("the number of rows of Y does not match the number of rows of X");
    if(!all(status%in% c(0,1))) stop("status has to be a binary variable of 1 and 0")
  }else{
    if(length(Y) != nrow(X))  stop("length of Y does not match the number of rows of X");
  }

  if(response=="binary"){
    if(!all(Y%in% c(0,1))) stop("Y has to be a binary variable of 1 and 0.")
    if(robust) message("robust methods are not available for ", response, " response.")
  }
  if(alpha.i>1 | alpha.i<0) stop("alpha.i should be between 0 and 1")
  # if(is.null(initiation)){
  #   if(response == "survival") initiation = "zero" else initiation = "elnet"
  # }

  if(is.null(r)) r = 5
  alpha = alpha.i # temporary

  out=switch (response,
              "binary" = LogitCD(X, Y, penalty, lamb.1, lamb.2, r, alpha, init=initiation, alpha.i,standardize),
              "continuous" = ContCD(X, Y, penalty, lamb.1, lamb.2, clv, r, alpha, init=initiation, alpha.i, robust, standardize, debugging),
              "survival" = SurvCD(X, Y0, status, penalty, lamb.1, lamb.2, clv, r, init=initiation, alpha.i, robust, standardize, debugging)
  )
  para = list(penalty=penalty, lamb.1=lamb.1, lamb.2=lamb.2, robust=robust)
  fit = list(call = this.call, coeff = out$b, Adj=out$Adj, para=para)
  class(fit) = "regnet"
  fit
}
