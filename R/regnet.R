#' fit a regression for given lambda with network-based regularization
#'
#' Fit a network-based regression for given values of lambda1 and lambda2.
#' Typical usage is to have the cv.regnet function compute the optimal lambdas, then provide them to the
#' regnet function. User can also choose to use MCP or Lasso regularization.
#'
#' @param X matrix of predictors withour intercept. Each row should be an observation vector. The column of 1 will be added to the X matric
#' by the programme as the intercept.
#' @param Y response vairable. For response="binary" should be a numeric vector with zeros and ones. For response="survival", y should be a
#' two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0'
#' indicating right censored.
#' @param response response type. regnet now supports two types of response: "binary" and "survival".
#' @param penalty penalty type. regnet provides three choices of penalty: "network", "mcp" and "lasso".
#' @param lamb.1 the tuning parameter \eqn{\lambda 1} that imposes sparsity.
#' @param lamb.2 the tuning parameter \eqn{\lambda 2} that controls the smoothness among coefficient profiles. Only needed for network penalty.
#' @param r the regularization parameter in MCP.
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
#'
#' @return the estimated coefficients vector.
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. BMC Genetics, 18(1):44.
#'
#' @seealso \code{\link{cv.regnet}}
#'
#' @examples
#' data(SurvExample)
#' X = rgn.surv$X
#' Y = rgn.surv$Y
#' clv = c(1:5) # variable 1 to 5 are clinical variables, we choose not to penalize them here.
#' lambda1 = 5*10^-4
#' lambda2 = 10^-3
#' b = regnet(X, Y, "survival", "network", lambda1, lambda2, clv=clv, robust=TRUE)
#' index = which(rgn.surv$beta != 0)
#' pos = which(b != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export

regnet <- function(X, Y, response=c("continuous", "binary", "survival"), penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL,
                   r=NULL, clv=NULL, initiation=NULL, alpha.i=1, robust=TRUE, standardize=TRUE)
{
  intercept = TRUE
  response = match.arg(response)
  penalty = match.arg(penalty)
  if(penalty != "network") lamb.2 = 0
  if(missing(lamb.1)) stop("Both lambda1 and lambda2 need to be provided")
  if(missing(lamb.2)) stop("Lambda2 needs to be provided for network method")

  this.call = match.call()
  if(response == "survival"){
    if(ncol(Y) != 2) stop("y should be a two-column matrix")
    Y0 = Y[,"time"]
    status = Y[,"status"]
    if(sum(Y0<=0)>0) stop("Survival times need to be positive")
  }
  if(alpha.i>1 | alpha.i<0) stop("alpha.i should be between 0 and 1")
  if(is.null(initiation)){
    if(response == "survival") initiation = "zero" else initiation = "elnet"
  }
  if(is.null(r)) r = 5
  alpha = alpha.i # temporary

  fit=switch (response,
              "binary" = LogitCD(X, Y, penalty, lamb.1, lamb.2, r, alpha, init=initiation, alpha.i,standardize),
              "survival" = SurvCD(X, Y0, status, penalty, lamb.1, lamb.2, clv, r, init=initiation, alpha.i, robust, standardize)
              # "continuous" = NULL
  )
  # fit$call = this.call
  # class(fit) = "regnet"
  fit
}
