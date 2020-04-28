#' print a cv.regnet object
#'
#' Print a summary of a cv.regnet object
#'
#' @param x cv.regnet object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{cv.regnet}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{cv.regnet}}
#' @export
print.cv.regnet=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  cat("\nLambda:\n")
  print(x$lambda, digits)
  cat("\nCV error:\n")
  print(x$mcvm, digits)
  cat("\nPenalty:\n")
  print(x$penalty)

  # print(cbind(Df=x$df,"%Dev"=signif(x$dev.ratio,digits),Lambda=signif(x$lambda,digits)))
}



#' print a regnet object
#'
#' Print a summary of a regnet object
#'
#' @param x regnet object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{regnet}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{regnet}}
#' @export
print.regnet=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  cat("\nCoefficients:\n")
  print(x$coeff, digits)
  cat("Class:\n")
  print(class(x))
}
