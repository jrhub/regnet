#' Example datasets for demonstrating the features of regnet
#'
#' Example datasets for demonstrating the features of regnet.
#'
#' @docType data
#' @keywords datasets
#' @name rgn
#' @aliases rgn.logi rgn.surv rgn.tcga rgn.htr
#' @usage data("LogisticExample")
#' data("SurvExample")
#' data("ContExample")
#' data("HeteroExample")
#' @format "LogisticExample", "SurvExample" and "HeteroExample" are simulated data. Each data includes three main components: X, Y, and beta; beta is a vector of the
#' true coefficients used to generate Y.
#'
#' "ContExample" is a subset of the skin cutaneous melanoma data from the Cancer Genome Atlas (TCGA). The response variable Y is
#' the log-transformed Breslow’s depth. X is a matrix of gene expression data.
#'
#' @examples
#' data("LogisticExample")
#' lapply(rgn.logi, class)
NULL
