#ifndef RUNLOGISTIC_h
#define RUNLOGISTIC_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec RunLogit(arma::mat const &x, arma::vec const &y, double lamb1, double lamb2, arma::vec b, double r, arma::mat const &a, arma::vec const &triRowAbsSums, int p, double alpha, char method);

arma::vec RunNet(arma::mat& x, arma::vec& y, double lamb1, double lamb2, arma::vec b, double r, arma::mat& a, int p);
arma::vec RunMCP(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double r, int p);
arma::vec RunElastic(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double alpha, int p);

#endif
