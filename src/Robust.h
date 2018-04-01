#ifndef ROBUST_h
#define ROBUST_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec LadNet(arma::mat& x, arma::vec& y, double lam1, double lam2, arma::vec b, double r, arma::mat& a, int n, int p);
arma::vec LadMCP(arma::mat& x, arma::vec& y, double lam1, arma::vec b, double r, int n, int p);
arma::vec LadLasso(arma::mat& x, arma::vec& y, double lam1, arma::vec b, int n, int p);

#endif