#ifndef UTILITIES_h
#define UTILITIES_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double Soft(double z, double lambda);
double validation_LS(arma::mat& x, arma::vec& y, arma::vec& b);
double validation_LAD(arma::mat& x, arma::vec& y, arma::vec& b);
double validation_logit(arma::mat& x0, arma::vec& y0, arma::vec& b);

#endif