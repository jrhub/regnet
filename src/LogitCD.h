#ifndef LOGITCD_h
#define LOGITCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double Soft(double z, double lambda);
arma::vec Network(arma::mat& x, arma::vec& y, double lam1, double lam2, arma::vec b, double r, arma::mat& a, int n, int p);
arma::vec MCP(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double r, int n, int p);
arma::vec Elastic(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double alpha, int n, int p);

#endif
