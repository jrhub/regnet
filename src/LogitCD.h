#ifndef LOGITCD_h
#define LOGITCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec Network(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec b, double r, arma::mat const &a, arma::vec const &u, int n, int p);

arma::vec Network(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p);

arma::vec MCP(arma::mat const &x, arma::vec const &y, double lambda, arma::vec b, double r, int n, int p);

arma::vec Elastic(arma::mat const &x, arma::vec const &y, double lambda, arma::vec b, double alpha, int n, int p);

#endif
