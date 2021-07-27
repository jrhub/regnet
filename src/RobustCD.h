#ifndef ROBUSTCD_h
#define ROBUSTCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

void LadNet(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec &b, arma::mat &Wg, double r, arma::mat const &a, int n, int p, bool debugging);
void LadNet(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec &b, double r, arma::mat const &a, int n, int p, bool debugging);

void LadMCP(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, double r, int n, int p, bool debugging);
void LadLasso(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, int n, int p, bool debugging);

#endif