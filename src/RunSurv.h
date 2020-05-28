#ifndef RUNSURV_h
#define RUNSURV_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec RunSurv_robust(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method);

void RunSurv_robust_inplace(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec &bc, arma::vec &bg, double r, arma::mat const &a, int p, int pc, char method);

arma::vec RunSurv(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method);

#endif