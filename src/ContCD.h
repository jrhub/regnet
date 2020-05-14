#ifndef CONTCD_h
#define CONTCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec ContNet(arma::mat  const &x, arma::vec  const &y, double lam1, double lam2, arma::vec b, double r, arma::mat  const &a, int n, int p);
arma::vec ContMCP(arma::mat  const &x, arma::vec  const &y, double lam1, arma::vec b, double r, int n, int p);
arma::vec ContLasso(arma::mat  const &x, arma::vec  const &y, double lam1, arma::vec b, int n, int p);

#endif