#ifndef CONTCD_h
#define CONTCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

void ContNet(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec &b, double r, arma::mat const &a, arma::vec const &u, int n, int p);

void ContMCP(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, double r, arma::vec const &inp, int n, int p);

void ContLasso(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, arma::vec const &inp, int n, int p);

#endif