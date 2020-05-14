#ifndef QR_h
#define QR_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec QRWMR(arma::mat const &x, arma::vec const &y, arma::vec b);

#endif