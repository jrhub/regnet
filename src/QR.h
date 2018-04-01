#ifndef QR_h
#define QR_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec QRWMR(arma::mat& x, arma::vec& y, arma::vec b);
arma::vec fastLm(const arma::vec & y, const arma::mat & X);

#endif