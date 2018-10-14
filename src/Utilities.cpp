#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;


double Soft(double z, double lambda){
  if(z > lambda) return(z - lambda);
  else if(z < -lambda) return(z + lambda);
  else return(0);
}

//Non-robust MSE
double validation_LS(arma::mat& x, arma::vec& y, arma::vec& b){
	//double mse = arma::accu(pow(y - x*b, 2))/y.n_elem;
    double mse = arma::accu(arma::square(y - x*b));
	return(mse);
}

double validation_LAD(arma::mat& x, arma::vec& y, arma::vec& b){
	//double lad = accu(arma::abs(y - x*b))/y.n_elem;
    double lad = accu(arma::abs(y - x*b));
	return(lad);
}

double validation_logit(arma::mat& x0, arma::vec& y0, arma::vec& b){
    arma::vec yi = 1/(1 + arma::exp(x0*b*-1));
    arma::uvec y = (yi > 0.5);
    double mc = arma::accu(arma::abs(y0 - y));
	return(mc);
}
coeff3.subvec
vec fastLm(const arma::vec & y, const arma::mat & X){
    arma::colvec coef = arma::solve(X, y);    
    return coef;
}