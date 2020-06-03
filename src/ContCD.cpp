#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"ContCD.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;

void ContNet(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec &b, double r, arma::mat const &a, arma::vec const &u, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
	// double inp = arma::dot(x.col(m),x.col(m))/n, net1 = 0, net2 = 0;
	double net1 = 0;
    if(m < p-1){
      net1 = lam2 * as_scalar(a.submat(m, m+1, m, p-1) * b.subvec(m+1, p-1));
      // net2 = lam2 * arma::accu(arma::abs(a.row(m).subvec(m+1, p-1)));
    }
    // double z = arma::accu(x.col(m) % t)/n + net1;
	double z = arma::dot(x.col(m), t)/n + net1;
    // double u = inp + net2;
    if(std::abs(z) > (r*lam1*u(m))){
      b(m) = z/u(m);
    }else{
      b(m) = Soft(z, lam1)/(u(m) - 1/r);
    }
    t -= x.col(m) * b(m);
  }
  // return(b);
}

void ContMCP(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, double r, arma::vec const &inp, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    // double inp = arma::accu(square(x.col(m)))/n;
    // double z = arma::accu(x.col(m) % t)/n;
	double z = arma::dot(x.col(m), t)/n;

    if(std::abs(z) > (r*lam1*inp(m))){
      b(m) = z/inp(m);
    }else{
      b(m) = Soft(z, lam1)/(inp(m) - 1/r);
    }
    t -= x.col(m) * b(m);
  }
  // return(b);
}

void ContLasso(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, arma::vec const &inp, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    // double inp = arma::accu(square(x.col(m)))/n;
    // double z = arma::accu(x.col(m) % t)/n;
	double z = arma::dot(x.col(m), t)/n;

    b(m) = Soft(z, lam1)/inp(m);
    t -= x.col(m) * b(m);
  }
  // return(b);
}
