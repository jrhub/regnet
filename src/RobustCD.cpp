#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RobustCD.h"

using namespace Rcpp;
using namespace arma;

//Robust Network
void LadNet(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec &b, arma::mat &Wg, double r, arma::mat const &a, int n, int p, bool debugging)
{
  arma::vec t = y - x * b;
  arma::uvec index;
  
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    arma::vec u(n+p-m, fill::zeros);
    u.subvec(0,n-1) = t/x.col(m);
    
	Wg(n,m) = std::max(0.0, lam1 - std::abs(b(m))/r);
	
    if(m < p-1){
	  u.subvec(n+1, n+p-m-1) = arma::sign(a.submat(m, m+1, m, p-1)).t() % b.subvec(m+1, p-1);
	  Wg.submat(n+1, m, n+p-m-1, m) = arma::abs(a.submat(m, m+1, m, p-1)).t() * lam2;
    }
	
	if(debugging) u.replace(datum::nan, 0);
    index = arma::sort_index(u);
    double halfWeight = arma::accu(Wg.col(m))*0.5, SUM = 0;
    int j = 0;
    do{
	  SUM += Wg(index(j),m);
      j++;
    }while(SUM <= halfWeight);
    b(m) = u(index(j-1));
    t -= x.col(m) * b(m);
  }
}

void LadNet(arma::mat const &x, arma::vec const &y, double lam1, double lam2, arma::vec &b, double r, arma::mat const &a, int n, int p, bool debugging)
{
  arma::vec t = y - x * b;
  arma::uvec index;
  
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    arma::vec u(n+p-m, fill::zeros), w(n+p-m, fill::zeros);
    u.subvec(0,n-1) = t/x.col(m);
    
    w.subvec(0, n-1) = arma::abs(x.col(m))/n;
	w(n) = lam1 - std::abs(b(m))/r;
	if(w(n)<0) w(n) = 0;
	
    if(m < p-1){
	  u.subvec(n+1, n+p-m-1) = arma::sign(a.row(m).subvec(m+1, p-1)).t() % b.subvec(m+1, p-1);
      w.subvec(n+1, n+p-m-1) = arma::abs(a.row(m).subvec(m+1, p-1)).t() * lam2;
    }

	if(debugging) u.replace(datum::nan, 0);
    index = arma::sort_index(u);
    double TotalWeight = arma::accu(w), SUM = 0;
    int j = 0;
    do{
	  SUM += w(index(j))/TotalWeight;
      j++;
    }while(SUM <= 0.5);
    b(m) = u(index(j-1));
    t -= x.col(m) * b(m);
  }
}


//Robust MCP
void LadMCP(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, double r, int n, int p, bool debugging)
{
  arma::vec t = y - x * b;
  arma::uvec index;
  arma::vec u(n+1, fill::zeros), w(n+1, fill::zeros);
  
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    u.subvec(0,n-1) = t/x.col(m);
	if(debugging) u.replace(datum::nan, 0);
    // u(n) = 0;
    
    w.subvec(0, n-1) = arma::abs(x.col(m))/n;
    w(n) = lam1 - std::abs(b(m))/r;
	if(w(n)<0) w(n) = 0;
    
    index = arma::sort_index(u);
    double TotalWeight = arma::accu(w), SUM = 0;
    int j = 0;
    do{
      SUM += w(index(j))/TotalWeight;
      j++;
    }while(SUM <= 0.5);
    b(m) = u(index(j-1));
    
    t -= x.col(m) * b(m);
  }
  // return(b);
}

//Robust Lasso
void LadLasso(arma::mat const &x, arma::vec const &y, double lam1, arma::vec &b, int n, int p, bool debugging)
{
  arma::vec t = y - x * b;
  arma::uvec index;
  arma::vec u(n+1, fill::zeros), w(n+1, fill::zeros);
  
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    u.subvec(0,n-1) = t/x.col(m);
	if(debugging) u.replace(datum::nan, 0);
    // u(n) = 0;
    
    w.subvec(0, n-1) = arma::abs(x.col(m))/n;
    w(n) = lam1;
    
    index = arma::sort_index(u);
    // arma::vec _w = w(index), _u = u(index);
    double TotalWeight = arma::accu(w), SUM = 0;
    int j = 0;
    do{
      SUM += w(index(j))/TotalWeight;
      j++;
    }while(SUM <= 0.5);
    b(m) = u(index(j-1));
    
    t -= x.col(m) * b(m);
  }
  // return(b);
}
