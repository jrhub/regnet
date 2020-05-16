// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"ContCD.h"
#include"Robust.h"
#include"QR.h"
#include"RunCont.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec RunCont(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, 
					arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method)
{
  int count = 0, n = xc.n_rows;
  arma::vec bnew(p, fill::none), yc, yg; // bc = bc0, bg = bg0;

  while(count < 20){
	yc = y - xg * bg;
	bc = fastLm(yc, xc);
	yg = y - xc * bc;
    if(method == 'n'){
      bnew = ContNet(xg, yg, lamb1, lamb2, bg, r, a, n, p);
    }else if(method == 'm'){
      bnew = ContMCP(xg, yg, lamb1, bg, r, n, p);
    }else{
      bnew = ContLasso(xg, yg, lamb1, bg, n, p);
    }
    double diff = arma::accu(arma::abs(bg - bnew))/(arma::accu(bg != 0)+0.1);
	//std::cout << "diff: " << diff <<std::endl;
    if(diff < 0.001) break;
    else{
      bg = bnew;
      count++;
    }
  }
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bnew;
  return(b1);
}


// [[Rcpp::export]]
arma::vec RunCont_robust(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, 
					arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method)
{
  int count = 0, n = xc.n_rows;
  arma::vec bnew(p, fill::none), yc, yg; // bc = bc0, bg = bg0;

  while(count < 20){
	yc = y - xg * bg;
	bc = QRWMR(xc, yc, bc);
	yg = y - xc * bc;
    if(method == 'n'){
      bnew = LadNet(xg, yg, lamb1, lamb2, bg, r, a, n, p);
    }else if(method == 'm'){
      bnew = LadMCP(xg, yg, lamb1, bg, r, n, p);
    }else{
      bnew = LadLasso(xg, yg, lamb1, bg, n, p);
    }
    double diff = arma::accu(arma::abs(bg - bnew))/(arma::accu(bg != 0)+0.1);
	//std::cout << "diff: " << diff <<std::endl;
    if(diff < 0.001) break;
    else{
      bg = bnew;
      count++;
    }
  }
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bnew;
  return(b1);
}









