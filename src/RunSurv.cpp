#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunSurv.h"
#include"QR.h"
#include"RobustCD.h"
#include"ContCD.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export]]
arma::vec RunSurv_robust(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, char method, bool debugging)
{
  int count = 0, n = xc.n_rows;
  arma::vec bold(p, fill::none), yc, yg;
  arma::mat const wc = arma::abs(xc)/n;
  arma::vec const totalWeights = arma::sum(wc, 0).t();
  arma::mat Wg(n+p, p, fill::zeros);
  Wg.rows(0,n-1) = arma::abs(xg)/n;

  while(count < 20){
	yc = y - xg * bg;
	// bc = QRWMR(xc, yc, bc);
	QRWMR(xc, yc, bc, wc, totalWeights);
	yg = y - xc * bc;
	bold = bg;
    if(method == 'n'){
      LadNet(xg, yg, lamb1, lamb2, bg, Wg, r, a, n, p, debugging);
    }else if(method == 'm'){
      LadMCP(xg, yg, lamb1, bg, r, n, p, debugging);
    }else{
      LadLasso(xg, yg, lamb1, bg, n, p, debugging);
    }
    double diff = arma::accu(arma::abs(bg - bold))/(arma::accu(bg != 0)+0.1);
	// double diff = arma::accu(arma::abs(bg - bold));
	// RcppThread::Rcout << "diff: " << diff <<std::endl;
    if(diff < 0.001) break;
    else{
      count++;
    }
  }
  // Rcpp::Rcout << "count: " << count << "\n";
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bg;
  return(b1);
}


void RunSurv_robust_warm(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec &bc, arma::vec &bg, double r, arma::mat const &a, int p, int pc, char method)
{
  int count = 0, n = xc.n_rows;
  arma::vec bold(p, fill::none), yc, yg; // bc = bc0, bg = bg0;
  arma::mat const wc = arma::abs(xc)/n;
  arma::vec const totalWeights = arma::sum(wc, 0).t();
  arma::mat Wg(n+p, p, fill::zeros);
  Wg.rows(0,n-1) = arma::abs(xg)/n;

  while(count < 20){
	yc = y - xg * bg;
	// bc = QRWMR(xc, yc, bc);
	QRWMR(xc, yc, bc, wc, totalWeights);
	yg = y - xc * bc;
	bold = bg;
    if(method == 'n'){
      LadNet(xg, yg, lamb1, lamb2, bg, Wg, r, a, n, p, false);
	  // LadNet(xg, yg, lamb1, lamb2, bg, r, a, n, p);
    }else if(method == 'm'){
      LadMCP(xg, yg, lamb1, bg, r, n, p, false);
    }else{
      LadLasso(xg, yg, lamb1, bg, n, p, false);
    }
    double diff = arma::accu(arma::abs(bg - bold))/(arma::accu(bg != 0)+0.1);
	// RcppThread::Rcout << "diff: " << diff <<std::endl;
    if(diff < 0.001) break;
    else{
      count++;
    }
  }
}


// [[Rcpp::export]]
arma::vec RunSurv(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, double lamb1, double lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, arma::vec const &triRowAbsSums, int p, int pc, char method)
{
  int count = 0, n = xc.n_rows;
  arma::vec bold(p, fill::none), yc, yg, u; // bc = bc0, bg = bg0;
  arma::vec const inp = arma::sum(arma::square(xg),0).t()/n;
  if(method == 'n') u = inp + lamb2 * triRowAbsSums;

  while(count < 20){
	yc = y - xg * bg;
	bc = fastLm(yc, xc);
	yg = y - xc * bc;
	bold = bg;
    if(method == 'n'){
      ContNet(xg, yg, lamb1, lamb2, bg, r, a, u, n, p);
    }else if(method == 'm'){
      ContMCP(xg, yg, lamb1, bg, r, inp, n, p);
    }else{
      ContLasso(xg, yg, lamb1, bg, inp, n, p);
    }
    double diff = arma::accu(arma::abs(bg - bold))/(arma::accu(bg != 0)+0.1);
    if(diff < 0.001) break;
    else{
      count++;
    }
  }
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bg;
  return(b1);
}
