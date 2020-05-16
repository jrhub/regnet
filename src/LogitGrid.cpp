// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunLogistic.h"
#include"Utilities.h"
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;

// [[Rcpp::export()]]
arma::mat LogitGrid(arma::mat& x, arma::vec& y, arma::mat& x2, arma::vec& y2, arma::vec lamb1, arma::vec lamb2, arma::vec b, double r, arma::mat& a, int p, double alpha, char method, int ncores)
{
	arma::vec btmp;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	RcppThread::checkUserInterrupt();
	
	omp_set_num_threads(ncores);
	#pragma omp parallel for collapse(2) private(btmp)
	for(unsigned int j=0; j<lamb2.n_elem ; j++){
		for(unsigned int i=0; i<lamb1.n_elem ; i++){
            if(method == 'n'){
                btmp = RunNet(x, y, lamb1(i), lamb2(j), b, r, a, p);
            }else if(method == 'm'){
                btmp = RunMCP(x, y, lamb1(i), b, r, p);
            }else{
                btmp = RunElastic(x, y, lamb1(i), b, alpha, p);
            }
            CVM(i, j) = validation_logit(x2, y2, btmp);
		}
	}
    return CVM;
}