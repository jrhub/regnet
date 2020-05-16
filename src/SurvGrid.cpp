// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"SurvCD.h"
#include"Utilities.h"
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;
          
// [[Rcpp::export()]]
arma::mat NetGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec lamb1, arma::vec lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, bool robust, int ncores)
{
	arma::vec bnet;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	
	omp_set_num_threads(ncores);
    // #pragma omp parallel for collapse(2) private(bnet)
	for(unsigned int j=0; j<lamb2.n_elem ; j++){
		#pragma omp parallel for private(bnet)
		for(unsigned int i=0; i<lamb1.n_elem ; i++){
			bnet = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, robust);
            //bnet = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), bnet.subvec(0, pc-1), bnet.subvec(pc, pc+p-1), r, a, p, pc, robust);
			if(robust){
				CVM(i, j) = validation_LAD(x2, y2, bnet);
			}else{
				CVM(i, j) = validation_LS(x2, y2, bnet);
			}
		}
		RcppThread::checkUserInterrupt();
	}
    return CVM;
}


// [[Rcpp::export()]]
arma::vec MCPGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec lamb1, arma::vec bc, arma::vec bg, double r, int p, int pc, bool robust, int ncores)
{
	// arma::vec b((pc+p), fill::zeros);
	arma::vec b;
    arma::vec CVM(lamb1.n_elem, fill::zeros);
	RcppThread::checkUserInterrupt();
	
	omp_set_num_threads(ncores);
    #pragma omp parallel for private(b)
    for(unsigned int i=0; i<lamb1.n_elem ; i++){
        b = RunMCPSurv(xc, xg, y, lamb1(i), bc, bg, r, p, pc, robust);
        if(robust){
            CVM(i) = validation_LAD(x2, y2, b);
        }else{
            CVM(i) = validation_LS(x2, y2, b);
        }
    }
    return CVM;
}

// [[Rcpp::export()]]
arma::vec LassoGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec lamb1, arma::vec bc, arma::vec bg, int p, int pc, bool robust, int ncores)
{
	arma::vec b;
    arma::vec CVM(lamb1.n_elem, fill::zeros);
	RcppThread::checkUserInterrupt();
	
	omp_set_num_threads(ncores);
    #pragma omp parallel for private(b)
    for(unsigned int i=0; i<lamb1.n_elem ; i++){
        b = RunLassoSurv(xc, xg, y, lamb1(i), bc, bg, p, pc, robust);
        if(robust){
            CVM(i) = validation_LAD(x2, y2, b);
        }else{
            CVM(i) = validation_LS(x2, y2, b);
        }
    }
    return CVM;
}
