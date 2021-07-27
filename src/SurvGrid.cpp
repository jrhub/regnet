// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunSurv.h"
#include"Utilities.h"
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
arma::mat SurvCV(arma::mat const &Xc, arma::mat const &Xg, arma::vec const &Y, unsigned int folds, arma::vec lamb1, arma::vec lamb2, arma::vec const &bc0, arma::vec const &bg0, double r, arma::mat const &a, int p, int pc, bool robust, char method, int ncores, bool debugging)
{
	int L = std::floor(Xc.n_rows/(double)folds), m = Xc.n_rows%folds, start=0;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros), xc, xg, xc2, xg2;
	arma::vec btmp, bc, bg, y, y2, triRowAbsSums;
	Rcpp::IntegerVector sample_seq = Rcpp::rev(Rcpp::seq(0,Xc.n_rows-1)), temp;

	for(unsigned int f=0; f<folds ; f++){
		if(f<(folds-m)){
			temp = Rcpp::seq(start, start+L-1);
			xc2 = Xc.rows(start, start+L-1);
			xg2 = Xg.rows(start, start+L-1);
			y2 = Y.subvec(start, start+L-1);
			start += L;
		}else{
			temp = Rcpp::seq(start, start+L);
			xc2 = Xc.rows(start, start+L);
			xg2 = Xg.rows(start, start+L);
			y2 = Y.subvec(start, start+L);
			start += L+1;
		}
		
		arma::uvec train = Rcpp::as<arma::uvec>(Rcpp::setdiff(sample_seq, temp));
		xc = Xc.rows(train);
		xg = Xg.rows(train);
		y = Y.elem(train);
		
		RcppThread::Rcout << "CrossValidation: " << f+1 << "/" << folds << "\n";
		// RcppThread::Rcout << "train: " << train.t() << "\n";
		
		// for(unsigned int j=0; j<lamb2.n_elem ; j++){
			// bc = bc0; bg = bg0;
			// for(unsigned int i=0; i<lamb1.n_elem ; i++){
				// if(robust){
					// RunSurv_robust_warm(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, method);
					// CVM(i, j) += validation_LAD(xc2, xg2, y2, bc, bg);
				// }else{
					// btmp = RunSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, method);
					// CVM(i, j) += validation_LS(xc2, xg2, y2, btmp, p, pc);
				// }
				// RcppThread::checkUserInterrupt();
			// }
		// }
		
		if(robust){
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				RunSurv_robust_warm(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, method);
				CVM(i, j) += validation_LAD(xc2, xg2, y2, bc, bg);
				RcppThread::checkUserInterrupt();
			}
		}
		}else{
			if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
			for(unsigned int j=0; j<lamb2.n_elem ; j++){
				for(unsigned int i=0; i<lamb1.n_elem ; i++){
					btmp = RunSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, triRowAbsSums, p, pc, method);
					CVM(i, j) += validation_LS(xc2, xg2, y2, btmp, p, pc);
					RcppThread::checkUserInterrupt();
				}
			}
		}
	
	}
    return CVM;
}


// [[Rcpp::export()]]
arma::mat SurvGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec const &lamb1, arma::vec const &lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, bool robust, char method, bool debugging)
{
	arma::vec btmp, triRowAbsSums;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	
	if(robust){
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				btmp = RunSurv_robust(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, method, debugging);
				CVM(i, j) = validation_LAD(x2, y2, btmp);
				RcppThread::checkUserInterrupt();
			}
		}
	}else{
		if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				btmp = RunSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, triRowAbsSums, p, pc, method);
				CVM(i, j) = validation_LS(x2, y2, btmp);
				RcppThread::checkUserInterrupt();
			}
		}
	}
	
    return CVM;
}

// [[Rcpp::export()]]
arma::mat SurvGrid_MC(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec const &lamb1, arma::vec const &lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, bool robust, char method, int ncores, bool debugging)
{
	arma::vec btmp, triRowAbsSums;
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	omp_set_num_threads(ncores);
	
	if(robust){
		#pragma omp parallel for collapse(2) private(btmp)
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				btmp = RunSurv_robust(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, method, debugging);
				CVM(i, j) = validation_LAD(x2, y2, btmp);
			}
		}
	}else{
		if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
		#pragma omp parallel for collapse(2) private(btmp)
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				btmp = RunSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, triRowAbsSums, p, pc, method);
				CVM(i, j) = validation_LS(x2, y2, btmp);
			}
		}
	}
	
    return CVM;
}

          
// // [[Rcpp::export()]]
// arma::mat NetGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec lamb1, arma::vec lamb2, arma::vec bc, arma::vec bg, double r, arma::mat const &a, int p, int pc, bool robust, int ncores)
// {
	// arma::vec bnet;
    // arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	
	// omp_set_num_threads(ncores);
    // // #pragma omp parallel for collapse(2) private(bnet)
	// for(unsigned int j=0; j<lamb2.n_elem ; j++){
		// #pragma omp parallel for private(bnet)
		// for(unsigned int i=0; i<lamb1.n_elem ; i++){
			// bnet = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, robust);
            // //bnet = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), bnet.subvec(0, pc-1), bnet.subvec(pc, pc+p-1), r, a, p, pc, robust);
			// if(robust){
				// CVM(i, j) = validation_LAD(x2, y2, bnet);
			// }else{
				// CVM(i, j) = validation_LS(x2, y2, bnet);
			// }
		// }
		// RcppThread::checkUserInterrupt();
	// }
    // return CVM;
// }


// // [[Rcpp::export()]]
// arma::vec MCPGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec lamb1, arma::vec bc, arma::vec bg, double r, int p, int pc, bool robust, int ncores)
// {
	// // arma::vec b((pc+p), fill::zeros);
	// arma::vec b;
    // arma::vec CVM(lamb1.n_elem, fill::zeros);
	// RcppThread::checkUserInterrupt();
	
	// omp_set_num_threads(ncores);
    // #pragma omp parallel for private(b)
    // for(unsigned int i=0; i<lamb1.n_elem ; i++){
        // b = RunMCPSurv(xc, xg, y, lamb1(i), bc, bg, r, p, pc, robust);
        // if(robust){
            // CVM(i) = validation_LAD(x2, y2, b);
        // }else{
            // CVM(i) = validation_LS(x2, y2, b);
        // }
    // }
    // return CVM;
// }

// // [[Rcpp::export()]]
// arma::vec LassoGrid(arma::mat const &xc, arma::mat const &xg, arma::vec const &y, arma::mat const &x2, arma::vec const &y2, arma::vec lamb1, arma::vec bc, arma::vec bg, int p, int pc, bool robust, int ncores)
// {
	// arma::vec b;
    // arma::vec CVM(lamb1.n_elem, fill::zeros);
	// RcppThread::checkUserInterrupt();
	
	// omp_set_num_threads(ncores);
    // #pragma omp parallel for private(b)
    // for(unsigned int i=0; i<lamb1.n_elem ; i++){
        // b = RunLassoSurv(xc, xg, y, lamb1(i), bc, bg, p, pc, robust);
        // if(robust){
            // CVM(i) = validation_LAD(x2, y2, b);
        // }else{
            // CVM(i) = validation_LS(x2, y2, b);
        // }
    // }
    // return CVM;
// }
