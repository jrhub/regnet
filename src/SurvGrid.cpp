
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunSurv.h"
#include"Utilities.h"

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
		
		Rcpp::Rcout << "CrossValidation: " << f+1 << "/" << folds << "\n";
		
		if(robust){
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				RunSurv_robust_warm(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, method);
				CVM(i, j) += validation_LAD(xc2, xg2, y2, bc, bg);
				Rcpp::checkUserInterrupt();
			}
		}
		}else{
			if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
			for(unsigned int j=0; j<lamb2.n_elem ; j++){
				for(unsigned int i=0; i<lamb1.n_elem ; i++){
					btmp = RunSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, triRowAbsSums, p, pc, method);
					CVM(i, j) += validation_LS(xc2, xg2, y2, btmp, p, pc);
					Rcpp::checkUserInterrupt();
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
				Rcpp::checkUserInterrupt();
			}
		}
	}else{
		if(method == 'n') triRowAbsSums = TriRowAbsSums(a);
		for(unsigned int j=0; j<lamb2.n_elem ; j++){
			for(unsigned int i=0; i<lamb1.n_elem ; i++){
				btmp = RunSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, triRowAbsSums, p, pc, method);
				CVM(i, j) = validation_LS(x2, y2, btmp);
				Rcpp::checkUserInterrupt();
			}
		}
	}
	
    return CVM;
}

/* 
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

 */