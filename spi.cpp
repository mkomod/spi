#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rbf(arma::vec x, arma::vec y, double sigma)
{
    int n = x.n_rows;
    int p = y.n_rows;
    arma::mat res = arma::mat(n, p, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < p; ++j) {
	    res(i, j) = exp(-pow(x(i) - y(j), 2) / (2.0 * sigma * sigma));
	}
    }
    return res;
}

