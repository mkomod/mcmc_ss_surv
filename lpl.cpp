#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_PL(arma::mat X, arma::vec b, arma::uvec Y_sorted, arma::uvec Y_failure)
{
    arma::vec xb = X * b;
    double a = max(xb);

    double risk = 0.0;
    double tot = 0.0;
    int ind_last_failure = X.n_rows;
    int ind_curr_failure = 0;

    for (int i = (Y_failure.n_rows - 1); i >= 0; --i) {
	ind_curr_failure = Y_failure(i);
	auto R = arma::span(ind_curr_failure, ind_last_failure - 1);

	// compute the denom and num
	risk += sum(exp(xb(Y_sorted(R)) - a));
	tot += xb(Y_sorted(ind_curr_failure)) - (a + log(risk));

	ind_last_failure = ind_curr_failure;
    }
    return tot;
}

// [[Rcpp::export]]
double log_Laplace(double beta, double lambda)
{
    return log(lambda) + log(0.5) - lambda * std::abs(beta);
}

// [[Rcpp::export]] 
double sigmoid(double x)
{
    double res = 0.0;
    if (x >= 0) {
	res = 1.0 / (1.0 + exp(-x));	
    } else {
	res = exp(x) / (1.0 + exp(x));
    }
    return res;
}
