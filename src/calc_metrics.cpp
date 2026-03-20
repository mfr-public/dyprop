// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List calc_metrics(arma::mat X_clr, IntegerVector geneA_idx,
                  IntegerVector geneB_idx, arma::vec pseudotime, double h_opt,
                  double lambda_reg = 1e-5) {

  // X_clr must be strictly Genes x Time
  int n_events = geneA_idx.length();
  int n_time = pseudotime.n_elem;

  // Output containers (Events x Time)
  arma::mat Phi_traces(n_events, n_time);
  arma::mat Rho_traces(n_events, n_time);

  // Pre-calculate Gaussian Weights Matrix W (Time x Time)
  // W(i, j) represents the weight of cell j when the kernel is centered at cell
  // i
  arma::mat W(n_time, n_time);
  for (int i = 0; i < n_time; ++i) {
    double t0 = pseudotime(i);
    arma::rowvec w =
        arma::exp(-arma::pow(pseudotime.t() - t0, 2) / (2.0 * h_opt * h_opt));
    W.row(i) = w / arma::sum(w); // Ensure weights sum directly to 1
  }

// Vectorized computation of variances per event trajectory via OpenMP
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int e = 0; e < n_events; ++e) {
    // R indices are 1-based, converting to C++ 0-based indexing
    int gA = geneA_idx(e) - 1;
    int gB = geneB_idx(e) - 1;

    arma::rowvec valA = X_clr.row(gA);
    arma::rowvec valB = X_clr.row(gB);
    arma::rowvec lr = valA - valB;

    for (int t_idx = 0; t_idx < n_time; ++t_idx) {
      arma::rowvec w = W.row(t_idx);

      // Expected Value (Mu)
      double mu_A = arma::sum(w % valA);
      double mu_B = arma::sum(w % valB);
      double mu_lr = arma::sum(w % lr);

      // Expected Variance
      double var_A = arma::sum(w % arma::square(valA - mu_A));
      double var_B = arma::sum(w % arma::square(valB - mu_B));
      double var_lr = arma::sum(w % arma::square(lr - mu_lr));

      double denom = var_A + var_B + 1e-12; // Safety guard preventing NaN

      // Dynamic Instability
      Phi_traces(e, t_idx) = var_lr / (denom + lambda_reg);

      // Dynamic Coupling (Proportionality)
      Rho_traces(e, t_idx) = 1.0 - (var_lr / denom);
    }
  }

  return List::create(Named("Phi") = Phi_traces, Named("Rho") = Rho_traces);
}
