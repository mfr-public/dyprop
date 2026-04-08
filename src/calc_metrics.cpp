// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

struct MetricsWorker : public RcppParallel::Worker {
  const arma::mat& X_clr;
  const arma::mat& W_t;
  const arma::mat& Mu_X;
  const arma::mat& Var_X;
  const RcppParallel::RVector<int> geneA_idx;
  const RcppParallel::RVector<int> geneB_idx;
  double lambda_reg;
  int batch_size;
  int n_events;
  int n_time;

  // Output pointer to arma memory
  arma::uword* class_id_ptr;

  MetricsWorker(const arma::mat& X_clr, const arma::mat& W_t,
                const arma::mat& Mu_X, const arma::mat& Var_X,
                const IntegerVector& geneA_idx, const IntegerVector& geneB_idx,
                double lambda_reg, int batch_size, int n_events, int n_time,
                arma::uword* class_id_ptr)
      : X_clr(X_clr), W_t(W_t), Mu_X(Mu_X), Var_X(Var_X),
        geneA_idx(geneA_idx), geneB_idx(geneB_idx),
        lambda_reg(lambda_reg), batch_size(batch_size),
        n_events(n_events), n_time(n_time), class_id_ptr(class_id_ptr) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t b = begin; b < end; ++b) {
      int start_e = b * batch_size;
      int end_e = std::min(n_events, start_e + batch_size);
      int K = end_e - start_e;

      arma::mat lr_sq_mat(K, n_time);

      for (int k = 0; k < K; ++k) {
        int e = start_e + k;
        int gA = geneA_idx[e] - 1;
        int gB = geneB_idx[e] - 1;
        lr_sq_mat.row(k) = arma::square(X_clr.row(gA) - X_clr.row(gB));
      }

      // Level 3 BLAS GEMM (Matrix-Matrix product). 
      arma::mat E_lr_sq_mat = lr_sq_mat * W_t;

      for (int k = 0; k < K; ++k) {
        int e = start_e + k;
        int gA = geneA_idx[e] - 1;
        int gB = geneB_idx[e] - 1;

        arma::rowvec E_lr_sq = E_lr_sq_mat.row(k);
        arma::rowvec mu_lr = Mu_X.row(gA) - Mu_X.row(gB);
        arma::rowvec var_lr = E_lr_sq - arma::square(mu_lr);

        arma::rowvec denom = Var_X.row(gA) + Var_X.row(gB) + 1e-12; // Safety guard
        arma::rowvec phi_vec = var_lr / (denom + lambda_reg);
        arma::rowvec rho_vec = 1.0 - (var_lr / denom);

        double max_phi = phi_vec.max();
        double min_rho = rho_vec.min();

        if (max_phi < 0.2) {
          class_id_ptr[e] = 0; // Homeostasis
        } else {
          if (min_rho < 0.5) {
            class_id_ptr[e] = 2; // Decoupling
          } else {
            class_id_ptr[e] = 1; // Switch
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
List calc_metrics(arma::mat X_clr, IntegerVector geneA_idx,
                  IntegerVector geneB_idx, arma::vec pseudotime, double h_opt,
                  double lambda_reg = 1e-5) {

  int n_events = geneA_idx.length();
  int n_time = pseudotime.n_elem;
  int n_genes = X_clr.n_rows;

  arma::uvec class_id(n_events);

  // Pre-calculate Gaussian Weights Matrix W (Time x Time)
  arma::mat W(n_time, n_time);
  for (int i = 0; i < n_time; ++i) {
    double t0 = pseudotime(i);
    arma::rowvec w =
        arma::exp(-arma::pow(pseudotime.t() - t0, 2) / (2.0 * h_opt * h_opt));
    W.row(i) = w / arma::sum(w);
  }

  arma::mat W_t = W.t(); // Transpose once

  // Pre-compute Expected Values (Mu) and Expected Variances (Var) for ALL genes
  arma::mat Mu_X = X_clr * W_t;
  arma::mat Mu_X_sq = arma::square(X_clr) * W_t;
  arma::mat Var_X = Mu_X_sq - arma::square(Mu_X);

  int batch_size = 512;
  int n_batches = std::ceil((double)n_events / batch_size);

  // Dispatch Level 3 BLAS dynamically across available cores
  MetricsWorker worker(X_clr, W_t, Mu_X, Var_X, geneA_idx, geneB_idx, 
                       lambda_reg, batch_size, n_events, n_time, class_id.memptr());
  
  RcppParallel::parallelFor(0, n_batches, worker);

  return List::create(Named("Class_ID") = class_id);
}
