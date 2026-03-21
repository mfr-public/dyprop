// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// HELPER: Inline Winsorization and Z-Score
// -----------------------------------------------------------------------------
void winsorize_and_zscore(arma::rowvec &y) {
  int n = y.n_elem;
  if (n < 3)
    return; // safeguard

  // Find 1st and 99th percentiles (O(N) using nth_element)
  int idx_low = std::max(0, (int)(0.01 * n));
  int idx_high = std::min(n - 1, (int)(0.99 * n));

  arma::rowvec y_copy = y;
  std::nth_element(y_copy.begin(), y_copy.begin() + idx_low, y_copy.end());
  double val_low = y_copy[idx_low];

  std::nth_element(y_copy.begin(), y_copy.begin() + idx_high, y_copy.end());
  double val_high = y_copy[idx_high];

  // Clip outliers robustly
  y.elem(arma::find(y < val_low)).fill(val_low);
  y.elem(arma::find(y > val_high)).fill(val_high);

  // Standardize to Z-Scores
  double mu = arma::mean(y);
  double sigma = arma::stddev(y);
  if (sigma > 1e-12) {
    y = (y - mu) / sigma;
  } else {
    y.zeros();
  }
}

// -----------------------------------------------------------------------------
// HELPER: Generalized Logistic Function (Switch)
// -----------------------------------------------------------------------------
arma::rowvec logistic_curve(const arma::vec &t, double tau, double epsilon) {
  return (1.0 / (1.0 + arma::exp(-(t - tau) / epsilon))).t();
}

// -----------------------------------------------------------------------------
// HELPER: Gaussian Pulse Function (Transient)
// -----------------------------------------------------------------------------
arma::rowvec gaussian_curve(const arma::vec &t, double tau, double sigma_w) {
  return arma::exp(-arma::pow(t - tau, 2) / (2.0 * sigma_w * sigma_w)).t();
}

// -----------------------------------------------------------------------------
// HELPER: Quadratic Interpolation (Brent's Method)
// -----------------------------------------------------------------------------
double interpolate_peak(double left, double center, double right) {
  double numerator = left - right;
  double denominator = 2.0 * (left - 2.0 * center + right);
  if (std::abs(denominator) < 1e-9)
    return 0.0;
  return 0.5 * (numerator / denominator);
}

// -----------------------------------------------------------------------------
// MAIN ENGINE: Vectorized Template Scan (Direct X_clr block processing)
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List scan_sequences(arma::mat X_clr, arma::vec time_vec, arma::vec tau_grid,
                    arma::vec epsilon_grid, double min_score = 0.5) {

  // X_clr must be oriented as Genes x Time
  int n_genes = X_clr.n_rows;
  int n_time = X_clr.n_cols;

  // Total archetypes = tau_grid size * epsilon_grid size * 2 structures
  int n_tau = tau_grid.n_elem;
  int n_eps = epsilon_grid.n_elem;
  int n_archetypes = n_tau * n_eps * 2;

  arma::mat M(n_archetypes, n_time);
  arma::vec meta_tau(n_archetypes);
  arma::vec meta_eps(n_archetypes);
  std::vector<int> meta_type(n_archetypes); // 1 = Logistic, 2 = Gaussian

  int idx = 0;
  for (int i = 0; i < n_tau; ++i) {
    for (int j = 0; j < n_eps; ++j) {
      // Archetype A: Logistic Switch
      M.row(idx) = logistic_curve(time_vec, tau_grid(i), epsilon_grid(j));
      meta_tau(idx) = tau_grid(i);
      meta_eps(idx) = epsilon_grid(j);
      meta_type[idx] = 1;
      idx++;

      // Archetype B: Gaussian Pulse
      M.row(idx) = gaussian_curve(time_vec, tau_grid(i), epsilon_grid(j));
      meta_tau(idx) = tau_grid(i);
      meta_eps(idx) = epsilon_grid(j);
      meta_type[idx] = 2;
      idx++;
    }
  }

  // Row-wise normalize M dictionary for strict Pearson correlation translation
  for (int i = 0; i < M.n_rows; ++i) {
    arma::rowvec r = M.row(i);
    double mu = arma::mean(r);
    double sigma = arma::stddev(r);
    if (sigma > 1e-12)
      M.row(i) = (r - mu) / sigma;
  }

  // Setup thread-safe accumulator containers
  int n_threads = 1;
#ifdef _OPENMP
  n_threads = omp_get_max_threads();
#endif

  std::vector<std::vector<int>> res_gene1(n_threads);
  std::vector<std::vector<int>> res_gene2(n_threads);
  std::vector<std::vector<double>> res_tau(n_threads);
  std::vector<std::vector<double>> res_eps(n_threads);
  std::vector<std::vector<double>> res_score(n_threads);
  std::vector<std::vector<int>> res_type(n_threads);

  int update_interval =
      std::max(1, (n_genes - 1) / 40); // 2.5% progression ticks
  int progress_counter = 0;

  // Execute Combinatoric Map-Reduce Pattern
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < n_genes - 1; ++i) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif

    // Track thread-safe progress efficiently on outer loop limits
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      progress_counter++;
      if (progress_counter % update_interval == 0) {
        Rcpp::Rcout << ".";
      }
    }

    for (int j = i + 1; j < n_genes; ++j) {
      // Generate exact pairwise difference dynamically
      arma::rowvec y = X_clr.row(i) - X_clr.row(j);

      // Inline CoDa Robustness
      winsorize_and_zscore(y);

      // Strict Pearson Correlation Vector computation
      arma::vec cor = (M * y.t()) / (n_time - 1);

      // Discover matching optimal topology
      uword best_idx;
      double best_val = cor.max(best_idx);

      // Core Threshold Pruning Logic
      // If the signal matches purely random noise (low R), safely obliterate it
      // from RAM
      if (best_val >= min_score) {
        double t_grid = meta_tau(best_idx);
        double e_grid = meta_eps(best_idx);
        int m_type = meta_type[best_idx];

        double t_final = t_grid;

        // Continuous boundary layer interpolation
        int stride = 2 * n_eps; // Steps between similar models at adjacent Taus
        if (best_idx >= stride && best_idx < (n_archetypes - stride)) {
          double val_left = cor(best_idx - stride);
          double val_right = cor(best_idx + stride);
          double step = tau_grid(1) - tau_grid(0);
          double offset = interpolate_peak(val_left, best_val, val_right);
          t_final = t_grid + (offset * step);
        }

        // Retain valid event parameters
        res_gene1[tid].push_back(i + 1); // Expose 1-indexed to R
        res_gene2[tid].push_back(j + 1);
        res_tau[tid].push_back(t_final);
        res_eps[tid].push_back(e_grid);
        res_score[tid].push_back(best_val);
        res_type[tid].push_back(m_type);
      }
    }
  }

  // Reduce logic: linearizing memory
  std::vector<int> out_gene1;
  std::vector<int> out_gene2;
  std::vector<double> out_tau;
  std::vector<double> out_eps;
  std::vector<double> out_score;
  std::vector<int> out_type;

  for (int t = 0; t < n_threads; ++t) {
    out_gene1.insert(out_gene1.end(), res_gene1[t].begin(), res_gene1[t].end());
    out_gene2.insert(out_gene2.end(), res_gene2[t].begin(), res_gene2[t].end());
    out_tau.insert(out_tau.end(), res_tau[t].begin(), res_tau[t].end());
    out_eps.insert(out_eps.end(), res_eps[t].begin(), res_eps[t].end());
    out_score.insert(out_score.end(), res_score[t].begin(), res_score[t].end());
    out_type.insert(out_type.end(), res_type[t].begin(), res_type[t].end());
  }

  return List::create(Named("GeneA_Idx") = out_gene1,
                      Named("GeneB_Idx") = out_gene2, Named("Tau") = out_tau,
                      Named("Epsilon") = out_eps, Named("Score") = out_score,
                      Named("Model_Type") = out_type);
}