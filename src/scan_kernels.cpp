// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

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
// HELPER: Bicor Standardization (L2 Normalized robust weights)
// -----------------------------------------------------------------------------
inline arma::rowvec compute_bicor_zscore(const arma::rowvec& y) {
    int n = y.n_elem;
    
    // 1. Calculate Median and MAD
    double med = arma::median(y);
    arma::rowvec abs_dev = arma::abs(y - med);
    double mad = arma::median(abs_dev);
    
    if (mad == 0.0) mad = 1e-8; 
    
    // 2. Calculate the Biweights
    arma::rowvec u = (y - med) / (9.0 * mad);
    arma::rowvec weights = arma::zeros<arma::rowvec>(n);
    
    for(int i = 0; i < n; ++i) {
        if (std::abs(u[i]) < 1.0) {
            double one_minus_u2 = 1.0 - (u[i] * u[i]);
            weights[i] = one_minus_u2 * one_minus_u2; // (1 - u^2)^2
        } else {
            weights[i] = 0.0;
        }
    }
    
    // 3. Transform to Topological Projection
    arma::rowvec y_centered_weighted = (y - med) % weights;
    
    // 4. L2 Normalization (incorporates n-1 natively)
    double sum_sq = arma::sum(arma::square(y_centered_weighted));
    
    if (sum_sq < 1e-12) return arma::zeros<arma::rowvec>(n);
    
    return y_centered_weighted / std::sqrt(sum_sq);
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
// HELPER: RcppParallel Worker
// -----------------------------------------------------------------------------
struct ScanWorker : public RcppParallel::Worker {
  const arma::mat &X_clr;
  const arma::mat &M;
  const arma::vec &time_vec_z;
  const arma::vec &tau_grid;
  const arma::vec &epsilon_grid;
  const arma::vec &meta_tau;
  const arma::vec &meta_eps;
  const std::vector<int> &meta_type;
  double min_score;
  double min_var_delta;
  int method;
  int n_genes;
  int n_time;
  int n_archetypes;

  // Lock-free result containers partitioned by GeneA index `i`
  std::vector<std::vector<int>> &res_gene1;
  std::vector<std::vector<int>> &res_gene2;
  std::vector<std::vector<double>> &res_tau;
  std::vector<std::vector<double>> &res_eps;
  std::vector<std::vector<double>> &res_r2_lin;
  std::vector<std::vector<double>> &res_r2_sig;
  std::vector<std::vector<double>> &res_tau_dec;
  std::vector<std::vector<double>> &res_var_delta;
  std::vector<std::vector<int>> &res_type;

  ScanWorker(const arma::mat &X_clr, const arma::mat &M, const arma::vec &time_vec_z,
             const arma::vec &tau_grid, const arma::vec &epsilon_grid,
             const arma::vec &meta_tau, const arma::vec &meta_eps,
             const std::vector<int> &meta_type, double min_score, double min_var_delta,
             int method, int n_genes, int n_time, int n_archetypes,
             std::vector<std::vector<int>> &res_gene1,
             std::vector<std::vector<int>> &res_gene2,
             std::vector<std::vector<double>> &res_tau,
             std::vector<std::vector<double>> &res_eps,
             std::vector<std::vector<double>> &res_r2_lin,
             std::vector<std::vector<double>> &res_r2_sig,
             std::vector<std::vector<double>> &res_tau_dec,
             std::vector<std::vector<double>> &res_var_delta,
             std::vector<std::vector<int>> &res_type)
      : X_clr(X_clr), M(M), time_vec_z(time_vec_z), tau_grid(tau_grid), epsilon_grid(epsilon_grid),
        meta_tau(meta_tau), meta_eps(meta_eps), meta_type(meta_type),
        min_score(min_score), min_var_delta(min_var_delta), method(method), n_genes(n_genes), 
        n_time(n_time), n_archetypes(n_archetypes), res_gene1(res_gene1), 
        res_gene2(res_gene2), res_tau(res_tau), res_eps(res_eps), 
        res_r2_lin(res_r2_lin), res_r2_sig(res_r2_sig), res_tau_dec(res_tau_dec), 
        res_var_delta(res_var_delta), res_type(res_type) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      for (int j = i + 1; j < n_genes; ++j) {
        arma::rowvec y_raw = X_clr.row(i) - X_clr.row(j);
        
        // Z-score strictly isolated for Var_Delta decoupling physics
        arma::rowvec y_z = y_raw;
        winsorize_and_zscore(y_z);
        
        arma::rowvec y_processed;
        if (method == 1) {
            y_processed = compute_bicor_zscore(y_raw);
        } else {
            y_processed = y_z / std::sqrt(n_time - 1.0); // L2 transform for dot product
        }

        // 1. Phase Transition fit (Sigmoid) via strict L2 dot product correlation
        arma::vec cor = (M * y_processed.t());
        arma::uword best_idx = 0;
        double best_val = -1.0;
        for (arma::uword k = 0; k < cor.n_elem; ++k) {
          double a = std::abs(cor[k]);
          if (a > best_val) {
            best_val = a;
            best_idx = k;
          }
        }
        double r2_sigmoid = best_val * best_val;
        
        // 2. Linear Drift fit (Line)
        double r_lin = arma::dot(y_processed, time_vec_z);
        double r2_linear = r_lin * r_lin;

        // 3. O(1) Decoupling Math (Variance Explosion/Collapse)
        // Strictly uses raw, continuous variance arrays (y_z) mapped chronologically
        std::vector<double> cumsum_y(n_time + 1, 0.0);
        std::vector<double> cumsum_y2(n_time + 1, 0.0);
        for (int k = 0; k < n_time; ++k) {
            cumsum_y[k+1] = cumsum_y[k] + y_z[k];
            cumsum_y2[k+1] = cumsum_y2[k] + y_z[k] * y_z[k];
        }
        
        double max_var_delta = 0.0;
        double best_tau_dec = 0.0;
        
        for (int t = 0; t < tau_grid.n_elem; ++t) {
            double tau_val = tau_grid(t);
            int split_k = std::max(2, std::min(n_time - 2, (int)std::round(tau_val * n_time)));
            
            // Safe Variance Pre
            double sum1 = cumsum_y[split_k];
            double sum_sq1 = cumsum_y2[split_k];
            double var_pre = (sum_sq1 - (sum1 * sum1) / split_k) / (split_k - 1.0);
            if (var_pre < 1e-6) var_pre = 1e-6; 
            
            // Safe Variance Post
            int n_post = n_time - split_k;
            double sum2 = cumsum_y[n_time] - sum1;
            double sum_sq2 = cumsum_y2[n_time] - sum_sq1;
            double var_post = (sum_sq2 - (sum2 * sum2) / n_post) / (n_post - 1.0);
            if (var_post < 1e-6) var_post = 1e-6; 
            
            double ratio = std::max(var_pre / var_post, var_post / var_pre);
            if (ratio > max_var_delta) {
                max_var_delta = ratio;
                best_tau_dec = tau_val;
            }
        }

        // 4. Sparse Boolean Bottleneck
        bool keep = (r2_sigmoid >= min_score || r2_linear >= min_score || max_var_delta >= min_var_delta);

        if (keep) {
          double t_grid = meta_tau(best_idx);
          double e_grid = meta_eps(best_idx);
          int m_type = meta_type[best_idx];
          double t_final = t_grid;

          int stride = 2 * epsilon_grid.n_elem;
          if (best_idx >= stride && best_idx < (n_archetypes - stride)) {
            double val_left = cor(best_idx - stride);
            double val_right = cor(best_idx + stride);
            double step = tau_grid(1) - tau_grid(0);

            double numerator = std::abs(val_left) - std::abs(val_right);
            double denominator = 2.0 * (std::abs(val_left) - 2.0 * std::abs(cor(best_idx)) + std::abs(val_right));
            double offset = 0.0;
            if (std::abs(denominator) > 1e-9) {
              offset = 0.5 * (numerator / denominator);
            }
            t_final = t_grid + (offset * step);
          }

          res_gene1[i].push_back(static_cast<int>(i) + 1);
          res_gene2[i].push_back(j + 1);
          res_tau[i].push_back(t_final);
          res_eps[i].push_back(e_grid);
          res_r2_sig[i].push_back(r2_sigmoid);
          res_r2_lin[i].push_back(r2_linear);
          res_var_delta[i].push_back(max_var_delta);
          res_tau_dec[i].push_back(best_tau_dec);
          res_type[i].push_back(m_type);
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
// MAIN ENGINE: Vectorized Template Scan (Direct X_clr block processing)
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List scan_sequences(arma::mat X_clr, arma::vec time_vec, arma::vec tau_grid,
                    arma::vec epsilon_grid, double min_score = 0.5, double min_var_delta = 5.0,
                    int method = 0) {

  int n_genes = X_clr.n_rows;
  int n_time = X_clr.n_cols;

  // Standardization of the linear target explicitly locked to regression method
  arma::vec time_vec_z;
  if (method == 1) {
      time_vec_z = compute_bicor_zscore(time_vec.t()).t();
  } else {
      time_vec_z = time_vec;
      double t_mu = arma::mean(time_vec_z);
      double t_sig = arma::stddev(time_vec_z);
      if (t_sig > 1e-12) {
          time_vec_z = (time_vec_z - t_mu) / t_sig;
      } else {
          time_vec_z.zeros();
      }
      time_vec_z = time_vec_z / std::sqrt(n_time - 1.0); // Standardize to L2 Pearson
  }

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
      M.row(idx) = logistic_curve(time_vec, tau_grid(i), epsilon_grid(j));
      meta_tau(idx) = tau_grid(i);
      meta_eps(idx) = epsilon_grid(j);
      meta_type[idx] = 1;
      idx++;

      M.row(idx) = gaussian_curve(time_vec, tau_grid(i), epsilon_grid(j));
      meta_tau(idx) = tau_grid(i);
      meta_eps(idx) = epsilon_grid(j);
      meta_type[idx] = 2;
      idx++;
    }
  }

  for (int i = 0; i < M.n_rows; ++i) {
    if (method == 1) {
        M.row(i) = compute_bicor_zscore(M.row(i));
    } else {
        arma::rowvec r = M.row(i);
        double mu = arma::mean(r);
        double sigma = arma::stddev(r);
        if (sigma > 1e-12) {
          M.row(i) = (r - mu) / sigma;
        }
        M.row(i) = M.row(i) / std::sqrt(n_time - 1.0); // Standardize to L2 Pearson
    }
  }

  std::vector<std::vector<int>> res_gene1(n_genes);
  std::vector<std::vector<int>> res_gene2(n_genes);
  std::vector<std::vector<double>> res_tau(n_genes);
  std::vector<std::vector<double>> res_eps(n_genes);
  
  std::vector<std::vector<double>> res_r2_lin(n_genes);
  std::vector<std::vector<double>> res_r2_sig(n_genes);
  std::vector<std::vector<double>> res_tau_dec(n_genes);
  std::vector<std::vector<double>> res_var_delta(n_genes);
  
  std::vector<std::vector<int>> res_type(n_genes);

  ScanWorker worker(X_clr, M, time_vec_z, tau_grid, epsilon_grid, meta_tau, meta_eps,
                    meta_type, min_score, min_var_delta, method, n_genes, n_time, n_archetypes,
                    res_gene1, res_gene2, res_tau, res_eps, 
                    res_r2_lin, res_r2_sig, res_tau_dec, res_var_delta, res_type);

  RcppParallel::parallelFor(0, n_genes - 1, worker);

  int total_size = 0;
  for (int i = 0; i < n_genes; ++i) {
    total_size += res_gene1[i].size();
  }

  std::vector<int> out_gene1, out_gene2, out_type;
  std::vector<double> out_tau, out_eps, out_r2_lin, out_r2_sig, out_tau_dec, out_var_delta;

  out_gene1.reserve(total_size);
  out_gene2.reserve(total_size);
  out_tau.reserve(total_size);
  out_eps.reserve(total_size);
  out_r2_lin.reserve(total_size);
  out_r2_sig.reserve(total_size);
  out_tau_dec.reserve(total_size);
  out_var_delta.reserve(total_size);
  out_type.reserve(total_size);

  for (int i = 0; i < n_genes; ++i) {
    if (!res_gene1[i].empty()) {
      out_gene1.insert(out_gene1.end(), res_gene1[i].begin(), res_gene1[i].end());
      out_gene2.insert(out_gene2.end(), res_gene2[i].begin(), res_gene2[i].end());
      out_tau.insert(out_tau.end(), res_tau[i].begin(), res_tau[i].end());
      out_eps.insert(out_eps.end(), res_eps[i].begin(), res_eps[i].end());
      out_r2_lin.insert(out_r2_lin.end(), res_r2_lin[i].begin(), res_r2_lin[i].end());
      out_r2_sig.insert(out_r2_sig.end(), res_r2_sig[i].begin(), res_r2_sig[i].end());
      out_tau_dec.insert(out_tau_dec.end(), res_tau_dec[i].begin(), res_tau_dec[i].end());
      out_var_delta.insert(out_var_delta.end(), res_var_delta[i].begin(), res_var_delta[i].end());
      out_type.insert(out_type.end(), res_type[i].begin(), res_type[i].end());
    }
  }

  return List::create(Named("GeneA_Idx") = out_gene1,
                      Named("GeneB_Idx") = out_gene2, 
                      Named("Tau") = out_tau,
                      Named("Epsilon") = out_eps, 
                      Named("R2_Linear") = out_r2_lin,
                      Named("R2_Sigmoid") = out_r2_sig,
                      Named("Tau_Decouple") = out_tau_dec,
                      Named("Var_Delta") = out_var_delta,
                      Named("Model_Type") = out_type);
}