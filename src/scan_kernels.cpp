// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <omp.h>

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// HELPER: Row-wise Z-Score Normalization
// Ensures that Dot Product == Pearson Correlation
// -----------------------------------------------------------------------------
arma::mat normalize_rows(arma::mat X) {
    int n_rows = X.n_rows;
    int n_cols = X.n_cols;
    
    arma::mat X_norm = X;
    
    // Parallelize row normalization
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_rows; ++i) {
        arma::rowvec r = X.row(i);
        double mu = arma::mean(r);
        double sigma = arma::stddev(r);
        
        if (sigma > 1e-12) {
            X_norm.row(i) = (r - mu) / sigma;
        } else {
            X_norm.row(i).zeros(); // Handle constant rows (variance = 0)
        }
    }
    return X_norm;
}

// -----------------------------------------------------------------------------
// HELPER: Generalized Logistic Function (The Boundary Function)
// -----------------------------------------------------------------------------
arma::rowvec logistic_curve(arma::vec t, double tau, double epsilon) {
    // Formula: 1 / (1 + exp( -(t - tau)/epsilon ))
    // We center it to mean=0 later, so amplitude doesn't matter here (0 to 1 is fine)
    arma::rowvec y = 1.0 / (1.0 + arma::exp(-(t - tau) / epsilon)).t();
    return y;
}

// -----------------------------------------------------------------------------
// HELPER: Quadratic Interpolation (Brent's Method)
// Returns the offset from the grid point (-0.5 to +0.5 of a step)
// -----------------------------------------------------------------------------
double interpolate_peak(double left, double center, double right) {
    double numerator = left - right;
    double denominator = 2.0 * (left - 2.0 * center + right);
    
    if (std::abs(denominator) < 1e-9) return 0.0; // Flat peak or error
    return 0.5 * (numerator / denominator);
}

// -----------------------------------------------------------------------------
// MAIN ENGINE: Vectorized Template Scan
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List scan_sequences(arma::mat Y_raw, arma::vec time_vec, 
                   arma::vec tau_grid, arma::vec epsilon_grid) {
    
    // 1. DYNAMIC DICTIONARY GENERATION ----------------------------------------
    // Total archetypes = tau_grid size * epsilon_grid size
    int n_tau = tau_grid.n_elem;
    int n_eps = epsilon_grid.n_elem;
    int n_archetypes = n_tau * n_eps;
    int n_time = time_vec.n_elem;
    
    arma::mat M(n_archetypes, n_time);
    arma::vec meta_tau(n_archetypes);
    arma::vec meta_eps(n_archetypes);
    
    int idx = 0;
    for (int i = 0; i < n_tau; ++i) {
        for (int j = 0; j < n_eps; ++j) {
            M.row(idx) = logistic_curve(time_vec, tau_grid(i), epsilon_grid(j));
            meta_tau(idx) = tau_grid(i);
            meta_eps(idx) = epsilon_grid(j);
            idx++;
        }
    }
    
    // 2. NORMALIZATION (The "Pearson" Step) -----------------------------------
    // Crucial: Z-score both Data (Y) and Dictionary (M) so Dot Prod = Correlation
    arma::mat Y = normalize_rows(Y_raw);
    arma::mat M_norm = normalize_rows(M);
    
    // 3. VECTORIZED CONVOLUTION (BLAS Level 3) --------------------------------
    // R (Pairs x Archetypes) = Y (Pairs x Time) * M_norm.t (Time x Archetypes)
    arma::mat Scores = Y * M_norm.t();
    
    // 4. ARGMAX & INTERPOLATION (Parallelized) --------------------------------
    int n_pairs = Y.n_rows;
    
    // Output containers
    arma::vec best_tau(n_pairs);
    arma::vec best_eps(n_pairs);
    arma::vec max_score(n_pairs);
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_pairs; ++i) {
        arma::rowvec s = Scores.row(i);
        
        // Find best discrete match
        uword best_idx;
        double best_val = s.max(best_idx);
        
        // Retrieve grid parameters
        double t_grid = meta_tau(best_idx);
        double e_grid = meta_eps(best_idx);
        
        // REFINEMENT: Quadratic Interpolation for Tau
        // We only interpolate if the best match is NOT at the edge of the grid
        // and if the grid structure allows it (neighbors must share same epsilon)
        // For simplicity, we assume dictionary is ordered by Tau then Epsilon
        
        double t_final = t_grid;
        
        // Check neighbors in the Tau dimension (stride = n_eps)
        // This assumes the loop structure: Tau (outer), Epsilon (inner)
        // If we are not at the first or last Tau
        if (best_idx >= n_eps && best_idx < (n_archetypes - n_eps)) {
            
            double val_left  = s(best_idx - n_eps); // Same epsilon, prev tau
            double val_right = s(best_idx + n_eps); // Same epsilon, next tau
            
            // Calculate step size (assuming uniform grid)
            double step = tau_grid(1) - tau_grid(0); 
            
            // Apply offset
            double offset = interpolate_peak(val_left, best_val, val_right);
            t_final = t_grid + (offset * step);
        }
        
        best_tau(i) = t_final;
        best_eps(i) = e_grid;
        max_score(i) = best_val;
    }
    
    return List::create(
        Named("tau") = best_tau,
        Named("epsilon") = best_eps,
        Named("score") = max_score
    );
}