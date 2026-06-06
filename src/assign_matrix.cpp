// src/assign_matrix.cpp

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void scatter_assign_mat(List geno_list, IntegerVector rows, IntegerVector cols, List vals) {

  int n = rows.size();
  int n_samples = geno_list.size();

  // --- cache matrices (same idea as vec version) ---
  std::vector<CharacterMatrix> cache(n_samples);

  for (int i = 0; i < n_samples; i++) {
    cache[i] = as<CharacterMatrix>(geno_list[i]);
  }

  // --- scatter ---
  for (int i = 0; i < n; i++) {

    int r = rows[i] - 1;  // sample index
    int c = cols[i] - 1;  // locus index

    CharacterMatrix& mat = cache[r];
    CharacterVector v = vals[i];  // length = n_rows_per_cell

    int kmax = v.size();

    // --- fast path for diploid (common case) ---
    if (kmax == 2) {
      mat(0, c) = v[0];
      mat(1, c) = v[1];
    } else {
      for (int k = 0; k < kmax; k++) {
        mat(k, c) = v[k];
      }
    }
  }
}
