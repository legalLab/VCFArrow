// src/assign_vector.cpp

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void scatter_assign_vec(List geno_list, IntegerVector rows, IntegerVector cols, CharacterVector vals) {

  int n = rows.size();

  std::vector<CharacterVector> cache(geno_list.size());

  for (int i = 0; i < geno_list.size(); i++) {
    cache[i] = geno_list[i];
  }

  for (int i = 0; i < n; i++) {
    int r = rows[i] - 1;
    int c = cols[i] - 1;

    cache[r][c] = vals[i];
  }
}
