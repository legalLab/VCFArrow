// src/sample_counts.cpp
//
// Per-sample genotype counts for one chunk, in a single pass (used by
// vcf_stats() and .scan_vcf_gt()).  Equivalent to tabulating these R masks
// by sample:
//   called  = !(is.na(a1) | is.na(a2))
//   het     = called & a1 != a2
//   hom_ref = called & a1 == 0 & a2 == 0
//   hom_alt = called & a1 == 1 & a2 == 1
//   dp_ok   = !is.na(DP)             (dp_n, and dp_sum = sum of DP)

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// s: 1-based sample index of each row (rows with NA s are skipped)
// dp: may be empty (no DP counts)
// [[Rcpp::export]]
List sample_counts_cpp(const IntegerVector& s, const IntegerVector& a1,
                       const IntegerVector& a2, const NumericVector& dp,
                       int n_samples) {
  const R_xlen_t n = s.size();
  const bool use_dp = dp.size() > 0;
  if (a1.size() != n || a2.size() != n || (use_dp && dp.size() != n))
    stop("sample_counts_cpp: input lengths differ");

  IntegerVector total(n_samples), called(n_samples), het(n_samples),
    hom_ref(n_samples), hom_alt(n_samples), dp_n(n_samples);
  NumericVector dp_sum(n_samples);

  for (R_xlen_t i = 0; i < n; ++i) {
    const int si = s[i];
    if (si == NA_INTEGER || si < 1 || si > n_samples) continue;
    const int k = si - 1;
    ++total[k];
    const int x = a1[i], y = a2[i];
    if (x != NA_INTEGER && y != NA_INTEGER) {
      ++called[k];
      if (x != y) ++het[k];
      else if (x == 0) ++hom_ref[k];
      else if (x == 1) ++hom_alt[k];
    }
    if (use_dp && !std::isnan(dp[i])) {
      ++dp_n[k];
      dp_sum[k] += dp[i];
    }
  }

  return List::create(
    Named("total") = total, Named("called") = called, Named("het") = het,
    Named("hom_ref") = hom_ref, Named("hom_alt") = hom_alt,
    Named("dp_n") = dp_n, Named("dp_sum") = dp_sum
  );
}


// match(ids, row_ids) given pos = .row_id_pos(row_ids): pos[id], or NA when
// id is NA, out of range, or absent (pos == 0).  One pass, no temporaries.
// [[Rcpp::export]]
IntegerVector match_row_id_cpp(const IntegerVector& ids, const IntegerVector& pos) {
  const R_xlen_t n = ids.size(), m = pos.size();
  IntegerVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int id = ids[i];
    int p = NA_INTEGER;
    if (id != NA_INTEGER && id >= 1 && id <= m) {
      p = pos[id - 1];
      if (p == 0) p = NA_INTEGER;
    }
    out[i] = p;
  }
  return out;
}
