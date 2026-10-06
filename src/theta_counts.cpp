// src/theta_counts.cpp
//
// Watterson's theta and pi sums for one chunk, in a single pass over the
// chunk's genotype rows (used by vcf_theta()).  For each variant and group the
// ALT alleles and observed alleles are counted; then, for each (group,
// variant) and for each variant over all groups, with nobs observed alleles:
//   valid   = nobs / 2 >= 2                       (at least 2 called samples)
//   p       = alt / nobs
//   pi      = 2 * p * (1 - p)
//   contrib = 1 / H(nobs - 1) if 0 < p < 1        (Watterson's a1 for nobs
//             chromosomes; H(k) = sum of 1/i for i = 1..k)
// and the sums of contrib, pi and valid are returned.  Same definitions as
// .theta_sums() in R/vcf_theta.R.

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

static inline void add_site(int alt, int nobs, const NumericVector& harmonic,
                            double* contrib, double* pi, double* n_valid) {
  if (nobs < 4) return;                      // n_called = nobs / 2 < 2
  const double p = static_cast<double>(alt) / nobs;
  *pi += 2.0 * p * (1.0 - p);
  *n_valid += 1.0;
  if (p > 0.0 && p < 1.0) {
    if (nobs - 1 > harmonic.size()) stop("theta_chunk_cpp: harmonic numbers too short");
    *contrib += 1.0 / harmonic[nobs - 2];     // H(nobs - 1), 1-based
  }
}

// var: 1-based variant index within the chunk of each row (1..n_var)
// grp: 1-based group of each row's sample (rows with grp 0 or NA are skipped)
// a1, a2: allele codes (0 = REF, 1 = ALT, NA = missing)
// harmonic: H(1), H(2), ..., at least up to H(2 * number of samples - 1)
// [[Rcpp::export]]
List theta_chunk_cpp(const IntegerVector& var, const IntegerVector& grp,
                     const IntegerVector& a1, const IntegerVector& a2,
                     int n_var, int n_pops, const NumericVector& harmonic) {
  const R_xlen_t n = var.size();
  if (grp.size() != n || a1.size() != n || a2.size() != n)
    stop("theta_chunk_cpp: input lengths differ");

  // counts per (group, variant), group-major within each variant
  std::vector<int> alt(static_cast<size_t>(n_pops) * n_var, 0),
    nobs(static_cast<size_t>(n_pops) * n_var, 0);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int v = var[i], g = grp[i];
    if (v == NA_INTEGER || v < 1 || v > n_var) continue;
    if (g == NA_INTEGER || g < 1 || g > n_pops) continue;
    const size_t k = static_cast<size_t>(v - 1) * n_pops + (g - 1);
    const int x = a1[i], y = a2[i];
    if (x != NA_INTEGER) { ++nobs[k]; if (x == 1) ++alt[k]; }
    if (y != NA_INTEGER) { ++nobs[k]; if (y == 1) ++alt[k]; }
  }

  NumericMatrix group(n_pops, 3);            // contrib, pi, n_valid per group
  double t_contrib = 0, t_pi = 0, t_valid = 0, nobs_sum = 0;
  for (int v = 0; v < n_var; ++v) {
    int alt_t = 0, nobs_t = 0;
    for (int g = 0; g < n_pops; ++g) {
      const size_t k = static_cast<size_t>(v) * n_pops + g;
      add_site(alt[k], nobs[k], harmonic, &group(g, 0), &group(g, 1), &group(g, 2));
      alt_t += alt[k];
      nobs_t += nobs[k];
    }
    add_site(alt_t, nobs_t, harmonic, &t_contrib, &t_pi, &t_valid);
    nobs_sum += nobs_t;
  }
  colnames(group) = CharacterVector::create("contrib", "pi", "n_valid");
  return List::create(
    _["group"] = group,
    _["total"] = NumericVector::create(_["contrib"] = t_contrib, _["pi"] = t_pi,
                                       _["n_valid"] = t_valid),
    _["nobs"] = nobs_sum);
}
