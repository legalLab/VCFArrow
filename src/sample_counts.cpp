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
#include <climits>
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


// ── Inbreeding coefficient F (vcf_stats()) ──────────────────────────────────
//
// Method-of-moments F per sample, as vcftools --het and plink --het:
//   F = (O - E) / (N - E)
// over the biallelic SNVs (ok[vpos]) at which the sample is called, at least
// 2 samples of its reference set are called (n >= 4 called alleles) and the
// reference set is polymorphic (0 < p < 1; a monomorphic site would add 1 to
// N, O and E alike, leaving F unchanged; vcftools skips them too):
//   N += 1;  O += 1 if homozygous;  E += 1 - n / (n - 1) * 2p(1 - p)
// where p is the ALT allele frequency and n the number of called alleles at
// the site in the reference set: all samples (n_groups = 1) or the sample's
// group (grp, 1-based; samples with grp < 1 are skipped).
//
// Inputs (rows of one chunk): s = 1-based sample index, vpos = 1-based
// position of the row's variant in the variant table, a1/a2 = allele codes
// (0 = REF, 1 = ALT, NA = missing).

// Allele counts per (group, variant) for the variants of one chunk, over the
// range of positions it covers: element (vpos - lo) * n_groups + (g - 1).
// [[Rcpp::export]]
List f_allele_counts_cpp(const IntegerVector& s, const IntegerVector& vpos,
                         const IntegerVector& a1, const IntegerVector& a2,
                         const LogicalVector& ok, const IntegerVector& grp,
                         int n_groups) {
  const R_xlen_t n = s.size();
  if (vpos.size() != n || a1.size() != n || a2.size() != n)
    stop("f_allele_counts_cpp: input lengths differ");
  int lo = INT_MAX, hi = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const int v = vpos[i];
    if (v == NA_INTEGER) continue;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  const int nv = hi >= lo ? hi - lo + 1 : 0;
  IntegerVector alt((R_xlen_t)nv * n_groups), nobs((R_xlen_t)nv * n_groups);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int v = vpos[i], si = s[i], x = a1[i], y = a2[i];
    if (v == NA_INTEGER || si == NA_INTEGER || x == NA_INTEGER || y == NA_INTEGER) continue;
    if (v > ok.size() || ok[v - 1] != TRUE) continue;
    const int g = n_groups == 1 ? 1 : grp[si - 1];
    if (g < 1 || g > n_groups) continue;
    const R_xlen_t j = (R_xlen_t)(v - lo) * n_groups + (g - 1);
    alt[j] += x + y;
    nobs[j] += 2;
  }
  return List::create(_["lo"] = nv ? lo : 1, _["n_var"] = nv, _["alt"] = alt, _["nobs"] = nobs);
}

// Per-sample N, O and E for the rows of one chunk, given allele counts per
// (group, variant) laid out as in f_allele_counts_cpp() from position lo on
// (the chunk's own counts, or counts over all chunks with lo = 1).
// [[Rcpp::export]]
List f_hom_cpp(const IntegerVector& s, const IntegerVector& vpos,
               const IntegerVector& a1, const IntegerVector& a2,
               const LogicalVector& ok, const IntegerVector& grp,
               const IntegerVector& alt, const IntegerVector& nobs, int lo,
               int n_groups, int n_samples) {
  const R_xlen_t n = s.size(), m = alt.size();
  if (vpos.size() != n || a1.size() != n || a2.size() != n || nobs.size() != m)
    stop("f_hom_cpp: input lengths differ");
  IntegerVector n_sites(n_samples), obs_hom(n_samples);
  NumericVector exp_hom(n_samples);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int v = vpos[i], si = s[i], x = a1[i], y = a2[i];
    if (v == NA_INTEGER || si == NA_INTEGER || x == NA_INTEGER || y == NA_INTEGER) continue;
    if (si < 1 || si > n_samples || v > ok.size() || ok[v - 1] != TRUE) continue;
    const int g = n_groups == 1 ? 1 : grp[si - 1];
    if (g < 1 || g > n_groups) continue;
    const R_xlen_t j = (R_xlen_t)(v - lo) * n_groups + (g - 1);
    if (j < 0 || j >= m) stop("f_hom_cpp: variant outside the allele counts");
    const int nn = nobs[j];
    if (nn < 4) continue;                      // fewer than 2 called samples
    if (alt[j] == 0 || alt[j] == nn) continue; // monomorphic in the reference set
    const double p = static_cast<double>(alt[j]) / nn;
    const int k = si - 1;
    ++n_sites[k];
    if (x == y) ++obs_hom[k];
    exp_hom[k] += 1.0 - 2.0 * p * (1.0 - p) * nn / (nn - 1.0);
  }
  return List::create(_["n_sites"] = n_sites, _["obs_hom"] = obs_hom, _["exp_hom"] = exp_hom);
}
