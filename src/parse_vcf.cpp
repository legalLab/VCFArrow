// src/parse_vcf.cpp

#include <Rcpp.h>
using namespace Rcpp;

// ---- fast double parser (much faster than atof) ----
inline double fast_atof(const char* p, const char* end) {
  if (p >= end || *p == '.') return NA_REAL;

  double val = 0.0;
  while (p < end && *p >= '0' && *p <= '9') {
    val = val * 10.0 + (*p - '0');
    ++p;
  }

  if (p < end && *p == '.') {
    ++p;
    double frac = 0.1;
    while (p < end && *p >= '0' && *p <= '9') {
      val += (*p - '0') * frac;
      frac *= 0.1;
      ++p;
    }
  }

  return val;
}

// ---- compare token without allocating ----
inline bool token_eq(const char* start, const char* end, const char* str) {
  size_t len = end - start;
  size_t i = 0;
  for (; i < len && str[i]; i++) {
    if (start[i] != str[i]) return false;
  }
  return (i == len && str[i] == '\0');
}

// [[Rcpp::export]]
List parse_vcf_cpp(CharacterVector lines, int nsamples) {

  int n = lines.size();

  IntegerMatrix a1(n, nsamples);
  IntegerMatrix a2(n, nsamples);
  LogicalMatrix phased(n, nsamples);

  NumericMatrix DP(n, nsamples);
  NumericMatrix GQ(n, nsamples);

  CharacterMatrix variants(n, 7);
  CharacterVector info(n);
  CharacterVector format(n);
  CharacterMatrix fmt(n, nsamples);

  for (int i = 0; i < n; i++) {

    const char* ptr = CHAR(lines[i]);
    const char* field_start = ptr;
    const char* field_end;

    // ---- parse first 9 fields ----
    const char* fields[9];

    for (int k = 0; k < 9; k++) {
      field_start = ptr;
      while (*ptr && *ptr != '\t') ptr++;
      fields[k] = field_start;

      if (*ptr == '\t') {
        *const_cast<char*>(ptr) = '\0';  // safe: R makes a copy-on-write buffer
        ptr++;
      }
    }

    // store fixed columns
    for (int k = 0; k < 7; k++) variants(i,k) = fields[k];
    info[i] = fields[7];
    format[i] = fields[8];

    // ---- parse FORMAT keys (pointer-based) ----
    int dp_pos = -1, gq_pos = -1, hq_pos = -1;

    int key_idx = 0;
    const char* p = fields[8];
    const char* key_start = p;

    while (true) {
      if (*p == ':' || *p == '\0') {
        if (token_eq(key_start, p, "DP")) dp_pos = key_idx;
        else if (token_eq(key_start, p, "GQ")) gq_pos = key_idx;
        else if (token_eq(key_start, p, "HQ")) hq_pos = key_idx;

        key_idx++;
        if (*p == '\0') break;
        key_start = p + 1;
      }
      p++;
    }

    // ---- samples ----
    for (int j = 0; j < nsamples; j++) {

      field_start = ptr;
      while (*ptr && *ptr != '\t') ptr++;
      field_end = ptr;

      fmt(i,j) = std::string(field_start, field_end);

      // ---- parse FORMAT values inline ----
      int val_idx = 0;
      const char* v = field_start;
      const char* val_start = v;

      double dp_val = NA_REAL;
      double gq_val = NA_REAL;

      // GT parsing
      if (field_end - field_start >= 3) {
        char c1 = field_start[0];
        char sep = field_start[1];
        char c2 = field_start[2];

        phased(i,j) = (sep == '|');

        if (c1 == '.' || c2 == '.') {
          a1(i,j) = NA_INTEGER;
          a2(i,j) = NA_INTEGER;
        } else {
          a1(i,j) = c1 - '0';
          a2(i,j) = c2 - '0';
        }
      } else {
        a1(i,j) = NA_INTEGER;
        a2(i,j) = NA_INTEGER;
        phased(i,j) = false;
      }

      while (true) {

        if (v == field_end || *v == ':') {

          if (val_idx == dp_pos) {
            dp_val = fast_atof(val_start, v);
          }

          if (val_idx == gq_pos) {
            gq_val = fast_atof(val_start, v);
          }

          if (v == field_end) break;

          val_idx++;
          val_start = v + 1;
        }
        v++;
      }

      DP(i,j) = dp_val;
      GQ(i,j) = gq_val;

      if (*ptr == '\t') ptr++;
    }
  }

  return List::create(
    Named("variants") = variants,
    Named("info") = info,
    Named("format") = format,
    Named("fmt") = fmt,
    Named("a1") = a1,
    Named("a2") = a2,
    Named("phased") = phased,
    Named("DP") = DP,
    Named("GQ") = GQ
  );
}
