// src/write_vcf_arrow.cpp

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <zlib.h>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

using namespace Rcpp;

// ---- buffer helpers ----
inline void buf_char(std::vector<char>& b, char c) {
  b.push_back(c);
}

inline void buf_cstr(std::vector<char>& b, const char* s) {
  if (!s) { b.push_back('.'); return; }
  while (*s) b.push_back(*s++);
}

inline void buf_int(std::vector<char>& b, int v) {
  char tmp[16];
  int len = snprintf(tmp, sizeof(tmp), "%d", v);
  b.insert(b.end(), tmp, tmp + len);
}

// writes one VCF field from an R character SEXP, substituting "." for NA/empty
inline void buf_rstr(std::vector<char>& b, SEXP sx) {
  if (sx == NA_STRING) { b.push_back('.'); return; }
  const char* s = CHAR(sx);
  if (!s || s[0] == '\0') { b.push_back('.'); return; }
  buf_cstr(b, s);
}

// [[Rcpp::export]]
void write_vcf_chunk_cpp(
    std::string output_file,
    CharacterVector chrom,
    IntegerVector pos,
    CharacterVector id,
    CharacterVector ref,
    CharacterVector alt,
    CharacterVector qual,
    CharacterVector filter_col,
    CharacterVector info,
    CharacterVector format_col,
    CharacterVector fmt_vec, // flat, row-major: variant i, sample j → fmt_vec[i*n_samples + j]
    int n_samples,
    bool gzip = false
    ) {

  int n_chroms = chrom.size();

  if (fmt_vec.size() != (R_xlen_t)(n_chroms * n_samples))
    stop("fmt_vec length must equal nrow * nsamples");

  // output setup (always append: header already written)
  std::ofstream out;
  gzFile gz = nullptr;

  if (gzip) {
    gz = gzopen(output_file.c_str(), "ab");
    if (!gz) stop("Cannot open gzip output: " + output_file);
  } else {
    out.open(output_file, std::ios::app | std::ios::binary);
    if (!out.is_open()) stop("Cannot open output: " + output_file);
  }

  std::vector<char> buf;
  buf.reserve(1 << 20);  // 1 MB

  auto flush = [&]() {
    if (buf.empty()) return;
    if (gzip) gzwrite(gz, buf.data(), buf.size());
    else out.write(buf.data(), buf.size());
    buf.clear();
  };

  for (int i = 0; i < n_chroms; i++) {

    buf_rstr(buf, chrom[i]); buf_char(buf, '\t');
    buf_int (buf, pos[i]); buf_char(buf, '\t');
    buf_rstr(buf, id[i]); buf_char(buf, '\t');
    buf_rstr(buf, ref[i]); buf_char(buf, '\t');
    buf_rstr(buf, alt[i]); buf_char(buf, '\t');
    buf_rstr(buf, qual[i]); buf_char(buf, '\t');
    buf_rstr(buf, filter_col[i]); buf_char(buf, '\t');
    buf_rstr(buf, info[i]); buf_char(buf, '\t');
    buf_rstr(buf, format_col[i]);

    // samples: flat row-major → offset i*n_samples
    for (int j = 0; j < n_samples; j++) {
      buf_char(buf, '\t');
      buf_rstr(buf, fmt_vec[i * n_samples + j]);
    }
    buf_char(buf, '\n');

    if (buf.size() > (1u << 20)) flush();
  }

  flush();
  if (gzip) gzclose(gz);
  else out.close();
}
