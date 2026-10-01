// src/write_vcf_arrow.cpp

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <zlib.h>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

#include "arrow_c_data.h"

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

// ---- per-sample FORMAT strings from an Arrow string array ----
//
// The array is handed over from R through the Arrow C Data Interface
// (Array$export_to_c()), so the strings are read from Arrow's buffers and
// never become R strings.  Handles utf8 ("u", int32 offsets) and large_utf8
// ("U", int64 offsets).
struct ArrowStrings {
  const uint8_t* valid = nullptr;
  const int32_t* off32 = nullptr;
  const int64_t* off64 = nullptr;
  const char* data = nullptr;
  int64_t offset = 0, length = 0;

  ArrowStrings(const ArrowArray* a, const ArrowSchema* s) {
    const bool large = std::strcmp(s->format, "U") == 0;
    if (!large && std::strcmp(s->format, "u") != 0)
      stop("write_vcf_chunk_cpp: fmt must be a string array (got format '%s')", s->format);
    if (a->n_buffers != 3) stop("write_vcf_chunk_cpp: unexpected string array layout");
    valid = static_cast<const uint8_t*>(a->buffers[0]);
    if (large) off64 = static_cast<const int64_t*>(a->buffers[1]);
    else off32 = static_cast<const int32_t*>(a->buffers[1]);
    data = static_cast<const char*>(a->buffers[2]);
    offset = a->offset;
    length = a->length;
  }

  // writes element k, substituting "." for null/empty
  void write(std::vector<char>& b, int64_t k) const {
    const int64_t i = offset + k;
    if (valid && !(valid[i >> 3] & (1u << (i & 7)))) { b.push_back('.'); return; }
    const int64_t lo = off32 ? off32[i] : off64[i];
    const int64_t hi = off32 ? off32[i + 1] : off64[i + 1];
    if (hi <= lo) { b.push_back('.'); return; }
    b.insert(b.end(), data + lo, data + hi);
  }
};

// Empty ArrowArray / ArrowSchema structs for R to export an Arrow array into
// (Array$export_to_c()), as external pointers that free them when collected.
// [[Rcpp::export]]
List arrow_c_alloc_cpp() {
  ArrowArray* a = new ArrowArray();
  ArrowSchema* s = new ArrowSchema();
  a->release = nullptr;
  s->release = nullptr;
  SEXP pa = PROTECT(R_MakeExternalPtr(a, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(pa, arrow_array_finalizer, TRUE);
  SEXP ps = PROTECT(R_MakeExternalPtr(s, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ps, arrow_schema_finalizer, TRUE);
  List out = List::create(Named("array") = pa, Named("schema") = ps);
  UNPROTECT(2);
  return out;
}

// Appends VCF lines for one chunk.  fmt_array / fmt_schema: external pointers
// (see arrow_c_alloc_cpp()) holding the per-sample FORMAT strings, flat and
// row-major: variant i, sample j -> element i * n_samples + j.  The array is
// released once written.
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
    SEXP fmt_array,
    SEXP fmt_schema,
    int n_samples,
    bool gzip = false
    ) {

  ArrowArray* fa = static_cast<ArrowArray*>(R_ExternalPtrAddr(fmt_array));
  ArrowSchema* fs = static_cast<ArrowSchema*>(R_ExternalPtrAddr(fmt_schema));
  if (!fa || !fs || !fa->release || !fs->release)
    stop("write_vcf_chunk_cpp: fmt is not an exported Arrow array");
  const ArrowStrings fmt(fa, fs);

  int n_chroms = chrom.size();

  if (fmt.length != (int64_t) n_chroms * n_samples)
    stop("fmt length must equal nrow * nsamples");

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
      fmt.write(buf, (int64_t) i * n_samples + j);
    }
    buf_char(buf, '\n');

    if (buf.size() > (1u << 20)) flush();
  }

  flush();
  if (gzip) gzclose(gz);
  else out.close();

  // free the Arrow buffers now rather than when R collects the pointers
  fa->release(fa);
  fs->release(fs);
}
