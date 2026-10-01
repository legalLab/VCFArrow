// src/read_vcf_arrow.cpp
//
// Chunked VCF reader for read_vcf().
//
// Lines are read with zlib (plain or gzip-compressed input) into a buffer
// owned here, and each chunk of variants is parsed directly into Arrow column
// buffers.  The long genotype table (one row per variant x sample) is handed
// to R as an Arrow record batch through the Arrow C Data Interface, so its
// string columns (sample, fmt) never become R strings: creating and later
// converting ~n_variants x n_samples R strings dominated reading time.

#include <Rcpp.h>
#include <zlib.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
using namespace Rcpp;


// Arrow C Data Interface structs and pointer finalizers
#include "arrow_c_data.h"


// ═══════════════════════════════════════════════════════════════════════════════
// II.  Line reader
// ═══════════════════════════════════════════════════════════════════════════════

struct VcfReader {
  gzFile gz;
  std::vector<char> buf;
  std::size_t start, end;
  bool eof;

  explicit VcfReader(const std::string& path)
    : gz(nullptr), buf(1u << 22), start(0), end(0), eof(false) {
    gz = gzopen(path.c_str(), "rb");  // reads uncompressed files transparently
    if (!gz) Rcpp::stop("Cannot open '%s' for reading.", path.c_str());
    gzbuffer(gz, 1u << 20);
  }
  ~VcfReader() { if (gz) gzclose(gz); }

  // Next line without its terminator (LF or CRLF), as [*p, *p + *n).  The
  // line is NUL-terminated and lives in this reader's own mutable buffer
  // (the parser writes NULs into it) until the next call.  False at EOF.
  bool next_line(char** p, std::size_t* n) {
    for (;;) {
      char* s = buf.data() + start;
      char* nl = static_cast<char*>(std::memchr(s, '\n', end - start));
      if (nl || (eof && end > start)) {
        char* e = nl ? nl : buf.data() + end;
        if (!nl && end == buf.size()) {  // final line fills the buffer: room for NUL
          buf.resize(buf.size() + 1);
          s = buf.data() + start;
          e = buf.data() + end;
        }
        start = nl ? static_cast<std::size_t>(nl - buf.data()) + 1 : end;
        if (e > s && e[-1] == '\r') --e;
        *e = '\0';
        *p = s;
        *n = static_cast<std::size_t>(e - s);
        return true;
      }
      if (eof) return false;
      // keep the partial line, then refill
      if (start > 0) {
        std::memmove(buf.data(), buf.data() + start, end - start);
        end -= start;
        start = 0;
      }
      if (end == buf.size()) buf.resize(buf.size() * 2);
      std::size_t room = buf.size() - end;
      int got = gzread(gz, buf.data() + end,
                       static_cast<unsigned>(room > (unsigned) INT_MAX ? INT_MAX : room));
      if (got < 0) {
        int err;
        const char* msg = gzerror(gz, &err);
        Rcpp::stop("Error reading VCF: %s", msg);
      }
      if (got == 0) eof = true;
      end += static_cast<std::size_t>(got);
    }
  }
};

static void reader_finalizer(SEXP x) {
  VcfReader* r = static_cast<VcfReader*>(R_ExternalPtrAddr(x));
  if (r) {
    delete r;
    R_ClearExternalPtr(x);
  }
}

static VcfReader* get_reader(SEXP x) {
  VcfReader* r = static_cast<VcfReader*>(R_ExternalPtrAddr(x));
  if (!r) Rcpp::stop("VCF reader is closed.");
  return r;
}

// [[Rcpp::export]]
SEXP vcf_open_cpp(std::string path) {
  VcfReader* r = new VcfReader(path);
  SEXP x = PROTECT(R_MakeExternalPtr(r, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(x, reader_finalizer, TRUE);
  UNPROTECT(1);
  return x;
}

// [[Rcpp::export]]
void vcf_close_cpp(SEXP handle) { reader_finalizer(handle); }

// Header lines up to and including #CHROM.
// status: 0 = ok, 1 = EOF before #CHROM, 2 = data line before #CHROM
// [[Rcpp::export]]
List vcf_read_header_cpp(SEXP handle) {
  VcfReader* r = get_reader(handle);
  std::vector<std::string> lines;
  int status = 1;
  char* p;
  std::size_t n;
  while (r->next_line(&p, &n)) {
    if (n >= 2 && p[0] == '#' && p[1] == '#') {
      lines.emplace_back(p, n);
    } else if (n >= 6 && std::strncmp(p, "#CHROM", 6) == 0) {
      lines.emplace_back(p, n);
      status = 0;
      break;
    } else {
      status = 2;
      break;
    }
  }
  CharacterVector header(lines.size());
  for (std::size_t i = 0; i < lines.size(); ++i)
    header[i] = Rf_mkCharLen(lines[i].data(), static_cast<int>(lines[i].size()));
  return List::create(Named("header") = header, Named("status") = status);
}


// ═══════════════════════════════════════════════════════════════════════════════
// III.  Arrow column builders
// ═══════════════════════════════════════════════════════════════════════════════

struct BitmapBuilder {
  std::vector<uint8_t> bits;
  void set(int64_t i, bool v) {
    std::size_t byte = static_cast<std::size_t>(i >> 3);
    if (byte >= bits.size()) bits.resize(byte + 1, 0);
    if (v) bits[byte] |= static_cast<uint8_t>(1u << (i & 7));
  }
};

template <typename T>
struct PrimitiveBuilder {
  std::vector<T> values;
  BitmapBuilder valid;
  int64_t null_count = 0;
  void append(T v, bool is_valid) {
    int64_t i = static_cast<int64_t>(values.size());
    values.push_back(v);
    valid.set(i, is_valid);
    if (!is_valid) ++null_count;
  }
};

struct StringBuilder {
  std::vector<int32_t> offsets{0};
  std::vector<char> data;
  void append(const char* s, std::size_t n) {
    if (data.size() + n > static_cast<std::size_t>(INT32_MAX))
      Rcpp::stop("Genotype strings of one chunk exceed 2 GB; use a smaller chunk_size.");
    data.insert(data.end(), s, s + n);
    offsets.push_back(static_cast<int32_t>(data.size()));
  }
};

// Owns every buffer of one record batch until Arrow releases it.
struct BatchData {
  int64_t length = 0;
  PrimitiveBuilder<int32_t> row_id, a1, a2;
  BitmapBuilder phased;
  StringBuilder sample, fmt;
  PrimitiveBuilder<double> DP, GQ, ADR;

  static const int n_cols = 9;
  ArrowArray children[n_cols];
  ArrowArray* child_ptrs[n_cols];
  const void* child_buffers[n_cols][3];
  const void* parent_buffers[1];
};

static void release_child(ArrowArray* a) { a->release = nullptr; }  // freed by parent

static void release_batch(ArrowArray* a) {
  BatchData* d = static_cast<BatchData*>(a->private_data);
  for (int k = 0; k < BatchData::n_cols; ++k)
    if (d->children[k].release) d->children[k].release(&d->children[k]);
  delete d;
  a->release = nullptr;
}

static void set_child(BatchData* d, int k, int64_t null_count,
                      const void* b0, const void* b1, const void* b2, int n_buffers) {
  d->child_buffers[k][0] = b0;
  d->child_buffers[k][1] = b1;
  d->child_buffers[k][2] = b2;
  ArrowArray& c = d->children[k];
  c.length = d->length;
  c.null_count = null_count;
  c.offset = 0;
  c.n_buffers = n_buffers;
  c.n_children = 0;
  c.buffers = d->child_buffers[k];
  c.children = nullptr;
  c.dictionary = nullptr;
  c.release = release_child;
  c.private_data = nullptr;
  d->child_ptrs[k] = &c;
}

template <typename T>
static void set_primitive(BatchData* d, int k, PrimitiveBuilder<T>& b) {
  b.valid.bits.resize(static_cast<std::size_t>((d->length + 7) / 8), 0);
  set_child(d, k, b.null_count, b.null_count ? b.valid.bits.data() : nullptr,
            b.values.data(), nullptr, 2);
}

// Column schema (static strings; nothing to free but the structs)
static const char* col_names[BatchData::n_cols] =
  {".row_id", "sample", "a1", "a2", "phased", "fmt", "DP", "GQ", "ADR"};
static const char* col_formats[BatchData::n_cols] =
  {"i", "u", "i", "i", "b", "u", "g", "g", "g"};

struct SchemaData {
  ArrowSchema children[BatchData::n_cols];
  ArrowSchema* child_ptrs[BatchData::n_cols];
};

static void release_child_schema(ArrowSchema* s) { s->release = nullptr; }

static void release_schema(ArrowSchema* s) {
  SchemaData* d = static_cast<SchemaData*>(s->private_data);
  for (int k = 0; k < BatchData::n_cols; ++k)
    if (d->children[k].release) d->children[k].release(&d->children[k]);
  delete d;
  s->release = nullptr;
}

static SEXP export_batch(BatchData* d) {
  set_primitive(d, 0, d->row_id);
  set_child(d, 1, 0, nullptr, d->sample.offsets.data(), d->sample.data.data(), 3);
  set_primitive(d, 2, d->a1);
  set_primitive(d, 3, d->a2);
  d->phased.bits.resize(static_cast<std::size_t>((d->length + 7) / 8), 0);
  set_child(d, 4, 0, nullptr, d->phased.bits.data(), nullptr, 2);
  set_child(d, 5, 0, nullptr, d->fmt.offsets.data(), d->fmt.data.data(), 3);
  set_primitive(d, 6, d->DP);
  set_primitive(d, 7, d->GQ);
  set_primitive(d, 8, d->ADR);

  d->parent_buffers[0] = nullptr;
  ArrowArray* a = new ArrowArray;
  a->length = d->length;
  a->null_count = 0;
  a->offset = 0;
  a->n_buffers = 1;
  a->n_children = BatchData::n_cols;
  a->buffers = d->parent_buffers;
  a->children = d->child_ptrs;
  a->dictionary = nullptr;
  a->release = release_batch;
  a->private_data = d;

  SchemaData* sd = new SchemaData;
  for (int k = 0; k < BatchData::n_cols; ++k) {
    ArrowSchema& c = sd->children[k];
    c.format = col_formats[k];
    c.name = col_names[k];
    c.metadata = nullptr;
    c.flags = ARROW_FLAG_NULLABLE;
    c.n_children = 0;
    c.children = nullptr;
    c.dictionary = nullptr;
    c.release = release_child_schema;
    c.private_data = nullptr;
    sd->child_ptrs[k] = &c;
  }
  ArrowSchema* s = new ArrowSchema;
  s->format = "+s";
  s->name = "";
  s->metadata = nullptr;
  s->flags = 0;
  s->n_children = BatchData::n_cols;
  s->children = sd->child_ptrs;
  s->dictionary = nullptr;
  s->release = release_schema;
  s->private_data = sd;

  SEXP pa = PROTECT(R_MakeExternalPtr(a, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(pa, arrow_array_finalizer, TRUE);
  SEXP ps = PROTECT(R_MakeExternalPtr(s, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ps, arrow_schema_finalizer, TRUE);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, pa);
  SET_VECTOR_ELT(out, 1, ps);
  UNPROTECT(3);
  return out;
}


// ═══════════════════════════════════════════════════════════════════════════════
// IV.  Parsing
// ═══════════════════════════════════════════════════════════════════════════════

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

// Parse up to chunk_size variant lines, obtained from next_line(&p, &n)
// (NUL-terminated, in a mutable buffer).  Returns n = 0 when there are none.
// Per-variant fields come back as R vectors; the long genotype table as an
// Arrow record batch (list(array, schema) of external pointers) in
// variant-major row order, with .row_id = row_offset + 1, 2, ...
template <typename NextLine>
static List parse_chunk(NextLine next_line, int chunk_size,
                        const CharacterVector& samples, int row_offset) {
  const int nsamples = samples.size();

  std::vector<std::string> sample_names(nsamples);
  for (int j = 0; j < nsamples; ++j) sample_names[j] = as<std::string>(samples[j]);

  std::vector<std::vector<std::string>> fixed(9);
  BatchData* d = new BatchData;
  std::size_t n_cell = static_cast<std::size_t>(chunk_size) * nsamples;
  d->row_id.values.reserve(n_cell);
  d->a1.values.reserve(n_cell);
  d->a2.values.reserve(n_cell);
  d->DP.values.reserve(n_cell);
  d->GQ.values.reserve(n_cell);
  d->ADR.values.reserve(n_cell);
  d->sample.offsets.reserve(n_cell + 1);
  d->fmt.offsets.reserve(n_cell + 1);

  int n = 0;
  char* line;
  std::size_t line_len;

  try {
    while (n < chunk_size && next_line(&line, &line_len)) {

      char* ptr = line;
      const char* field_start;
      const char* field_end;

      // ---- parse first 9 fields (NUL-terminated in place) ----
      const char* fields[9];

      for (int k = 0; k < 9; k++) {
        field_start = ptr;
        while (*ptr && *ptr != '\t') ptr++;
        fields[k] = field_start;

        if (*ptr == '\t') {
          *ptr = '\0';  // our own buffer
          ptr++;
        }
      }

      for (int k = 0; k < 9; k++) fixed[k].emplace_back(fields[k]);

      // ---- parse FORMAT keys (pointer-based) ----
      int dp_pos = -1, gq_pos = -1, ad_pos = -1;

      int key_idx = 0;
      const char* p = fields[8];
      const char* key_start = p;

      while (true) {
        if (*p == ':' || *p == '\0') {
          if (token_eq(key_start, p, "DP")) dp_pos = key_idx;
          else if (token_eq(key_start, p, "GQ")) gq_pos = key_idx;
          else if (token_eq(key_start, p, "AD")) ad_pos = key_idx;

          key_idx++;
          if (*p == '\0') break;
          key_start = p + 1;
        }
        p++;
      }

      const int32_t row_id = row_offset + n + 1;

      // ---- samples ----
      for (int j = 0; j < nsamples; j++) {

        field_start = ptr;
        while (*ptr && *ptr != '\t') ptr++;
        field_end = ptr;

        const int64_t cell = d->length++;
        d->row_id.append(row_id, true);
        d->sample.append(sample_names[j].data(), sample_names[j].size());
        d->fmt.append(field_start, static_cast<std::size_t>(field_end - field_start));

        // ---- parse FORMAT values inline ----
        int val_idx = 0;
        const char* v = field_start;
        const char* val_start = v;

        double dp_val = NA_REAL;
        double gq_val = NA_REAL;
        double ad_ratio = NA_REAL;

        // GT parsing
        int a1 = NA_INTEGER, a2 = NA_INTEGER;
        bool phased = false;
        if (field_end - field_start >= 3) {
          char c1 = field_start[0];
          char sep = field_start[1];
          char c2 = field_start[2];

          phased = (sep == '|');

          if (c1 != '.' && c2 != '.') {
            a1 = c1 - '0';
            a2 = c2 - '0';
          }
        }
        d->a1.append(a1, a1 != NA_INTEGER);
        d->a2.append(a2, a2 != NA_INTEGER);
        d->phased.set(cell, phased);

        while (true) {

          if (v == field_end || *v == ':') {

            // --- DP ---
            if (val_idx == dp_pos) {
              dp_val = fast_atof(val_start, v);
            }

            // --- GQ ---
            if (val_idx == gq_pos) {
              gq_val = fast_atof(val_start, v);
            }

            // --- AD (compute ratio inline) ---
            if (val_idx == ad_pos) {

              const char* a = val_start;
              const char* b = val_start;

              // find comma
              while (b < v && *b != ',') b++;

              if (b < v) {
                double ref = fast_atof(a, b);
                double alt = fast_atof(b + 1, v);

                if (!R_IsNA(ref) && !R_IsNA(alt) && (ref + alt) > 0.0) {
                  ad_ratio = alt / (ref + alt);
                }
              }
            }

            if (v == field_end) break;

            val_idx++;
            val_start = v + 1;
          }
          v++;
        }

        d->DP.append(dp_val, !R_IsNA(dp_val));
        d->GQ.append(gq_val, !R_IsNA(gq_val));
        d->ADR.append(ad_ratio, !R_IsNA(ad_ratio));

        if (*ptr == '\t') ptr++;
      }
      n++;
    }
  } catch (...) {
    delete d;
    throw;
  }

  if (n == 0) {
    delete d;
    return List::create(Named("n") = 0);
  }

  CharacterMatrix variants(n, 7);
  CharacterVector info(n), format(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 7; ++k) variants(i, k) = fixed[k][i];
    info[i] = fixed[7][i];
    format[i] = fixed[8][i];
  }

  SEXP batch = PROTECT(export_batch(d));
  List out = List::create(
    Named("n") = n,
    Named("variants") = variants,
    Named("info") = info,
    Named("format") = format,
    Named("array") = VECTOR_ELT(batch, 0),
    Named("schema") = VECTOR_ELT(batch, 1)
  );
  UNPROTECT(1);
  return out;
}

// Serial reading: read and parse the next chunk_size lines.
// [[Rcpp::export]]
List vcf_read_chunk_cpp(SEXP handle, int chunk_size,
                        CharacterVector samples, int row_offset) {
  VcfReader* r = get_reader(handle);
  return parse_chunk([r](char** p, std::size_t* n) { return r->next_line(p, n); },
                     chunk_size, samples, row_offset);
}

// Parallel reading, step 1 (main process): the next chunk_size lines as raw
// bytes, each line terminated by '\n'.  Returns list(n, raw); n = 0 at EOF.
// [[Rcpp::export]]
List vcf_read_raw_cpp(SEXP handle, int chunk_size) {
  VcfReader* r = get_reader(handle);
  std::vector<char> bytes;
  int n = 0;
  char* p;
  std::size_t len;
  while (n < chunk_size && r->next_line(&p, &len)) {
    bytes.insert(bytes.end(), p, p + len);
    bytes.push_back('\n');
    ++n;
  }
  RawVector raw(bytes.size());
  if (!bytes.empty()) std::memcpy(RAW(raw), bytes.data(), bytes.size());
  return List::create(Named("n") = n, Named("raw") = raw);
}

// Parallel reading, step 2 (worker): parse the lines of vcf_read_raw_cpp().
// [[Rcpp::export]]
List vcf_parse_raw_cpp(RawVector raw, CharacterVector samples, int row_offset) {
  std::vector<char> buf(raw.begin(), raw.end());  // parsed in place
  std::size_t pos = 0;
  auto next_line = [&buf, &pos](char** p, std::size_t* n) {
    if (pos >= buf.size()) return false;
    char* s = buf.data() + pos;
    char* nl = static_cast<char*>(std::memchr(s, '\n', buf.size() - pos));
    if (!nl) Rcpp::stop("vcf_parse_raw_cpp: unterminated line");
    *nl = '\0';
    *p = s;
    *n = static_cast<std::size_t>(nl - s);
    pos = static_cast<std::size_t>(nl - buf.data()) + 1;
    return true;
  };
  // number of lines, which sizes the column buffers
  const int n_lines = static_cast<int>(std::count(buf.begin(), buf.end(), '\n'));
  return parse_chunk(next_line, n_lines, samples, row_offset);
}
