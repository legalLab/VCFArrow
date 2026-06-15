// vcf_wide.cpp
//
// Compile (standalone): Rcpp::sourceCpp("vcf_wide.cpp")
// In a package: place in src/, run devtools::document()
//
// Matrix layout convention (all individual-level writers):
//   rows = samples, cols = variants (n_samples x n_var, R column-major)
//   element (s, v)  →  raw offset  s + v * n_samples
//
// Population-level matrices (BayesScan, Treemix, Migrate):
//   rows = pops, cols = variants (n_pops x n_var)
//
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
using namespace Rcpp;


// ═══════════════════════════════════════════════════════════════════════════════
// I.  Buffered file writer
// ═══════════════════════════════════════════════════════════════════════════════

static const std::size_t IO_BUF = 1u << 22;   // 4 MiB

struct WFile {
  FILE* fp;
  std::vector<char> _buf;
  WFile(const std::string& path, bool append) : fp(nullptr), _buf(IO_BUF) {
    fp = std::fopen(path.c_str(), append ? "ab" : "wb");
    if (!fp) Rcpp::stop("Cannot open '%s' for writing.", path.c_str());
    std::setvbuf(fp, _buf.data(), _IOFBF, IO_BUF);
  }
  ~WFile() { if (fp) std::fclose(fp); }
  WFile(const WFile&) = delete;
  WFile& operator=(const WFile&) = delete;
};

static inline void wf(WFile& w, char c) { std::fputc(c, w.fp); }
static inline void wf(WFile& w, const char* s) { std::fputs(s, w.fp); }
static inline void wf_int(WFile& w, int v) {
  char tmp[12]; std::snprintf(tmp, sizeof(tmp), "%d", v); std::fputs(tmp, w.fp);
}
static inline void wf_allele(WFile& w, int a, const char* rp,
                             const char* ap, const char* miss) {
  if (a == NA_INTEGER) { std::fputs(miss, w.fp); return; }
  std::fputs((a == 0) ? rp : ap, w.fp);
}


// ═══════════════════════════════════════════════════════════════════════════════
// II.  Shared nucleotide helpers
// ═══════════════════════════════════════════════════════════════════════════════

// Nucleotide → 2-digit Genepop / Related code  (A=01 C=02 G=03 T=04)
static inline const char* nuc_to_digit(char c) {
  switch (c) {
    case 'A': return "01";
    case 'C': return "02";
    case 'G': return "03";
    case 'T': return "04";
    default: return "00";  // unknown / indel remnant
  }
}

// IUPAC ambiguity code for a heterozygous position (ref nucleotide, alt nucleotide)
static inline char iupac_het(char r, char a) {
  if (r > a) { char t = r; r = a; a = t; }  // canonical sort
  if (r == 'A') {
    if (a == 'C') return 'M';
    if (a == 'G') return 'R';
    if (a == 'T') return 'W';
  } else if (r == 'C') {
    if (a == 'G') return 'S';
    if (a == 'T') return 'Y';
  } else if (r == 'G' && a == 'T') {
    return 'K';
  }
  return 'N';
}


// ═══════════════════════════════════════════════════════════════════════════════
// III.  SmartSNP
// ═══════════════════════════════════════════════════════════════════════════════

//' Write the sample-names header for a SmartSNP file (creates / truncates)
// [[Rcpp::export]]
void write_smartsnp_header_cpp(const CharacterVector& samples,
                               const std::string& out_file) {
  WFile w(out_file, false);
  const int n = samples.size();
  for (int i = 0; i < n; ++i) { if (i) wf(w, ' '); wf(w, CHAR(STRING_ELT(samples, i))); }
  wf(w, '\n');
}

//' Encode (0/1/2/9) and append one chunk of variants to a SmartSNP file
//' a1_mat / a2_mat: integer matrices (n_samples x n_chunk_var)
// [[Rcpp::export]]
void write_smartsnp_chunk_cpp(const IntegerMatrix& a1_mat,
                              const IntegerMatrix& a2_mat,
                              const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, true);
  for (int v = 0; v < nv; ++v) {
    for (int s = 0; s < ns; ++s) {
      if (s) wf(w, ' ');
      const int v1 = a1_mat(s, v), v2 = a2_mat(s, v);
      if      (v1 == NA_INTEGER || v2 == NA_INTEGER) wf(w, '9');
      else if (v1 == 0 && v2 == 0) wf(w, '0');
      else if (v1 == 1 && v2 == 1) wf(w, '2');
      else wf(w, '1');
    }
    wf(w, '\n');
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// IV.  STRUCTURE
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a STRUCTURE input file from full-dataset integer matrices
//' method_int: 0 = Simple (S), 1 = FastStructure (F)
// [[Rcpp::export]]
void write_structure_cpp(const IntegerMatrix& a1_mat,
                         const IntegerMatrix& a2_mat,
                         const CharacterVector& samples,
                         const IntegerVector& group_ids,
                         int method_int,
                         const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, false);
  for (int s = 0; s < ns; ++s) {
    const char* sn = CHAR(STRING_ELT(samples, s));
    for (int hap = 0; hap < 2; ++hap) {
      wf(w, sn); wf(w, '\t'); wf_int(w, group_ids[s]);
      if (method_int == 1) wf(w, "\t0\t0\t0\t0");
      for (int v = 0; v < nv; ++v) {
        wf(w, '\t');
        const int a = (hap == 0) ? a1_mat(s, v) : a2_mat(s, v);
        if (a == NA_INTEGER) {
          wf(w, "-9"); continue;
        }
        wf(w, (a == 0) ? '0' : '1');
      }
      wf(w, '\n');
    }
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// V.  Arlequin
// ═══════════════════════════════════════════════════════════════════════════════

//' Write an Arlequin (.arp) input file
// [[Rcpp::export]]
void write_arlequin_cpp(const IntegerMatrix& a1_mat,
                        const IntegerMatrix& a2_mat,
                        const CharacterVector& REF,
                        const CharacterVector& ALT,
                        const CharacterVector& samples,
                        const CharacterVector& group_names,
                        const IntegerVector& group_sizes,
                        const std::string& out_file) {
  const int nv = a1_mat.ncol(), ng = group_names.size();
  WFile w(out_file, false);
  std::vector<const char*> rp(nv), ap(nv);
  for (int v = 0; v < nv; ++v) { rp[v] = CHAR(STRING_ELT(REF,v)); ap[v] = CHAR(STRING_ELT(ALT,v)); }

  wf(w, "[Profile]\n\nTitle = 'Generated by VCFArrow'\nNbSamples = "); wf_int(w, ng); wf(w, '\n');
  wf(w, "GenotypicData = 1\nLocusSeparator = WHITESPACE\nGameticPhase = 0\n");
  wf(w, "MissingData = '?'\nDataType = STANDARD\n\n[Data]\n[[Samples]]\n\n");
  wf(w, "#There are "); wf_int(w, nv); wf(w, " SNPs\n\n");

  int off = 0;
  for (int g = 0; g < ng; ++g) {
    const int gsz = group_sizes[g];
    wf(w, "SampleName = '"); wf(w, CHAR(STRING_ELT(group_names,g))); wf(w, "'\nSampleSize = ");
    wf_int(w, gsz); wf(w, "\nSampleData={\n");
    for (int i = 0; i < gsz; ++i) {
      const int s = off + i;
      wf(w, CHAR(STRING_ELT(samples,s))); wf(w, "\t1\t");
      for (int v = 0; v < nv; ++v) { if (v) wf(w,' '); wf_allele(w, a1_mat(s,v), rp[v], ap[v], "?"); }
      wf(w, '\n'); wf(w, "\t\t\t\t\t\t");
      for (int v = 0; v < nv; ++v) { if (v) wf(w,' '); wf_allele(w, a2_mat(s,v), rp[v], ap[v], "?"); }
      wf(w, '\n');
    }
    wf(w, "}\n\n"); off += gsz;
  }
  wf(w, "[[Structure]]\nStructureName = 'One Group'\nNbGroups = 1\n\nGroup = {\n");
  for (int g = 0; g < ng; ++g) { wf(w,"\t\t\""); wf(w, CHAR(STRING_ELT(group_names,g))); wf(w,"\"\n"); }
  wf(w, "}\n");
}


// ═══════════════════════════════════════════════════════════════════════════════
// VI.  Genepop
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a Genepop input file
// [[Rcpp::export]]
void write_genepop_cpp(const IntegerMatrix& a1_mat,
                       const IntegerMatrix& a2_mat,
                       const CharacterVector& REF,
                       const CharacterVector& ALT,
                       const CharacterVector& samples,
                       const CharacterVector& group_names,
                       const IntegerVector& group_sizes,
                       const CharacterVector& loci,
                       const std::string& out_file) {
  const int nv = a1_mat.ncol(), ng = group_names.size();
  WFile w(out_file, false);

  wf(w, "Title = 'Generated by VCFArrow'\n");
  for (int v = 0; v < nv; ++v) { if (v) wf(w,", "); wf(w, CHAR(STRING_ELT(loci,v))); }
  wf(w, '\n');

  int off = 0;
  for (int g = 0; g < ng; ++g) {
    const int gsz = group_sizes[g];
    wf(w, "pop "); wf(w, CHAR(STRING_ELT(group_names,g))); wf(w, '\n');
    for (int i = 0; i < gsz; ++i) {
      const int s = off + i;
      wf(w, CHAR(STRING_ELT(samples,s))); wf(w, ", ");
      for (int v = 0; v < nv; ++v) {
        if (v) wf(w, ' ');
        const char* r = CHAR(STRING_ELT(REF,v)), *a = CHAR(STRING_ELT(ALT,v));
        const int  v1 = a1_mat(s,v), v2 = a2_mat(s,v);
        if (v1 == NA_INTEGER || v2 == NA_INTEGER) { wf(w, "0000"); continue; }
        const char* c1 = nuc_to_digit((v1==0) ? r[0] : a[0]);
        const char* c2 = nuc_to_digit((v2==0) ? r[0] : a[0]);
        if (std::strcmp(c1,c2) > 0) { const char* t=c1; c1=c2; c2=t; }
        wf(w, c1); wf(w, c2);
      }
      wf(w, '\n');
    }
    off += gsz;
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// VII.  fineRADstructure
// ═══════════════════════════════════════════════════════════════════════════════

//' Write the sample-names header for a fineRADstructure file (creates / truncates)
// [[Rcpp::export]]
void write_fineradstructure_header_cpp(const CharacterVector& samples,
                                       const std::string& out_file) {
  WFile w(out_file, false);
  const int n = samples.size();
  for (int i = 0; i < n; ++i) { if (i) wf(w,'\t'); wf(w, CHAR(STRING_ELT(samples,i))); }
  wf(w, '\n');
}

//' Encode and append one chunk to a fineRADstructure file
//' REF and ALT here are per-chunk (length n_chunk_var).
// [[Rcpp::export]]
void write_fineradstructure_chunk_cpp(const IntegerMatrix& a1_mat,
                                      const IntegerMatrix& a2_mat,
                                      const CharacterVector& REF,
                                      const CharacterVector& ALT,
                                      const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, true);
  for (int v = 0; v < nv; ++v) {
    const char* rp = CHAR(STRING_ELT(REF,v)), *ap = CHAR(STRING_ELT(ALT,v));
    for (int s = 0; s < ns; ++s) {
      if (s) wf(w, '\t');
      const int v1 = a1_mat(s,v), v2 = a2_mat(s,v);
      if (v1 == NA_INTEGER || v2 == NA_INTEGER) { wf(w, "NA,NA"); continue; }
      wf(w, (v1==0) ? rp : ap); wf(w, ','); wf(w, (v2==0) ? rp : ap);
    }
    wf(w, '\n');
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// VIII.  Treemix
// ═══════════════════════════════════════════════════════════════════════════════

//' Write the population-names header for a Treemix file (creates / truncates)
// [[Rcpp::export]]
void write_treemix_header_cpp(const CharacterVector& pop_names,
                              const std::string& out_file) {
  WFile w(out_file, false);
  const int n = pop_names.size();
  for (int p = 0; p < n; ++p) { if (p) wf(w,' '); wf(w, CHAR(STRING_ELT(pop_names,p))); }
  wf(w, '\n');
}

//' Append one chunk of population allele counts to a Treemix file
//' ref_mat / alt_mat: integer matrices (n_pops x n_chunk_var)
// [[Rcpp::export]]
void write_treemix_chunk_cpp(const IntegerMatrix& ref_mat,
                             const IntegerMatrix& alt_mat,
                             const std::string& out_file) {
  const int np = ref_mat.nrow(), nv = ref_mat.ncol();
  WFile w(out_file, true);
  for (int v = 0; v < nv; ++v) {
    for (int p = 0; p < np; ++p) {
      if (p) wf(w,' ');
      wf_int(w, ref_mat(p,v)); wf(w,','); wf_int(w, alt_mat(p,v));
    }
    wf(w, '\n');
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// IX.  BayesScan
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a BayesScan input file
//' ref_mat / alt_mat / n_obs_mat: integer matrices (n_pops x n_var)
// [[Rcpp::export]]
void write_bayescan_cpp(const IntegerMatrix& ref_mat,
                        const IntegerMatrix& alt_mat,
                        const IntegerMatrix& n_obs_mat,
                        const CharacterVector& pop_names,
                        const std::string& out_file) {
  const int np = ref_mat.nrow(), nv = ref_mat.ncol();
  WFile w(out_file, false);
  wf(w,"[loci]="); wf_int(w,nv); wf(w,'\n'); wf(w,'\n');
  wf(w,"[populations]="); wf_int(w,np); wf(w,'\n');
  for (int p = 0; p < np; ++p) {
    wf(w,'\n'); wf(w,"[pop]="); wf_int(w,p+1); wf(w,'\n');
    for (int v = 0; v < nv; ++v) {
      wf_int(w,n_obs_mat(p,v)); wf(w,"\t2\t");
      wf_int(w,ref_mat(p,v)); wf(w,'\t');
      wf_int(w,alt_mat(p,v)); wf(w,'\n');
    }
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// X.  BayesAss
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a BayesAss3 (.immanc) input file
// [[Rcpp::export]]
void write_bayesass_cpp(const IntegerMatrix& a1_mat,
                        const IntegerMatrix& a2_mat,
                        const CharacterVector& REF,
                        const CharacterVector& ALT,
                        const CharacterVector& samples,
                        const CharacterVector& group_names,
                        const IntegerVector& group_sizes,
                        const CharacterVector& loci,
                        const std::string& out_file) {
  const int nv = a1_mat.ncol(), np = group_names.size();
  WFile w(out_file, false);
  std::vector<const char*> rp(nv), ap(nv);
  for (int v = 0; v < nv; ++v) { rp[v]=CHAR(STRING_ELT(REF,v)); ap[v]=CHAR(STRING_ELT(ALT,v)); }

  wf(w,"N\t"); wf_int(w,np); wf(w,'\t'); wf_int(w,nv); wf(w,'\n');
  for (int v = 0; v < nv; ++v) { if (v) wf(w,'\t'); wf(w,'1'); }
  wf(w,'\n');

  int off = 0;
  for (int g = 0; g < np; ++g) {
    const int gsz = group_sizes[g];
    wf_int(w, 2*gsz); wf(w,'\t'); wf(w, CHAR(STRING_ELT(group_names,g))); wf(w,'\n');
    for (int v = 0; v < nv; ++v) {
      wf(w, CHAR(STRING_ELT(loci,v))); wf(w,'\n');
      for (int i = 0; i < gsz; ++i) { wf_allele(w,a1_mat(off+i,v),rp[v],ap[v],"?"); wf(w,'\n'); }
      for (int i = 0; i < gsz; ++i) { wf_allele(w,a2_mat(off+i,v),rp[v],ap[v],"?"); wf(w,'\n'); }
    }
    off += gsz;
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// XI.  Migrate-N
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a Migrate-N allele-count input file
// [[Rcpp::export]]
void write_migrate_cpp(const IntegerMatrix& ref_mat,
                       const IntegerMatrix& alt_mat,
                       const IntegerMatrix& n_obs_mat,
                       const CharacterVector& pop_names,
                       const CharacterVector& loci,
                       const std::string& out_file) {
  const int np = ref_mat.nrow(), nv = ref_mat.ncol();
  WFile w(out_file, false);
  wf_int(w,np); wf(w,' '); wf_int(w,nv); wf(w,'\n');
  for (int p = 0; p < np; ++p) {
    int mx = 0;
    for (int v = 0; v < nv; ++v) if (n_obs_mat(p,v) > mx) mx = n_obs_mat(p,v);
    wf_int(w,mx); wf(w,' '); wf(w, CHAR(STRING_ELT(pop_names,p))); wf(w,'\n');
    for (int v = 0; v < nv; ++v) {
      wf_int(w,ref_mat(p,v)); wf(w,' '); wf_int(w,alt_mat(p,v)); wf(w,'\n');
    }
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// XII.  Related
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a Related input file
//'
//' One tab-delimited line per sample.  Each locus contributes two adjacent
//' tab-separated columns (allele1 code, allele2 code), interleaved across all
//' loci:  a1_L1 \\t a2_L1 \\t a1_L2 \\t a2_L2 ...
//' Nucleotide codes: A=01 C=02 G=03 T=04, canonical order (lower first).
//' Missing allele: NA.
// [[Rcpp::export]]
void write_related_cpp(const IntegerMatrix& a1_mat,
                       const IntegerMatrix& a2_mat,
                       const CharacterVector& REF,
                       const CharacterVector& ALT,
                       const CharacterVector& samples,
                       const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, false);

  for (int s = 0; s < ns; ++s) {
    if (s) wf(w, '\n');
    bool first = true;
    for (int v = 0; v < nv; ++v) {
      const int v1 = a1_mat(s,v), v2 = a2_mat(s,v);
      const char* rp = CHAR(STRING_ELT(REF,v));
      const char* ap = CHAR(STRING_ELT(ALT,v));

      // Resolve allele characters (NA stays as NA marker)
      const char* c1;
      const char* c2;
      if (v1 == NA_INTEGER) {
        c1 = "NA";
      } else {
        c1 = nuc_to_digit((v1==0) ? rp[0] : ap[0]);
        // map "00" (unknown nuc) to "NA"
        if (c1[0] == '0' && c1[1] == '0') c1 = "NA";
      }
      if (v2 == NA_INTEGER) {
        c2 = "NA";
      } else {
        c2 = nuc_to_digit((v2==0) ? rp[0] : ap[0]);
        if (c2[0] == '0' && c2[1] == '0') c2 = "NA";
      }

      // Canonical ordering: lower numeric code first (only when both known)
      if (c1[0] != 'N' && c2[0] != 'N' && std::strcmp(c1,c2) > 0) {
        const char* t = c1; c1 = c2; c2 = t;
      }

      if (!first) wf(w, '\t');
      wf(w, c1); wf(w, '\t'); wf(w, c2);
      first = false;
    }
  }
  wf(w, '\n');
}


// ═══════════════════════════════════════════════════════════════════════════════
// XIII.  Apparent
// ═══════════════════════════════════════════════════════════════════════════════

//' Write an Apparent input file
//'
//' One tab-delimited line per sample.
//' Format: key \\t locus1 \\t locus2 ...
//' Each locus cell: "allele1/allele2" nucleotide strings, canonical
//' (alphabetically lower allele first).  Missing: "?/?".
// [[Rcpp::export]]
void write_apparent_cpp(const IntegerMatrix& a1_mat,
                        const IntegerMatrix& a2_mat,
                        const CharacterVector& REF,
                        const CharacterVector& ALT,
                        const CharacterVector& samples,
                        const std::string& key,
                        const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, false);

  for (int s = 0; s < ns; ++s) {
    wf(w, key.c_str());
    for (int v = 0; v < nv; ++v) {
      wf(w, '\t');
      const int v1 = a1_mat(s,v), v2 = a2_mat(s,v);
      const char* rp = CHAR(STRING_ELT(REF,v));
      const char* ap = CHAR(STRING_ELT(ALT,v));
      if (v1 == NA_INTEGER || v2 == NA_INTEGER) { wf(w,"?/?"); continue; }
      const char* a1s = (v1==0) ? rp : ap;
      const char* a2s = (v2==0) ? rp : ap;
      // Canonical alphabetical order
      if (std::strcmp(a1s, a2s) > 0) { const char* t=a1s; a1s=a2s; a2s=t; }
      wf(w, a1s); wf(w, '/'); wf(w, a2s);
    }
    wf(w, '\n');
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// XIV.  FASTA
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a FASTA file (one line per sequence, IUPAC nucleotide codes)
//'
//' hom-ref → REF nucleotide
//' hom-alt → ALT nucleotide
//' het     → IUPAC ambiguity code
//' missing → ?
// [[Rcpp::export]]
void write_fasta_cpp(const IntegerMatrix& a1_mat,
                     const IntegerMatrix& a2_mat,
                     const CharacterVector& REF,
                     const CharacterVector& ALT,
                     const CharacterVector& samples,
                     const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, false);

  for (int s = 0; s < ns; ++s) {
    wf(w, '>'); wf(w, CHAR(STRING_ELT(samples,s))); wf(w, '\n');
    for (int v = 0; v < nv; ++v) {
      const int v1 = a1_mat(s,v), v2 = a2_mat(s,v);
      if (v1 == NA_INTEGER || v2 == NA_INTEGER) { wf(w,'?'); continue; }
      const char* rp = CHAR(STRING_ELT(REF,v));
      const char* ap = CHAR(STRING_ELT(ALT,v));
      const char r = rp[0], a = ap[0];
      if (v1 == 0 && v2 == 0) wf(w, r);
      else if (v1 == 1 && v2 == 1) wf(w, a);
      else wf(w, iupac_het(r, a));
    }
    wf(w, '\n');
  }
}


// ═══════════════════════════════════════════════════════════════════════════════
// XV.  Nexus / SNAPP  (shared writer, format_int controls encoding)
// ═══════════════════════════════════════════════════════════════════════════════

//' Write a NEXUS file (Nexus DNA or SNAPP encoding)
//'
//' format_int:
//'   0 = Nexus DNA — IUPAC nucleotide codes (same sequence content as FASTA)
//'   1 = SNAPP     — 0 (hom-ref) / 1 (het) / 2 (hom-alt) / ? (missing)
//'
//' Sample names are left-padded to align sequences in the MATRIX block.
// [[Rcpp::export]]
void write_nexus_cpp(const IntegerMatrix& a1_mat,
                     const IntegerMatrix& a2_mat,
                     const CharacterVector& REF,
                     const CharacterVector& ALT,
                     const CharacterVector& samples,
                     int format_int,
                     const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, false);

  // Compute max sample-name length for alignment padding
  int max_len = 0;
  for (int s = 0; s < ns; ++s) {
    int len = static_cast<int>(std::strlen(CHAR(STRING_ELT(samples,s))));
    if (len > max_len) max_len = len;
  }
  const int pad_width = max_len + 3;   // 3 spaces minimum gap

  // Header
  wf(w, "#NEXUS\n");
  wf(w, "[Data written by VCFArrow]\n");
  wf(w, "BEGIN DATA;\n");
  wf(w, "  DIMENSIONS NTAX="); wf_int(w,ns);
  wf(w, " NCHAR="); wf_int(w,nv); wf(w,";\n");
  wf(w, "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;\n");
  wf(w, "  MATRIX\n");

  for (int s = 0; s < ns; ++s) {
    wf(w, "    ");
    const char* sname = CHAR(STRING_ELT(samples,s));
    wf(w, sname);
    // Right-pad name to pad_width
    int name_len = static_cast<int>(std::strlen(sname));
    for (int k = name_len; k < pad_width; ++k) wf(w, ' ');

    for (int v = 0; v < nv; ++v) {
      const int v1 = a1_mat(s,v), v2 = a2_mat(s,v);
      if (v1 == NA_INTEGER || v2 == NA_INTEGER) { wf(w,'?'); continue; }

      if (format_int == 1) {
        // SNAPP: 0/1/2
        if (v1 == 0 && v2 == 0) wf(w,'0');
        else if (v1 == 1 && v2 == 1) wf(w,'2');
        else wf(w,'1');
      } else {
        // Nexus DNA: IUPAC
        const char r = CHAR(STRING_ELT(REF,v))[0];
        const char a = CHAR(STRING_ELT(ALT,v))[0];
        if (v1 == 0 && v2 == 0) wf(w, r);
        else if (v1 == 1 && v2 == 1) wf(w, a);
        else wf(w, iupac_het(r, a));
      }
    }
    wf(w, '\n');
  }
  wf(w, "  ;\nEND;\n");
}


// ═══════════════════════════════════════════════════════════════════════════════
// XVI.  EIGENSTRAT  .geno writer
// ═══════════════════════════════════════════════════════════════════════════════

//' Write the sample-names header for an EIGENSTRAT .geno file (creates / truncates)
//'
//' EIGENSTRAT .geno is variant-major: one line per variant, samples as a
//' concatenated string of digits (no delimiter) in the same order as the .ind
//' file.  Encoding: 0 = hom-ref, 1 = het, 2 = hom-alt, 9 = missing.
//' This function creates (or truncates) the file; subsequent calls to
//' write_eigenstrat_chunk_cpp() append to it.
// [[Rcpp::export]]
void write_eigenstrat_geno_header_cpp(const std::string& out_file) {
  // Nothing to write as a header — just create / truncate the file so the
  // first chunk call can safely open it in append mode.
  WFile w(out_file, false);   // create/truncate only
  (void)w;   // immediately closed by destructor
}

//' Encode and append one chunk of variants to an EIGENSTRAT .geno file
//'
//' a1_mat / a2_mat: integer matrices (n_samples x n_chunk_var), rows = samples,
//' cols = variants.  Writes one line per variant: all samples concatenated
//' without any separator.
// [[Rcpp::export]]
void write_eigenstrat_chunk_cpp(const IntegerMatrix& a1_mat,
                                const IntegerMatrix& a2_mat,
                                const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, true);    // append

  for (int v = 0; v < nv; ++v) {
    for (int s = 0; s < ns; ++s) {
      const int v1 = a1_mat(s, v), v2 = a2_mat(s, v);
      if (v1 == NA_INTEGER || v2 == NA_INTEGER) wf(w, '9');
      else if (v1 == 0 && v2 == 0) wf(w, '0');
      else if (v1 == 1 && v2 == 1) wf(w, '2');
      else wf(w, '1');
    }
    wf(w, '\n');
  }
}


 
// ═══════════════════════════════════════════════════════════════════════════════
// XVII.  sNMF  .geno writer
// ═══════════════════════════════════════════════════════════════════════════════
 
//' Append one chunk of variants to an sNMF .geno file
//'
//' Opens \code{out_file} in append mode and writes one line per variant.
//' Sample genotype codes are concatenated with NO separator, e.g. "01120".
//' Encoding: 0 = hom-ref, 1 = het, 2 = hom-alt, 9 = missing.
//' No header line (sNMF format).
//'
//' The caller is responsible for creating / truncating the file before the
//' first chunk call (e.g. via \code{file.create(out_file)} in R).
//'
//' a1_mat / a2_mat: integer matrices (n_samples x n_chunk_var).
// [[Rcpp::export]]
void write_snmf_cpp(const IntegerMatrix& a1_mat,
                    const IntegerMatrix& a2_mat,
                    const std::string& out_file) {
  const int ns = a1_mat.nrow(), nv = a1_mat.ncol();
  WFile w(out_file, true);  // append — caller must truncate before first chunk
 
  for (int v = 0; v < nv; ++v) {
    for (int s = 0; s < ns; ++s) {
      const int v1 = a1_mat(s, v), v2 = a2_mat(s, v);
      if (v1 == NA_INTEGER || v2 == NA_INTEGER) wf(w, '9');
      else if (v1 == 0 && v2 == 0) wf(w, '0');
      else if (v1 == 1 && v2 == 1) wf(w, '2');
      else wf(w, '1');
    }
    wf(w, '\n');
  }
}
