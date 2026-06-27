# vcf2formats.R
#
# Requires: vcf_wide.cpp compiled via Rcpp::sourceCpp("vcf_wide.cpp")
#           or placed in src/ of a package.
#
# Architecture
# ──────────────────────────────────────────────────────────────────────────────
# All functions bypass the Arrow query engine and read feather chunks directly
# with arrow::read_feather() — the same approach as write_vcf().
#
# Memory modes
#   Chunk-by-chunk O(chunk):               SmartSNP, Treemix*, fineRADstructure
#   Accumulate individuals O(S × V × 4B):  Structure, Arlequin, Genepop,
#                                          BayesAss, Genlight, Related,
#                                          Apparent, FASTA, Nexus, SNAPP
#   Accumulate pop counts  O(P × V × 4B):  BayesScan, Migrate, Treemix*
#   (* Treemix uses chunk path when inc_missing = TRUE, pop-count path otherwise)
#
# S = n_samples, V = n_var, P = n_pops
# For the example dataset (80 × 1.57M):  O(S×V) ≈ 1 GB per integer matrix


# ── Shared setup ──────────────────────────────────────────────────────────────
#
# Returns a named list:
#   $samples        character       sample names, group-contiguous order
#   $samples_groups character       group label per sample (parallel to $samples)
#   $group_ids      integer         1-based group id per sample
#   $group_names    character       ordered group labels  (== keep_groups)
#   $group_sizes    integer         samples per group
#   $n_pops         integer
#   $n_samples      integer
#   $P              integer matrix  n_pops × n_samples  (pop-membership 0/1)
#   $valid_row_ids  integer         .row_id values surviving biallelic/non-indel
#   $var_col_index  named int       .row_id (as char) → 1-based column position
#   $variants       data.frame      filtered, arranged by .row_id  (REF, ALT, …)
#   $loci           character       variant IDs or "CHROM_POS"
#   $feather_files  character       chunk paths in natural numeric order
#   $n_var          integer

.vcf_export_setup <- function(vcf_arrow, keep_groups) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  all_samples <- vcf_arrow@samples
  all_groups <- vcf_arrow@groups

  if (is.null(keep_groups)) keep_groups <- unique(all_groups)

  if (any(is.na(all_groups)))
    cli::cli_abort("Some samples have no group assignment – use set_vcf_groups()")
  if (any(!(keep_groups %in% all_groups)))
    cli::cli_abort("Some requested groups do not exist in the VCFArrow object")

  # Samples in group-contiguous order, respecting keep_groups ordering
  keep_mask <- all_groups %in% keep_groups
  samples <- all_samples[keep_mask]
  samples_groups <- all_groups[keep_mask]
  grp_order <- order(match(samples_groups, keep_groups))
  samples <- samples[grp_order]
  samples_groups <- samples_groups[grp_order]

  group_names <- keep_groups
  group_sizes <- as.integer(table(factor(samples_groups, levels = keep_groups)))
  group_ids <- as.integer(factor(samples_groups, levels = keep_groups))
  n_pops <- length(group_names)
  n_samples <- length(samples)

  # Pop-membership matrix  P[p, s] = 1 iff sample s belongs to pop p
  P <- matrix(0L, nrow = n_pops, ncol = n_samples)
  for (p in seq_len(n_pops)) {
    P[p, group_ids == p] <- 1L
    }

  # Variant filter: biallelic SNPs only
  valid_row_ids <- vcf_arrow@variants |>
    dplyr::filter(is_biallelic, !is_indel) |>
    dplyr::arrange(.row_id) |>
    dplyr::pull(.row_id)

  var_col_index <- setNames(seq_along(valid_row_ids), as.character(valid_row_ids))

  variants <- vcf_arrow@variants |>
    dplyr::filter(.row_id %in% valid_row_ids) |>
    dplyr::arrange(.row_id)

  loci <- if (!all(is.na(variants$ID))) {
    variants$ID
  } else {
    paste0(variants$CHROM, "_", variants$POS)
  }

  # Feather files in natural numeric chunk order
  feather_files <- list.files(vcf_arrow@path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(feather_files) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")
  chunk_nums <- as.integer(stringr::str_extract(basename(feather_files), "\\d+"))
  feather_files <- feather_files[order(chunk_nums)]

  list(
    samples = samples,
    samples_groups = samples_groups,
    group_ids = group_ids,
    group_names = group_names,
    group_sizes = group_sizes,
    n_pops = n_pops,
    n_samples = n_samples,
    P = P,
    valid_row_ids = valid_row_ids,
    var_col_index = var_col_index,
    variants = variants, loci = loci,
    feather_files = feather_files,
    n_var = length(valid_row_ids)
  )
}


# ── Chunk reshape: individual level ───────────────────────────────────────────
# Returns list($a1, $a2: n_samples × n_chunk_var integer matrices,
#              $col_idx: global column indices for this chunk)
# or NULL if no valid rows in this chunk.

.reshape_chunk <- function(chunk, setup) {
  chunk <- chunk[chunk$.row_id %in% setup$valid_row_ids & chunk$sample %in% setup$samples, ]
  if (nrow(chunk) == 0L) return(NULL)
  sample_order <- match(chunk$sample, setup$samples)
  chunk <- chunk[order(chunk$.row_id, sample_order), ]
  chunk_row_ids <- unique(chunk$.row_id)
  n_chunk_var <- length(chunk_row_ids)
  col_idx <- unname(setup$var_col_index[as.character(chunk_row_ids)])
  list(
    a1 = matrix(chunk$a1, nrow = setup$n_samples, ncol = n_chunk_var),
    a2 = matrix(chunk$a2, nrow = setup$n_samples, ncol = n_chunk_var),
    col_idx = col_idx, chunk_row_ids = chunk_row_ids
  )
}


# ── Chunk aggregation: population level ───────────────────────────────────────
# Returns $ref, $alt, $nobs: n_pops × n_chunk_var integer matrices.

.pop_counts_from_chunk <- function(rc, setup) {
  a1 <- rc$a1; a2 <- rc$a2

  # NA-safe logical matrices
  a1_ref <- (!is.na(a1)) & (a1 == 0L)
  a2_ref <- (!is.na(a2)) & (a2 == 0L)
  a1_alt <- (!is.na(a1)) & (a1 == 1L)
  a2_alt <- (!is.na(a2)) & (a2 == 1L)
  a1_obs <- (!is.na(a1))
  a2_obs <- (!is.na(a2))

  list(
    ref = setup$P %*% (a1_ref + a2_ref),
    alt = setup$P %*% (a1_alt + a2_alt),
    nobs = setup$P %*% (a1_obs + a2_obs)  # FIXED: was !is.na(a1) + !is.na(a2)
  )
}

# ── Complete-variant mask (for inc_missing = FALSE) ───────────────────────────

.complete_var_mask <- function(nobs_full, setup) {
  colSums(nobs_full) == 2L * setup$n_samples
}


# ── Raw-byte accumulation helpers ─────────────────────────────────────────────
#
# Encode a1/a2 integer values (0L, 1L, NA_integer_) into raw bytes:
#   0L         → 0x00
#   1L         → 0x01
#   NA_integer_→ 0xFF
#
# This gives 4× memory reduction vs. integer matrices:
#   80 × 1,570,876 × 1 byte × 2 = 251 MB   (vs. ~1 GB with integer)
#
# The C++ writers receive IntegerMatrix, so .raw_to_int_mat() decodes on the
# way out.  Decoding one chunk at a time keeps peak memory low.

.alloc_raw_matrix <- function(nrow, ncol) {
  # rawMatrix is not a base type; use a raw vector and track dims ourselves.
  m <- raw(nrow * ncol)          # all 0x00 by default
  m[seq_len(nrow * ncol)] <- as.raw(0xFF)   # initialise to NA marker (0xFF)
  attr(m, "dim") <- c(nrow, ncol)
  m
}

# Write a chunk's integer matrices into raw accumulation buffers in-place.
# col_idx: 1-based column positions (variant columns for this chunk).

.fill_raw_matrix <- function(raw_mat, int_mat, col_idx, nrow) {
  # int_mat is nrow × n_chunk_var integer matrix.
  # Map:  0L → 0x00, 1L → 0x01, NA → 0xFF
  encoded <- int_mat
  encoded[is.na(encoded)] <- 0xFFL
  storage.mode(encoded) <- "integer"
  # Column-major fill: column v occupies rows [(col-1)*nrow + 1 .. col*nrow]
  for (k in seq_along(col_idx)) {
    base <- (col_idx[k] - 1L) * nrow
    raw_mat[base + seq_len(nrow)] <- as.raw(encoded[, k])
  }
  raw_mat   # returned (modified copy; R semantics)
}

# Decode a raw accumulation matrix to an integer matrix for a C++ writer.
# Converts 0xFF back to NA_integer_.

.raw_to_int_mat <- function(raw_mat) {
  int_mat <- as.integer(raw_mat)
  int_mat[int_mat == 255L] <- NA_integer_
  dim(int_mat) <- dim(raw_mat)
  int_mat
}


# ── Memory-aware accumulation loop ────────────────────────────────────────────
#
# Drop-in replacement for .accumulate_individuals() that uses raw byte storage
# and calls gc() after each chunk to hint Arrow to release pool memory.
# The gc() call does not guarantee OS page return with jemalloc, but with the
# system allocator (set in VCFArrow-lifecycle) it does.

.accumulate_individuals <- function(setup, label) {
  nrow <- setup$n_samples
  ncol <- setup$n_var

  a1_raw <- .alloc_raw_matrix(nrow, ncol)
  a2_raw <- .alloc_raw_matrix(nrow, ncol)

  cli::cli_alert_info(
    "Accumulating {label}: {ncol} variants x {nrow} samples \\
     ({.strong {format(round(2 * nrow * ncol / 1024^2), big.mark=',')}} MiB raw storage)"
  )
  cli::cli_progress_bar("Reading chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    chunk <- NULL          # drop Arrow Table reference before gc()
    gc(verbose = FALSE, full = FALSE)   # with system allocator: pages returned to OS

    if (!is.null(rc)) {
      a1_raw <- .fill_raw_matrix(a1_raw, rc$a1, rc$col_idx, nrow)
      a2_raw <- .fill_raw_matrix(a2_raw, rc$a2, rc$col_idx, nrow)
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # Decode raw → integer only when handing off to C++.
  # Both matrices are decoded at the same time, so peak overhead is
  # 2 × nrow × ncol × 4 bytes = the same as before, but only briefly.
  list(a1 = .raw_to_int_mat(a1_raw),
       a2 = .raw_to_int_mat(a2_raw))
}


# ── Population-count accumulation: lowmem (uint16) variant ────────────────────
#
# Why uint16 instead of the 1-byte raw used for individuals
# ───────────────────────────────────────────────────────────────────────────
# Allele matrices only ever hold {0, 1, NA} — three states, one byte is ample.
# Pop-count matrices hold SUMS: ref/alt/nobs range from 0 up to 2 × group_size
# (every sample in the population contributing 2 alleles).  A single raw byte
# overflows silently above 255, i.e. above ~127 samples per population — too
# easy to hit in real datasets (cohort studies, large cohorts) to risk.
# Two bytes (uint16, max 65,535) covers any realistic population size while
# still halving memory vs. a 4-byte integer matrix.
#
# Relative benefit vs. .accumulate_individuals_lowmem
# ───────────────────────────────────────────────────────────────────────────
# n_pops is normally far smaller than n_samples (tens vs. hundreds+), so the
# pop-count matrices (n_pops × n_var) are already much smaller than the
# individual-level ones (n_samples × n_var).  For the example dataset
# (80 samples, ~5-10 pops typical), pop matrices are already only ~10-20% the
# size of the individual matrices.  This function matters most for datasets
# with very large n_var (tens of millions of variants) or unusually many
# populations, where even the smaller pop-count matrices start to add up.
#
# Encoding: two parallel raw vectors (lo, hi) per matrix; value = lo + hi*256.
# Values are clamped to 65535 (with a one-time warning) rather than wrapping,
# since silent wraparound would corrupt downstream allele-frequency math.

.alloc_u16_pair <- function(n) {
  list(lo = raw(n), hi = raw(n))   # both zero-initialised; counts start at 0
}

.encode_u16 <- function(int_vec) {
  if (any(int_vec > 65535L, na.rm = TRUE)) {
    cli::cli_warn(
      "Population allele count exceeds 65,535 for at least one cell; \\
       clamping. This indicates an unusually large population \\
       ({.code > 32,767} samples) — verify {.arg group_sizes}."
    )
    int_vec <- pmin(int_vec, 65535L)
  }
  list(lo = as.raw(int_vec %% 256L),
       hi = as.raw(int_vec %/% 256L))
}

.decode_u16 <- function(lo, hi) {
  as.integer(lo) + as.integer(hi) * 256L
}

# Add a chunk's per-population counts into the running u16 accumulator for
# the relevant columns (col_idx).  Each column is read, added, and re-written
# — true accumulation, matching .accumulate_pops()'s "+=" semantics (needed
# in case the same variant's rows are ever split across chunks, e.g. after
# vcf_bind_sparse()).

.add_u16_columns <- function(pair, add_mat, col_idx, n_pops) {
  for (k in seq_along(col_idx)) {
    base <- (col_idx[k] - 1L) * n_pops
    idx <- base + seq_len(n_pops)
    cur <- .decode_u16(pair$lo[idx], pair$hi[idx])
    new <- cur + add_mat[, k]
    enc <- .encode_u16(new)
    pair$lo[idx] <- enc$lo
    pair$hi[idx] <- enc$hi
  }
  pair
}

# Drop-in low-memory replacement for .accumulate_pops()
#
# Uses uint16 (2-byte) raw-pair storage for ref/alt/nobs during accumulation,
# decoding to integer matrices only at the very end (the format C++ writers
# require).  Returns the same list(ref=, alt=, nobs=) structure.

.accumulate_pops_lowmem <- function(setup, label) {
  n_pops <- setup$n_pops
  n_var  <- setup$n_var

  ref_pair <- .alloc_u16_pair(n_pops * n_var)
  alt_pair <- .alloc_u16_pair(n_pops * n_var)
  nobs_pair <- .alloc_u16_pair(n_pops * n_var)

  cli::cli_alert_info(
    "Accumulating {label}: {n_var} variants x {n_pops} pops \\
     ({.strong {format(round(6 * n_pops * n_var / 1024^2), big.mark=',')}} MiB \\
     raw storage, vs {format(round(12 * n_pops * n_var / 1024^2), big.mark=',')} \\
     MiB with integer matrices)"
  )
  cli::cli_progress_bar("Reading chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    chunk <- NULL
    gc(verbose = FALSE, full = FALSE)

    if (!is.null(rc)) {
      pc <- .pop_counts_from_chunk(rc, setup)
      ref_pair  <- .add_u16_columns(ref_pair,  pc$ref,  rc$col_idx, n_pops)
      alt_pair  <- .add_u16_columns(alt_pair,  pc$alt,  rc$col_idx, n_pops)
      nobs_pair <- .add_u16_columns(nobs_pair, pc$nobs, rc$col_idx, n_pops)
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # Decode to integer matrices for the C++ writers.
  ref  <- matrix(.decode_u16(ref_pair$lo,  ref_pair$hi),  nrow = n_pops, ncol = n_var)
  alt  <- matrix(.decode_u16(alt_pair$lo,  alt_pair$hi),  nrow = n_pops, ncol = n_var)
  nobs <- matrix(.decode_u16(nobs_pair$lo, nobs_pair$hi), nrow = n_pops, ncol = n_var)

  list(ref = ref, alt = alt, nobs = nobs)
}

# ── Accumulation loop (shared boilerplate for population-level formats) ────────
# Returns list($ref, $alt, $nobs): n_pops × n_var integer matrices.

.accumulate_pops <- function(setup, label) {
  ref <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)
  alt <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)
  nobs <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)
  cli::cli_alert_info(
    "Accumulating {label}: {setup$n_var} variants x {setup$n_pops} pops \\
     ({.strong {format(round(12 * n_pops * n_var / 1024^2), big.mark=',')}} MiB \\
     raw storage)"
  )
  cli::cli_progress_bar("Reading chunk", total = length(setup$feather_files))
  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc)) {
      pc <- .pop_counts_from_chunk(rc, setup)
      ref[, rc$col_idx] <- ref[, rc$col_idx] + pc$ref
      alt[, rc$col_idx] <- alt[, rc$col_idx] + pc$alt
      nobs[, rc$col_idx] <- nobs[, rc$col_idx] + pc$nobs
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  list(ref = ref, alt = alt, nobs = nobs)
}
