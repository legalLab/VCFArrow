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
    samples = samples, samples_groups = samples_groups,
    group_ids = group_ids, group_names = group_names,
    group_sizes = group_sizes, n_pops = n_pops, n_samples = n_samples,
    P = P, valid_row_ids = valid_row_ids, var_col_index = var_col_index,
    variants = variants, loci = loci, feather_files = feather_files,
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
  list(
    ref = setup$P %*% ((!is.na(a1) & a1 == 0L) + (!is.na(a2) & a2 == 0L)),
    alt = setup$P %*% ((!is.na(a1) & a1 == 1L) + (!is.na(a2) & a2 == 1L)),
    nobs = setup$P %*% (!is.na(a1) + !is.na(a2))
  )
}


# ── Complete-variant mask (for inc_missing = FALSE) ───────────────────────────
.complete_var_mask <- function(nobs_full, setup) {
  colSums(nobs_full) == 2L * setup$n_samples
}


# ── Accumulation loop (shared boilerplate for individual-level formats) ───────
# Reads all feather chunks and fills pre-allocated a1_full / a2_full matrices.
# Returns list($a1, $a2).

.accumulate_individuals <- function(setup, label) {
  a1 <- matrix(NA_integer_, nrow = setup$n_samples, ncol = setup$n_var)
  a2 <- matrix(NA_integer_, nrow = setup$n_samples, ncol = setup$n_var)
  cli::cli_alert_info(
    "Accumulating {label}: {setup$n_var} variants x {setup$n_samples} samples"
  )
  cli::cli_progress_bar("Reading chunk", total = length(setup$feather_files))
  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc)) {
      a1[, rc$col_idx] <- rc$a1
      a2[, rc$col_idx] <- rc$a2
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  list(a1 = a1, a2 = a2)
}


# ── Accumulation loop (shared boilerplate for population-level formats) ───────
# Returns list($ref, $alt, $nobs): n_pops × n_var integer matrices.

.accumulate_pops <- function(setup, label) {
  ref  <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)
  alt  <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)
  nobs <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)
  cli::cli_alert_info(
    "Accumulating {label}: {setup$n_var} variants x {setup$n_pops} pops"
  )
  cli::cli_progress_bar("Reading chunk", total = length(setup$feather_files))
  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))
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
