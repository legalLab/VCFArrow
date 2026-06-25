#' @noRd

# vcf_filters.R
#
# ── Shared chunk-reading helpers ───────────────────────────────────────────────

.get_sorted_feather_files <- function(path) {
  files <- list.files(path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(files) == 0L) cli::cli_abort("No .arrow files found in {path}")
  nums <- as.integer(stringr::str_extract(basename(files), "\\d+"))
  files[order(nums)]
}

# Build the two fast-lookup structures used by every filter's chunk loop.
# Returns a list with $lv (logical vector for row_id membership) and
# $col_idx (named integer for row_id → 1-based position).

.vcf_filter_index <- function(vcf_arrow) {
  valid <- vcf_arrow@variants$.row_id
  max_id <- max(valid)
  lv <- logical(max_id + 1L)
  lv[valid] <- TRUE
  list(
    lv = lv,
    col_idx = setNames(seq_along(valid), as.character(valid)),
    n_var = length(valid),
    samples = vcf_arrow@samples
  )
}

# ── Shared chunk-reading helper ────────────────────────────────────────────────
#
# Reads every feather file directly, filtering each chunk to the VCFArrow's
# CURRENT variant set (vcf_arrow@variants$.row_id) and CURRENT sample list
# (vcf_arrow@samples) — the same filtering every exporter in this package
# applies, and the same fix that was needed in write_vcf() to avoid emitting
# already-filtered-out variants.
#
# Accumulates missingness counts in O(n_samples) memory.  Optionally also
# collects a bounded random subsample of (sample, DP) pairs for the coverage
# violin plot, controlled by max_points_per_sample so memory stays bounded
# regardless of dataset size.

.scan_vcf_gt <- function(vcf_arrow, collect_dp = FALSE,
                         max_points_per_sample = 5000L) {

  samples <- vcf_arrow@samples
  valid_row_ids <- vcf_arrow@variants$.row_id
  n_var <- length(valid_row_ids)
  n_samples <- length(samples)

  feather_files <- list.files(vcf_arrow@path, pattern = "\\.arrow$",
                              full.names = TRUE)
  if (length(feather_files) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")
  chunk_nums <- as.integer(stringr::str_extract(basename(feather_files), "\\d+"))
  feather_files <- feather_files[order(chunk_nums)]

  total_loci <- stats::setNames(integer(n_samples), samples)
  missing_n <- stats::setNames(integer(n_samples), samples)

  # Probability of keeping a row for the DP subsample: each sample has
  # exactly n_var rows total in the dense gt_long format, so a flat
  # per-row keep-probability of max_points_per_sample / n_var yields
  # ~max_points_per_sample retained rows per sample, independent of which
  # chunk they fall in.
  keep_p <- if (collect_dp) min(1, max_points_per_sample / max(n_var, 1L)) else 0
  dp_parts <- if (collect_dp) vector("list", length(feather_files)) else NULL

  cols <- if (collect_dp) c(".row_id", "sample", "a1", "a2", "DP")
  else c(".row_id", "sample", "a1", "a2")

  cli::cli_progress_bar("Scanning chunk", total = length(feather_files))

  for (i in seq_along(feather_files)) {
    chunk <- arrow::read_feather(feather_files[[i]], dplyr::all_of(cols))
    chunk <- chunk[chunk$.row_id %in% valid_row_ids &
                     chunk$sample %in% samples, , drop = FALSE]

    if (nrow(chunk) > 0L) {
      # Per-sample totals
      tt <- table(chunk$sample)
      total_loci[names(tt)] <- total_loci[names(tt)] + as.integer(tt)

      # Missing = either allele NA (consistent with .pop_counts_from_chunk()
      # elsewhere in this package, which treats a1/a2 symmetrically)
      miss_mask <- is.na(chunk$a1) | is.na(chunk$a2)
      if (any(miss_mask)) {
        mt <- table(chunk$sample[miss_mask])
        missing_n[names(mt)] <- missing_n[names(mt)] + as.integer(mt)
      }

      if (collect_dp) {
        dp_ok <- !is.na(chunk$DP)
        if (any(dp_ok)) {
          sub <- chunk[dp_ok, c("sample", "DP")]
          if (keep_p < 1) sub <- sub[stats::runif(nrow(sub)) < keep_p, , drop = FALSE]
          if (nrow(sub) > 0L) dp_parts[[i]] <- sub
        }
      }
    }

    chunk <- NULL   # release Arrow Table reference before next iteration
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  list(
    total_loci = total_loci,
    missing_n = missing_n,
    dp_df = if (collect_dp) do.call(rbind, Filter(Negate(is.null), dp_parts))
    else NULL
  )
}
