#' @title vcf_stats
#'
#' @description
#' Calculates basic stats of each samples from VCFArrow format data.
#' Includes average read depth per individual, missing data per individual,
#' Watterson's theta and pi.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param res_path -> directory where to write results
#' @param project -> base name of the project file
#' @param theta -> flag to perform theta and pi calculation, default FALSE (Boolean)
#'
#' @return table of statistics
#'
#' @details
#' This function calculates average read depth, heterozygosity
#' number of heterozygotes, number of reference and alternative homozygotes,
#' missing data and total number SNPs of each sample in an VCFArrow object.
#' Optionally calls vcf_theta() to get total and group Watterson's theta and pi.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_stats(vcf, res_path = tempdir(), project = "vaillantii",
#'           theta = TRUE)
#'
#' @export
#'

vcf_stats <- function(vcf_arrow, res_path, project, theta = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  if (any(is.na(vcf_arrow@groups)))
    cli::cli_alert_warning(
      "Missing group assignments. Use set_vcf_groups() to fill the 'groups' slot."
    )

  samples <- vcf_arrow@samples
  groups <- vcf_arrow@groups
  group_df <- tibble::tibble(sample = samples, group = groups)

  valid_row_ids <- vcf_arrow@variants$.row_id
  n_samples <- length(samples)

  ffiles <- list.files(vcf_arrow@path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(ffiles) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")
  chunk_nums <- as.integer(stringr::str_extract(basename(ffiles), "\\d+"))
  ffiles <- ffiles[order(chunk_nums)]

  # Per-sample accumulators (O(n_samples) memory — trivial) ────────────────────
  total_loci <- stats::setNames(integer(n_samples), samples)
  called_n <- stats::setNames(integer(n_samples), samples)
  het_n <- stats::setNames(integer(n_samples), samples)
  hom_ref_n <- stats::setNames(integer(n_samples), samples)
  hom_alt_n <- stats::setNames(integer(n_samples), samples)
  dp_n <- stats::setNames(integer(n_samples), samples)
  dp_sum <- stats::setNames(numeric(n_samples), samples)

  cli::cli_alert_info(
    "Computing per-sample stats: {length(valid_row_ids)} variant{?s} x \\
     {n_samples} sample{?s}, reading {length(ffiles)} chunk{?s} directly"
  )

  idx <- list(pos = .row_id_pos(valid_row_ids), samples = samples)

  parts <- .map_chunks(ffiles, .stats_chunk, shared = idx, label = "Scanning chunk")
  for (cnt in parts) {
    total_loci <- total_loci + cnt$total
    called_n <- called_n + cnt$called
    het_n <- het_n + cnt$het
    hom_ref_n <- hom_ref_n + cnt$hom_ref
    hom_alt_n <- hom_alt_n + cnt$hom_alt
    dp_n <- dp_n + cnt$dp_n
    dp_sum <- dp_sum + cnt$dp_sum
  }

  # Assemble per-sample stats (equivalent to the original collect() output) ────
  sample_stats <- tibble::tibble(
    sample = samples,
    read_depth = dp_sum[samples] / dp_n[samples],  # NaN if dp_n == 0, matches
    # mean(DP, na.rm=TRUE) on
    # an all-NA input
    homo_ref = hom_ref_n[samples],
    homo_alt = hom_alt_n[samples],
    hetero = het_n[samples],
    missing = total_loci[samples] - called_n[samples],
    non_missing = called_n[samples],
    total_loci = total_loci[samples]
  )
  sample_stats$heterozygosity <-
    sample_stats$hetero / (sample_stats$homo_ref + sample_stats$homo_alt + sample_stats$hetero)
  sample_stats$missing_p <- sample_stats$missing / sample_stats$total_loci

  # since calculation of theta and pi is computationally intensive
  # making it conditional
  if (theta) {
    # theta / pi
    theta_ <- vcf_theta(vcf_arrow)

    # final table
    out <- data.frame(
      sample = samples,
      read_depth = sample_stats$read_depth,
      heterozygosity = sample_stats$heterozygosity,
      heterozygotes = sample_stats$hetero,
      homozygotes = sample_stats$homo_ref + sample_stats$homo_alt,
      homozygotes_ref = sample_stats$homo_ref,
      homozygotes_alt = sample_stats$homo_alt,
      missing_p = sample_stats$missing_p,
      missing_n = sample_stats$missing,
      non_missing = sample_stats$non_missing,
      total_loci = sample_stats$total_loci,
      theta_total = theta_$theta_w,
      pi_total = theta_$pi
    )
    # explicit group join instead of relying on theta_$theta_g implicitly
    out <- dplyr::left_join(out, theta_$theta_g, by = "sample")
  } else {
    # final table
    out <- data.frame(
      sample = samples,
      read_depth = sample_stats$read_depth,
      heterozygosity = sample_stats$heterozygosity,
      heterozygotes = sample_stats$hetero,
      homozygotes = sample_stats$homo_ref + sample_stats$homo_alt,
      homozygotes_ref = sample_stats$homo_ref,
      homozygotes_alt = sample_stats$homo_alt,
      missing_p = sample_stats$missing_p,
      missing_n = sample_stats$missing,
      non_missing = sample_stats$non_missing,
      total_loci = sample_stats$total_loci
    )
    # since calculation of theta is skipped, add group information
    out <- dplyr::left_join(out, group_df, by = "sample")
  }

  out <- dplyr::arrange(out, group, sample)

  utils::write.table(
    out,
    file = file.path(res_path, paste0(project, "_stats.csv")),
    row.names = FALSE,
    quote = FALSE,
    sep = ","
  )

  invisible(vcf_arrow)
}

# Per-chunk, per-sample counts in one C++ pass (see src/sample_counts.cpp)
.stats_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2", "DP"), shared)
  sample_counts_cpp(ch$s, ch$a1, ch$a2, as.numeric(ch$DP), length(shared$samples))
}
