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
#'
#' @return table of statistics
#'
#' @details
#' This function calculates average read depth, heterozygosity
#' number of heterozygotes, number of reference and alternative homozygotes,
#' missing data and total number SNPs of each sample in an VCFArrow object.
#' It calls vcf_theta() to get total and group Watterson's theta and pi.
#'
#' @examples
#' vcf_stats(vcf_arrow = my_vcf, res_path = my_res_path, project = my_project)
#' vcf_stats(my_vcf, my_res_path, my_project)
#'

vcf_stats <- function(vcf_arrow, res_path, project) {

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

  feather_files <- list.files(vcf_arrow@path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(feather_files) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")
  chunk_nums <- as.integer(stringr::str_extract(basename(feather_files), "\\d+"))
  feather_files <- feather_files[order(chunk_nums)]

  # Per-sample accumulators (O(n_samples) memory — trivial) ────────────────────
  total_loci <- stats::setNames(integer(n_samples), samples)
  called_n <- stats::setNames(integer(n_samples), samples)
  het_n <- stats::setNames(integer(n_samples), samples)
  hom_ref_n <- stats::setNames(integer(n_samples), samples)
  hom_alt_n <- stats::setNames(integer(n_samples), samples)
  dp_n <- stats::setNames(integer(n_samples), samples)
  dp_sum <- stats::setNames(numeric(n_samples), samples)

  cli::cli_alert_info(
    "Computing per-sample stats: {length(valid_row_ids)} variants x \\
     {n_samples} samples, reading {length(feather_files)} chunk(s) directly"
  )
  cli::cli_progress_bar("Scanning chunk", total = length(feather_files))

  for (fpath in feather_files) {

    chunk <- arrow::read_feather(
      fpath, col_select = c(".row_id", "sample", "a1", "a2", "DP")
    )
    chunk <- chunk[chunk$.row_id %in% valid_row_ids &
                     chunk$sample %in% samples, , drop = FALSE]

    if (nrow(chunk) > 0L) {

      called <- !(is.na(chunk$a1) | is.na(chunk$a2))
      is_het <- called & (chunk$a1 != chunk$a2)
      is_hom_ref <- called & (chunk$a1 == 0 & chunk$a2 == 0)
      is_hom_alt <- called & (chunk$a1 == 1 & chunk$a2 == 1)

      tt <- table(chunk$sample)
      total_loci[names(tt)] <- total_loci[names(tt)] + as.integer(tt)

      if (any(called)) {
        ct <- table(chunk$sample[called])
        called_n[names(ct)] <- called_n[names(ct)] + as.integer(ct)
      }
      if (any(is_het)) {
        ht <- table(chunk$sample[is_het])
        het_n[names(ht)] <- het_n[names(ht)] + as.integer(ht)
      }
      if (any(is_hom_ref)) {
        rt <- table(chunk$sample[is_hom_ref])
        hom_ref_n[names(rt)] <- hom_ref_n[names(rt)] + as.integer(rt)
      }
      if (any(is_hom_alt)) {
        at <- table(chunk$sample[is_hom_alt])
        hom_alt_n[names(at)] <- hom_alt_n[names(at)] + as.integer(at)
      }

      dp_ok <- !is.na(chunk$DP)
      if (any(dp_ok)) {
        dpt <- table(chunk$sample[dp_ok])
        dp_n[names(dpt)] <- dp_n[names(dpt)] + as.integer(dpt)
        dps <- tapply(chunk$DP[dp_ok], chunk$sample[dp_ok], sum)
        dp_sum[names(dps)] <- dp_sum[names(dps)] + as.numeric(dps)
      }
    }

    chunk <- NULL  # release Arrow Table reference before next chunk
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

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

  # theta / pi
  theta <- vcf_theta(vcf_arrow)

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
    theta_total = theta$theta_w,
    pi_total = theta$pi
  )

  out <- cbind(out, theta$theta_g)

  # explicit group join instead of relying on theta$theta_g implicitly
  if (!"group" %in% names(out)) {
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
