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

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # ensure we have group info
  if (any(is.na(vcf_arrow@groups))) {
    cli::cli_alert_warning("Missing group assignments. Use set_vcf_groups() to fill the 'groups' slot.")
  }

  gt <- vcf_arrow@gt
  samples <- vcf_arrow@samples
  groups <- vcf_arrow@groups

  # per-sample stats
  sample_stats <- gt |>
    dplyr::mutate(
      called = !(is.na(a1) | is.na(a2)),
      is_het = (a1 != a2) & called,
      is_hom_ref = (a1 == 0 & a2 == 0),
      is_hom_alt = (a1 == 1 & a2 == 1)
    ) |>
    dplyr::group_by(sample) |>
    dplyr::summarise(
      read_depth = mean(DP, na.rm = TRUE),
      homo_ref = sum(is_hom_ref, na.rm = TRUE),
      homo_alt = sum(is_hom_alt, na.rm = TRUE),
      hetero = sum(is_het, na.rm = TRUE),
      miss = sum(!called),
      non_missing = sum(called),
      total_loci = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::collect()

  # derived stats
  sample_stats$heterozygosity <-
    sample_stats$hetero /
    (sample_stats$homo_ref + sample_stats$homo_alt + sample_stats$hetero)

  sample_stats$missing_p <-
    sample_stats$miss / sample_stats$total_loci

  sample_stats <- sample_stats[match(samples, sample_stats$sample), ]

  # theta (fully Arrow-native upstream)
  theta <- vcf_theta(gt, samples, groups)

  # final table
  stats <- data.frame(
    sample = samples,
    read_depth = sample_stats$read_depth,
    heterozygosity = sample_stats$heterozygosity,
    heterozygotes = sample_stats$hetero,
    homozygotes = sample_stats$homo_ref + sample_stats$homo_alt,
    homozygotes_ref = sample_stats$homo_ref,
    homozygotes_alt = sample_stats$homo_alt,
    missing_p = sample_stats$missing_p,
    missing = sample_stats$miss,
    non_missing = sample_stats$non_missing,
    total_loci = sample_stats$total_loci,
    theta_total = theta$theta_w,
    pi_total = theta$pi
  )

  stats <- cbind(stats, theta$theta_g)

  write.table(
    stats,
    file = paste0(res_path, project, "_stats.csv"),
    row.names = FALSE,
    quote = FALSE,
    sep = ","
  )

  invisible(vcf_arrow)
}
