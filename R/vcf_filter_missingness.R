#' @title vcf_filter_missingness
#'
#' @description
#' Remove loci above missingness threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal missingness threshold, default 0.1 (numeric)
#' @param verbose -> flag to report total % missing data after filtering, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci above missingness threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' Missingness is a locus focused metric, i.e. missing data per locus.
#'
#' @examples
#' vcf_filter_missingness(vcf_arrow = my_vcf, threshold = p_miss)
#' vcf_filter_missingness(my_vcf, p_miss, verbose = FALSE)
#' vcf_filter_missingness(my_vcf, p_miss)
#'

vcf_filter_missingness <- function(vcf_arrow, threshold = 0.1, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # compute per-variant missingness and select passing variants
  miss_tbl <- vcf_arrow@gt |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      n_miss = sum(is.na(a1)),
      n_total = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      p_miss = n_miss / n_total
    ) |>
    dplyr::filter(p_miss <= !!threshold) |>
    dplyr::select(.row_id)

  # collect only row_ids
  keep <- miss_tbl |>
    dplyr::collect() |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  # optional reporting
  if (verbose) {

    total_stats <- vcf_arrow@gt |>
      dplyr::summarise(
        total_na = sum(is.na(a1)),
        total_n = dplyr::n()
      ) |>
      dplyr::collect()

    p_missing <- total_stats$total_na / total_stats$total_n

    cli::cli_alert_info("Final % missing data in VCF is {sprintf('%.2f', 100 * p_missing)}%")

  }

  return(vcf_arrow)
}
