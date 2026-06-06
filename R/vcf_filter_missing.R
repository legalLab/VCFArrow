#' @title vcf_filter_missing
#'
#' @description
#' Remove samples with > % missing data from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal missing threshold, default 0.5 (numeric)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report filtering stats, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes samples from a VCFArrow object if they have
#' above threshold missing loci, returning a new VCFArrow object.
#' By default will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' By default will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_filter_missing(vcf_arrow = my_vcf, threshold = my_threshold)
#' vcf_filter_missing(my_vcf, my_threshold)
#' vcf_filter_missing(my_vcf)
#'

vcf_filter_missing <- function(vcf_arrow, threshold = 0.5,
                               f_invar = TRUE, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # compute missingness per sample
  samp_tbl <- vcf_arrow@gt |>
    dplyr::group_by(sample) |>
    dplyr::summarise(
      p_miss = mean(is.na(a1)),
      .groups = "drop"
    ) |>
    dplyr::collect()

  # samples to keep
  keep_samples <- samp_tbl$sample[samp_tbl$p_miss < threshold]

  if (length(keep_samples) == 0) {
    cli::cli_abort("All samples removed by filter")
  }

  # apply filter using unified API
  if (length(keep_samples) > 1) {
    vcf_arrow <- vcf_filter_columns(vcf_arrow, keep_samples, f_invar, verbose)
  } else {
    vcf_arrow <- vcf_filter_columns(vcf_arrow, keep_samples, f_invar = FALSE, verbose)
  }

  return(vcf_arrow)
}
