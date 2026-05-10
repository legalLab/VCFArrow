#' @title vcf_filter_rank
#'
#' @description
#' Remove loci below RANK threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal rank threshold, default 0.4 (numeric)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci below RANK threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' RANK is calculated in DiscoSNP-RAD (Gauthier et. al. 2020) and
#' registered as Pk in INFO.
#' RANK is calculated as sqrt(chi-sqr/n) of allele read counts, and
#' used for paralog detection -> very low rank values (<0.4) are
#' indicative of paralogs.
#'
#' @examples
#' vcf_filter_rank(vcf_arrow = my_vcf, rank = 0.4)
#' vcf_filter_rank(my_vcf, 0.4)
#' vcf_filter_rank(my_vcf)
#'

vcf_filter_rank <- function(vcf_arrow, threshold = 0.4, keep_na = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  rk <- vcf_arrow@variants$Rk

  # select passing variants
  if (all(is.na(rk))) {
    cli::cli_alert_warning("Rk field not present; returning unfiltered VCFArrow")
    return(vcf_arrow)
  }
  keep <- if (keep_na) {
    is.na(rk) | rk >= threshold
  } else {
    !is.na(rk) & rk >= threshold
  }

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
