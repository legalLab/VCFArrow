#' @title vcf_filter_quality
#'
#' @description
#' Remove loci below QUAL threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> QUAL threshold, default 30 (integer)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci below QUAL threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' QUAL is a Phred-scaled quality metric assessing whether any variant exists
#' at that position across all samples.
#' QUAL is different from GQ which represents the confidence in a
#' genotype call of a specific sample.
#' When QUAL not reported '.' (Dot), QUAL is assumed to be above threshold.
#'
#' @examples
#' vcf_filter_quality(vcf = my_vcf, threshold = 30)
#' vcf_filter_quality(my_vcf, 30)
#' vcf_filter_quality(my_vcf)
#'

vcf_filter_quality <- function(vcf_arrow, threshold = 30) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  qual <- suppressWarnings(as.numeric(vcf_arrow@variants$QUAL))

  # select passing variants
  if (any(is.na(qual))) {
    qual[is.na(qual)] <- threshold
    cli::cli_alert_warning("Some QUAL values are NA. Replacing with QUAL = {threshold} (threshold parameter)")
  }

  keep <- qual >= threshold

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
