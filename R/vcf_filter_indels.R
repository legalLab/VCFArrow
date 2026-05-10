#' @title vcf_filter_rank
#'
#' @description
#' Remove indels from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes indels from a VCFArrow object,
#' returning a new vcfR object.
#'
#' @examples
#' vcf_filter_quality(vcf_arrow = my_vcf)
#' vcf_filter_quality(my_vcf)
#'

vcf_filter_indels <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # select passing variants
  keep <- !vcf_arrow@variants$is_indel

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
