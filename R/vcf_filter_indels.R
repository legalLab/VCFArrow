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
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_quality(vcf_arrow = my_vcf)
#' vcf_filter_quality(my_vcf)
#'
#' @export
#'

vcf_filter_indels <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  # select passing variants
  keep <- !vcf_arrow@variants$is_indel

  cli::cli_alert_info(
    "Removed {length(keep)} / {idx$n_var} variants (Indels)"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
