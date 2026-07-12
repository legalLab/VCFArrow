#' @title vcf_filter_indels
#'
#' @description
#' Remove non-indels from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes non-indel loci from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_indels(vcf_arrow = my_vcf)
#' vcf_filter_indels(my_vcf)
#'
#' @export
#'

vcf_filter_indels <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  # select passing variants
  keep <- !vcf_arrow@variants$is_indel

  cli::cli_alert_info("Applying indel filter")

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variants (non-Indels)"
  )

  return(vcf_arrow)
}
