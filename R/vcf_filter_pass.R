#' @title vcf_filter_pass
#'
#' @description
#' Remove loci that did not PASS FILTER in a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci that did not PASS the FILTER in a VCFArrow object,
#' returning a new VCFArrow object.
#' PASS indicates that a variant has successfully passed all applied
#' quality control filters.
#' When no filters were applied and the FILTER field is '.' (Dot),
#' vcf_filter_pass() defaults to passing the locus.
#'
#' @examples
#' vcf_filter_pass(vcf_arrow = my_vcf)
#' vcf_filter_pass(my_vcf)
#'
#' @export
#'

vcf_filter_pass <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  # report if missing values in FILTER
  if (any(vcf_arrow@variants$FILTER == ".")) {
    cli::cli_alert_warning("Some FILTER values are NA. Defaulting to PASS")
  }

  # select passing variants
  keep <- vcf_arrow@variants$FILTER == "PASS" |
    vcf_arrow@variants$FILTER == "."

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variants (PASS)"
  )

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
