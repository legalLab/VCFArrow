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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_pass(vcf)
#'
#' @export
#'

vcf_filter_pass <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  cli::cli_alert_info("Applying PASS filter")

  # report if missing values in FILTER
  if (any(vcf_arrow@variants$FILTER == ".")) {
    cli::cli_alert_warning("Some FILTER values are NA. Defaulting to PASS")
  }

  # select passing variants
  keep <- vcf_arrow@variants$FILTER == "PASS" |
    vcf_arrow@variants$FILTER == "."

  cli::cli_alert_info(
    "Retained {sum(keep, na.rm = TRUE)} / {idx$n_var} variant{?s} (PASS)"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
