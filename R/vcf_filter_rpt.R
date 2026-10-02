#' @title vcf_filter_rpt
#'
#' @description
#' Remove loci above REPEAT threshold from a VCFArrow object
#'
#' @author Tomas Hrbek September 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal REPEAT threshold, loci above it are removed,
#'   default 0.5 (numeric)
#' @param keep_na -> retain loci without a REPEAT value, default FALSE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci above REPEAT threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' REPEAT is calculated by disco_haplotypes in DiscoSNP-RAD (Gauthier et. al. 2020)
#' and registered as RPT in INFO.
#' REPEAT is calculated as max of overlap, depth score and is 1 if in a 
#' giant-component bubble, and is used for paralog detection -> high repeat 
#' values (>0.5) are indicative of paralogs.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_rpt(vcf, threshold = 0.5)
#'
#' @export
#'

vcf_filter_rpt <- function(vcf_arrow, threshold = 0.5, keep_na = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  rpt <- vcf_arrow@variants$RPT
  
  # check if field present
  if (all(is.na(rpt))) {
    cli::cli_alert_warning("RPT field not present; returning unfiltered VCFArrow")
    return(vcf_arrow)
  }

  cli::cli_alert_info("Applying REPEAT filter")

  # select passing variants
  keep <- if (keep_na) {
    is.na(rpt) | rpt <= threshold
  } else {
    !is.na(rpt) & rpt <= threshold
  }

  cli::cli_alert_info(
    "Retained {sum(keep, na.rm = TRUE)} / {idx$n_var} variant{?s} (Repeat <= {threshold})"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
