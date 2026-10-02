#' @title vcf_filter_paralogs
#'
#' @description
#' Remove loci above REPEAT and below RANK threshold from a VCFArrow object
#'
#' @author Tomas Hrbek September 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold_rk -> decimal RANK threshold, loci below it are removed,
#'   default 0.2 (numeric)
#' @param threshold_rpt -> decimal REPEAT threshold, loci above it are removed,
#'   default 0.5 (numeric)
#' @param keep_na -> retain loci without a RANK or REPEAT value, default FALSE
#'   (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci below the RANK threshold and above the REPEAT
#' threshold from a VCFArrow object, returning a new VCFArrow object.
#' It applies vcf_filter_rank() and then vcf_filter_rpt().
#' RANK and REPEAT metrics focus on different paralog signals, and used
#' jointly maximize true positive and minimize false positive detections.
#' Simulations indicate that removing loci with Rk < 0.2 or RPT > 0.5 gives the
#' best results, with Rk < 0.4 or RPT > 0.5 being slightly more conservative.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_paralogs(vcf, threshold_rk = 0.2, threshold_rpt = 0.5)
#'
#' @export
#'

vcf_filter_paralogs <- function(vcf_arrow, threshold_rk = 0.2, 
                            threshold_rpt = 0.5, keep_na = FALSE) {
  
  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  
  rk <- vcf_arrow@variants$Rk
  rpt <- vcf_arrow@variants$RPT
  
  # check if field present
  if (all(is.na(rk)) & all(is.na(rpt))) {
    cli::cli_alert_warning("Neither the Rk nor RPT fields are present; 
                           not applying the parallog filter; 
                           returning unfiltered VCFArrow")
    return(vcf_arrow)
  }

  # use existing filters
  vcf_arrow <- vcf_filter_rank(vcf_arrow, threshold_rk, keep_na) |>
    vcf_filter_rpt(threshold_rpt, keep_na)

  return(vcf_arrow)
}
