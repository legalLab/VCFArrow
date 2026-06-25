#' @title vcf_extract_samples
#'
#' @description
#' Extract samples from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param samples -> individuals to retain/drop (character)
#' @param keep -> retain (TRUE) or drop (FALSE) individuals flag, default TRUE (Boolean)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report filtering stats, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes samples from a VCFArrow object,
#' returning a new VCFArrow object.
#' It uses the 'keep' flag to either keep or drop the samples in the list.
#' By default will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' By default will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_extract_samples(vcf_arrow = my_vcf, samples = my_samples, keep = TRUE, f_invar = TRUE, verbose = TRUE)
#' vcf_extract_samples(my_vcf, my_samples, TRUE, TRUE, TRUE)
#' vcf_extract_samples(my_vcf, my_samples)
#'
#' @export
#'

vcf_extract_samples <- function(vcf_arrow, samples, keep = TRUE, f_invar = TRUE, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  all_samples <- vcf_arrow@samples

  # keep only valid samples
  if (any(samples != "") && any(!(samples %in% all_samples))) {
    cli::cli_alert_info("Some samples do not exits in VCF")
  }

  # if empty → keep all - "keep everything" sentinel
  if (any(samples == "") || length(samples) == 0) {
    cli::cli_alert_info("No samples to keep or remove - keeping all samples")
    return(vcf_arrow)
  }

  # get final samples
  samples <- samples[samples %in% all_samples]

  # determine samples to keep
  keep_samples <- if (keep) samples else setdiff(all_samples, samples)

  # guard against silently producing a zero-sample VCFArrow.
  if (length(keep_samples) == 0L)
    cli::cli_abort(
      "All samples removed by vcf_extract_samples() — check {.arg samples} \\
       and {.arg keep}."
    )

  # apply filter using unified API
  if (length(keep_samples) > 1) {
    vcf_arrow <- .vcf_filter_columns(vcf_arrow, keep_samples, f_invar, verbose)
  } else {
    # f_invar forced FALSE with 1 sample: with a single individual, every
    # locus trivially has exactly one observed genotype state, so invariant
    # filtering would remove everything.
    vcf_arrow <- .vcf_filter_columns(vcf_arrow, keep_samples, f_invar = FALSE, verbose)
  }

  return(vcf_arrow)
}
