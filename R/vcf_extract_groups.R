#' @title vcf_extract_groups
#'
#' @description
#' Extract groups of samples from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param groups -> groups to retain/drop (character)
#' @param keep -> retain (TRUE) or drop (FALSE) individuals flag, default TRUE (Boolean)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report filtering stats, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes groups of samples from a VCFArrow object,
#' returning a new VCFArrow object.
#' It uses the 'keep' flag to either keep or drop the groups in the list.
#' By default will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' By default will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_extract_groups(vcf_arrow = my_vcf, groups = my_groups, keep = TRUE, f_invar = TRUE, verbose = TRUE)
#' vcf_extract_groups(my_vcf, my_groups, TRUE, TRUE, TRUE)
#' vcf_extract_groups(my_vcf, my_groups)
#'

vcf_extract_groups <- function(vcf_arrow, groups, keep = TRUE, f_invar = TRUE, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  all_samples <- vcf_arrow@samples
  all_groups <- vcf_arrow@groups

  # check if group assignments exist and groups are valid
  if (any(groups != "") && any(is.na(all_groups))) {
    cli::cli_abort("Some samples are not assigned to groups - use set_vcf_groups()")
  } else if (any(groups != "") && any(!(groups %in% all_groups))) {
    cli::cli_abort("Some selected groups do not exits in VCFArrow")
  }

  # if empty → keep all
  if (any(groups == "")) {
    cli::cli_alert_info("No groups of samples to keep or remove - keeping all samples")
    return(vcf_arrow)
  }

  # determine which samples belong to target populations
  if (keep) {
    keep_samples <- all_samples[all_groups %in% groups]
  } else {
    keep_samples <- all_samples[!(all_groups %in% groups)]
  }

  # apply filter using unified API
  if (length(keep_samples) > 1) {
    vcf_arrow <- vcf_filter_columns(vcf_arrow, keep_samples, f_invar, verbose)
  } else {
    vcf_arrow <- vcf_filter_columns(vcf_arrow, keep_samples, f_invar = FALSE, verbose)
  }

  return(vcf_arrow)
}
