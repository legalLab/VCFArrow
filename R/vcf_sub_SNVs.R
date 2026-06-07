#' @title vcf_sub_loci
#'
#' @description
#' Randomly subsets SNVs from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param n_loci -> number of SNVs to subset, default 10000 (integer)
#' @param seed -> random number generator seed, default NULL (Boolean)
#'
#' @return VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object to specific number of SNVs,
#' returning new VCFArrow object.
#' The seed for random number generator is automatically generated unless
#' otherwise specified.
#'
#' @examples
#' vcf_sub_SNVs(vcf_arrow = my_vcf, n_SNVs = n_SNVs, seed = my_seed)
#' vcf_sub_SNVs(my_vcf, n_SNVs, 42)
#' vcf_sub_SNVs(my_vcf)
#'

vcf_sub_SNVs <- function(vcf_arrow, n_SNVs = 10000, seed = NULL) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  variants <- vcf_arrow@variants
  n_vars <- nrow(variants)

  if (n_SNVs >= n_vars) {
    cli::cli_alert_warning("Number of SNVs to subsample ({n_SNVs}) is greater than number of variants available ({n_vars});
                           returning original VCFArrow object")
    return(vcf_arrow)
    break
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  keep <- sort(sample(variants$.row_id, n_SNVs))

  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
