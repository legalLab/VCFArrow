#' @title vcf_filter_rows
#'
#' @description
#' Unified API for filtering VCFArrow objects by row IDs.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep -> rows to keep (numeric or logical)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report removed samples and final % missing data, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes all rows in the 'keep' parameter from a VCFArrow object,
#' returning a new VCFArrow object.
#' Optionally will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' Optionally will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_filter_rank(vcf_arrow = my_vcf, keep = rows_to_keep)
#' vcf_filter_rank(my_vcf, rows_to_keep)
#'
#' @export
#'

vcf_filter_columns <- function(vcf_arrow, keep, f_invar = TRUE, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  s <- vcf_arrow@samples
  keep <- if (is.logical(keep)) {
    if (length(keep) != length(s)) cli::cli_abort("Logical filter length must match samples")
    s[keep]
  } else if (is.numeric(keep)) {
    s[keep]
  } else if (is.character(keep)) {
    intersect(s, keep)
  } else {
    cli::cli_abort("{.arg keep} must be character, logical, or numeric")
  }

  removed <- setdiff(s, keep)

  # Metadata update only — no @gt semi_join.
  idx <- match(keep, s)
  vcf_arrow@samples <- keep
  vcf_arrow@groups <- vcf_arrow@groups[idx]

  if (f_invar) vcf_arrow <- vcf_filter_invariant(vcf_arrow)

  if (verbose) {
    if (length(removed) > 0L)
      cli::cli_alert_info("Removed samples: {removed}")
    # Report current variant and sample counts from metadata — no @gt scan.
    cli::cli_alert_info(
      "Variants retained: {nrow(vcf_arrow@variants)} | \\
       Samples retained: {length(vcf_arrow@samples)}"
    )
  }

  return(vcf_arrow)
}
