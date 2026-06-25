#' @title vcf_filter_rows
#'
#' @description
#' Unified API for filtering VCFArrow objects by row IDs.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep -> rows to keep (numeric or logical)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes all rows in the 'keep' parameter from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_rank(vcf_arrow = my_vcf, keep = rows_to_keep)
#' vcf_filter_rank(my_vcf, rows_to_keep)
#'
#' @export
#'

.vcf_filter_rows <- function(vcf_arrow, keep) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  v <- vcf_arrow@variants
  keep <- if (is.logical(keep)) {
    if (length(keep) != nrow(v))
      cli::cli_abort("Logical filter length must match number of variants")
    v$.row_id[keep]
  } else if (is.numeric(keep)) {
    as.integer(keep)
  } else {
    cli::cli_abort("{.arg keep} must be logical or numeric (.row_id values)")
  }

  # Update variant metadata only.  @gt is NOT mutated — it stays as the
  # original open_dataset().  All read-time consumers filter via @variants$.row_id.
  vcf_arrow@variants <- v[v$.row_id %in% keep, , drop = FALSE]

  return(vcf_arrow)
}
