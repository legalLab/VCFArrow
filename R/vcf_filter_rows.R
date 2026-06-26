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
  keep_mask <- v$.row_id %in% keep

  vcf_arrow@variants <- v[keep_mask, , drop = FALSE]

  # Keep @info positionally parallel to @variants going forward.
  if (length(vcf_arrow@info) == nrow(v)) {
    vcf_arrow@info <- vcf_arrow@info[keep_mask]
  } else {
    # @info was ALREADY misaligned coming into this call — almost certainly
    # an object that was filtered before this fix existed. Subsetting it
    # here with keep_mask would silently extract the WRONG entries (mask is
    # sized for the current @variants, not whatever stale length @info is
    # at). Leave @info untouched and warn loudly instead of corrupting it
    # further. See the one-time repair snippet for how to fix such objects.
    cli::cli_warn(c(
      "@info length ({length(vcf_arrow@info)}) does not match @variants \\
       length ({nrow(v)}) BEFORE filtering — @info is already misaligned \\
       and was left untouched rather than risk corrupting it further.",
      "i" = "This object was likely filtered before the @info alignment \\
             fix was applied. See the one-time repair snippet."
    ))
  }

  return(vcf_arrow)
}
