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

vcf_filter_rows <- function(vcf_arrow, keep) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  v <- vcf_arrow@variants

  # normalize input → row_ids
  if (is.logical(keep)) {
    if (length(keep) != nrow(v)) {
      cli::cli_abort("Logical filter must match number of variants")
    }
    keep_ids <- v$.row_id[keep]
  } else if (is.numeric(keep)) {
    keep_ids <- keep
  } else {
    cli::cli_abort("The 'keep' parameter must be logical or numeric (.row_id)")
  }

  # enforce integer type
  keep_ids <- as.integer(keep_ids)

  # filter variants
  vcf_arrow@variants <- v[
    v$.row_id %in% keep_ids,
    , drop = FALSE
  ]

  # Arrow filter
  keep_tbl <- arrow::arrow_table(.row_id = keep_ids)

  vcf_arrow@gt <- vcf_arrow@gt |>
    dplyr::semi_join(keep_tbl, by = ".row_id")

  return(vcf_arrow)
}
