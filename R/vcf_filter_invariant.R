#' @title vcf_filter_invariant
#'
#' @description
#' Remove invariant loci from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes invariant loci from a VCFArrow object,
#' returning a new vcfR object.
#' This might be desirable after subsetting a VCFArrow object by individuals.
#'
#' @examples
#' vcf_filter_invariant(vcf_arrow = my_vcf)
#' vcf_filter_invariant(my_vcf)
#'

vcf_filter_invariant <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # compute number of unique genotypes per variant and select passing variants
  var_tbl <- vcf_arrow@gt |>
    dplyr::filter(!is.na(a1) & !is.na(a2)) |>
    dplyr::mutate(
      gt_str = paste0(a1, "/", a2)
    ) |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      n_unique = dplyr::n_distinct(gt_str),
      .groups = "drop"
    ) |>
    dplyr::filter(n_unique > 1) |>
    dplyr::select(.row_id)

  # collect row IDs only
  keep <- var_tbl |>
    dplyr::collect() |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
