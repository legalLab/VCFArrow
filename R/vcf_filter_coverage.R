#' @title vcf_filter_coverage
#'
#' @description
#' Remove genotypes below read DP threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> DP threshold, default 10 (integer)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes genotypes below a DP threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_coverage(vcf_arrow = my_vcf, threshold = 10)
#' vcf_filter_coverage(my_vcf, 10)
#' vcf_filter_coverage(my_vcf)
#'

vcf_filter_coverage <- function(vcf_arrow, threshold = 10) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # treat low coverage genotypes as missing and select passing variants
  var_tbl <- vcf_arrow@gt |>
    dplyr::filter(!is.na(a1) & !is.na(a2) & DP >= !!threshold) |>
    dplyr::mutate(
      gt_str = ifelse(a1 <= a2,
                      paste0(a1, "/", a2),
                      paste0(a2, "/", a1))
    ) |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      n_unique = dplyr::n_distinct(gt_str),
      .groups = "drop"
    ) |>
    dplyr::filter(n_unique > 1) |>
    dplyr::select(.row_id)

  # get row IDs
  keep <- var_tbl |>
    dplyr::collect() |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
