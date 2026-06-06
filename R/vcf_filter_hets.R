#' @title vcf_filter_hets
#'
#' @description
#' Remove loci above a heterozigosity threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> heterozigosity threshold, default 0.5 (numeric)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci above a heterozigosity threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' High heterozigosities are indicative of potential paralogs.
#'
#' @examples
#' vcf_filter_hets(vcf_arrow = my_vcf, threshold = 0.5)
#' vcf_filter_hets(my_vcf, 0.5)
#' vcf_filter_hets(my_vcf)
#'

vcf_filter_hets <- function(vcf_arrow, threshold = 0.5) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # compute heterozygosity per variant and select passing variants
  het_tbl <- vcf_arrow@gt |>
    dplyr::filter(!is.na(a1) & !is.na(a2)) |>
    dplyr::mutate(
      is_het = a1 != a2
    ) |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      n_het = sum(is_het),
      n_non_missing = dplyr::n(),
      het_rate = n_het / n_non_missing,
      .groups = "drop"
    ) |>
    dplyr::filter(het_rate < !!threshold) |>
    dplyr::select(.row_id)

  # collect only row_ids
  keep <- het_tbl |>
    dplyr::collect() |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
