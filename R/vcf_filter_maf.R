#' @title vcf_filter_maf
#'
#' @description
#' Remove loci below MAF threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal MAF threshold, default 0.05 (numeric)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci below MAF threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_maf(vcf_arrow = my_vcf, threshold = 0.05)
#' vcf_filter_maf(my_vcf, 0.05)
#' vcf_filter_maf(my_vcf)
#'

vcf_filter_maf <- function(vcf_arrow, threshold = 0.05) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # select passing variants
  maf_tbl <- vcf_arrow@gt |>
    # keep only called genotypes
    dplyr::filter(!is.na(a1) & !is.na(a2)) |>
    dplyr::mutate(
      allele_sum = a1 + a2
    ) |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      allele_sum = sum(allele_sum),
      n_called = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      af = allele_sum / (2 * n_called),
      maf = pmin(af, 1 - af)
    ) |>
    dplyr::filter(maf >= !!threshold) |>
    dplyr::select(.row_id)

  # collect only row_ids
  keep <- maf_tbl |>
    dplyr::collect() |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
