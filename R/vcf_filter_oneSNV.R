#' @title vcf_filter_oneSNV
#'
#' @description
#' Subset a VCFArrow object keeping only 1 SNV per locus
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param block_size -> size of linked SNV blocks, default 10000 bp (integer)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object keeping only 1 SNV per locus,
#' returning a new VCFArrow object.
#' Locus is defined as a different chromosome or a block of the
#' 'block_size' parameter value within a chromosome.
#' The first SNV independent of quality is taken (may modify this in the future).
#'
#' @examples
#' vcf_filter_oneSNV(vcf_arrow = my_vcf, block_size = 10000)
#' vcf_filter_oneSNV(my_vcf, 10000)
#' vcf_filter_oneSNV(my_vcf)
#'
#' @export
#'

vcf_filter_oneSNV <- function(vcf_arrow, block_size = 10000) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  cli::cli_alert_info("Applying unlinked SNV filter")

  # select variants
  keep <- vcf_arrow@variants |>
    dplyr::arrange(CHROM, POS) |>
    dplyr::group_by(CHROM) |>
    dplyr::mutate(
      block = ((POS - min(POS)) %/% block_size) + 1
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct(CHROM, block, .keep_all = TRUE) |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variants (unlinked SNVs)"
  )

  return(vcf_arrow)
}
