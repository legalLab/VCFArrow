#' @title vcf_filter_multiSNP
#'
#' @description
#' Subset a VCFArrow object keeping only loci with 2+ SNPs per locus
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param block_size -> size of linked SNP blocks, default 10000 bp (integer)
#' @param minSNP -> minimum linked block size, default 2 (integer)
#' @param maxSNP -> maximum number of selected linked SNPs per block, default 5 (integer)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object keeping only loci with between
#' min and max # of SNPs per locus, returning a new VCFArrow object.
#' Default min = 2 and max = 5 SNPs per locus
#' (recommended as input for fineRADstructure analyses).
#' Locus is defined as a different chromosome or a block of the
#' 'block_size' parameter value within a chromosome.
#'
#' @examples
#' vcf_filter_multiSNP(vcf_arrow = my_vcf, block_size = 10000, minSNP = 2, maxSNP = 5)
#' vcf_filter_multiSNP(my_vcf, 10000, 2, 5)
#' vcf_filter_multiSNP(my_vcf)
#'

vcf_filter_multiSNP <- function(vcf_arrow, block_size = 10000, minSNP = 2, maxSNP = 5) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # select variants
  keep <- vcf_arrow@variants |>
    dplyr::arrange(CHROM, POS) |>
    dplyr::group_by(CHROM) |>
    dplyr::mutate(
      first_pos = min(POS),
      block = ((POS - first_pos) %/% block_size) + 1
    ) |>
    dplyr::ungroup() |>
    # count SNPs per block
    dplyr::group_by(CHROM, block) |>
    dplyr::mutate(snps_in_block = dplyr::n()) |>
    # keep only blocks with enough SNPs
    dplyr::filter(snps_in_block >= minSNP) |>
    # rank SNPs within each block
    dplyr::arrange(CHROM, block, POS) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    # keep up to maxS SNPs per block
    dplyr::filter(rank <= maxSNP) |>
    dplyr::ungroup() |>
    dplyr::pull(.row_id)

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
