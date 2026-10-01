#' @title vcf_filter_multiSNV
#'
#' @description
#' Subset a VCFArrow object keeping only loci with 2+ SNVs per locus
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param block_size -> size of linked SNV blocks, default 10000 bp (integer)
#' @param minSNV -> minimum linked block size, default 2 (integer)
#' @param maxSNV -> maximum number of selected linked SNVs per block, default 5 (integer)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object keeping only loci with between
#' min and max # of SNVs per locus, returning a new VCFArrow object.
#' Default min = 2 and max = 5 SNVs per locus
#' (recommended as input for fineRADstructure analyses).
#' Locus is defined as a different chromosome or a block of the
#' 'block_size' parameter value within a chromosome.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_multiSNV(vcf, block_size = 10000, minSNV = 2, maxSNV = 5)
#'
#' @export
#'

vcf_filter_multiSNV <- function(vcf_arrow, block_size = 10000,
                                minSNV = 2, maxSNV = 5) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  cli::cli_alert_info("Applying linked SNV filter")

  # select variants
  b <- .snv_blocks(vcf_arrow@variants, block_size)
  # count SNVs per block, and rank SNVs within each block (by POS)
  snvs_in_block <- tabulate(b$run)[b$run]
  rank <- seq_along(b$run) - which(b$first)[b$run] + 1L
  # keep only blocks with enough SNVs, and up to maxSNV SNVs per block
  keep <- b$row_id[snvs_in_block >= minSNV & rank <= maxSNV]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} (linked SNVs)"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
