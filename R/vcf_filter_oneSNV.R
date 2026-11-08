#' @title vcf_filter_oneSNV
#'
#' @description
#' Subset a VCFArrow object keeping only 1 SNV per locus
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param block_size -> size of linked SNV blocks, default 10000 bp (integer)
#' @param by_quality -> keep the SNV with the highest QUAL of each block
#'   instead of the first, default FALSE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object keeping only 1 SNV per locus,
#' returning a new VCFArrow object.
#' Locus is defined as a different chromosome or a block of the
#' 'block_size' parameter value within a chromosome.
#' By default the first SNV (lowest POS) of each locus is kept, whatever its
#' quality.  With by_quality = TRUE the SNV with the highest QUAL is kept;
#' SNVs without QUAL ('.') rank below all others, and ties (or a locus
#' without any QUAL) go to the first SNV.  VCFs without QUAL values, such as
#' those of DiscoSnpRad, therefore give the same result as the default, with
#' a warning.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_oneSNV(vcf)
#' # keep the SNV with the highest QUAL of each locus
#' vcf_filter_oneSNV(vcf, by_quality = TRUE)
#'
#' @export
#'

vcf_filter_oneSNV <- function(vcf_arrow, block_size = 10000, by_quality = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)

  cli::cli_alert_info("Applying unlinked SNV filter")

  b <- .snv_blocks(vcf_arrow@variants, block_size)
  if (by_quality) {
    # select variants: the highest QUAL of each block; NA QUAL last, and the
    # stable sort leaves ties in POS order, so they go to the first variant
    qual <- suppressWarnings(as.numeric(vcf_arrow@variants$QUAL))
    qual <- qual[match(b$row_id, vcf_arrow@variants$.row_id)]
    if (all(is.na(qual)))
      cli::cli_alert_warning("No QUAL values: keeping the first SNV of each block.")
    o <- order(b$run, qual, decreasing = c(FALSE, TRUE), method = "radix", na.last = TRUE)
    keep <- b$row_id[sort(o[!duplicated(b$run[o])])]
  } else {
    # select variants: the first (lowest POS) variant of each block
    keep <- b$row_id[b$first]
  }

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} (unlinked SNVs{if (by_quality) ', highest QUAL' else ''})"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
