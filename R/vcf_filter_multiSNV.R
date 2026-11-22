#' @title vcf_filter_multiSNV
#'
#' @description
#' Subset a VCFArrow object keeping only loci with several linked SNVs
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param block_size -> size of linked SNV blocks, default 10000 bp (integer)
#' @param minSNV -> minimum number of SNVs per block; blocks with fewer are
#'   removed, default 2 (integer)
#' @param maxSNV -> maximum number of SNVs per block, default 5 (integer); Inf
#'   keeps every SNV of each block
#' @param over_max -> what to do with blocks with more than maxSNV SNVs:
#'   "remove" the whole block or "truncate" it to its first maxSNV SNVs (by
#'   position), default "remove" (character)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object keeping only linked SNVs, i.e.
#' blocks (loci) with at least minSNV SNVs, returning a new VCFArrow object.
#' A block is a CHROM, or within a CHROM a window of block_size bp counted
#' from its first SNV; with de novo RAD data, where each locus is its own
#' CHROM, a block is a locus.
#' Blocks with fewer than minSNV SNVs are removed. Blocks with more than
#' maxSNV SNVs are removed (over_max = "remove") or truncated to their first
#' maxSNV SNVs (over_max = "truncate"); blocks with minSNV to maxSNV SNVs are
#' kept whole.
#' With RAD loci and de novo SNV discovery, an unusually high number of SNVs
#' in a short locus often indicates misassembled paralogs, so such loci are
#' best removed (the default; fineRADstructure's Stacks2fineRAD.py script
#' likewise removes loci with more than a maximum number of SNPs). With SNVs
#' from whole genome data mapped to a reference, a block is an arbitrary
#' window and many SNVs indicate a variable region rather than an artefact;
#' over_max = "truncate" keeps the first, most tightly linked, maxSNV SNVs of
#' such blocks, which suits unphased data, where the haplotypes of loci with
#' several heterozygous SNVs are arbitrary, and maxSNV = Inf keeps all SNVs,
#' e.g. for phased data.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_multiSNV(vcf, block_size = 10000, minSNV = 2, maxSNV = 5)
#' # whole genome data: truncate blocks with many SNVs instead of removing them
#' vcf_filter_multiSNV(vcf, minSNV = 2, maxSNV = 5, over_max = "truncate")
#'
#' @export
#'

vcf_filter_multiSNV <- function(vcf_arrow, block_size = 10000,
                                minSNV = 2, maxSNV = 5,
                                over_max = c("remove", "truncate")) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  over_max <- match.arg(over_max)
  if (!is.numeric(minSNV) || length(minSNV) != 1L || is.na(minSNV) || minSNV < 1)
    cli::cli_abort("{.arg minSNV} must be a number >= 1")
  if (!is.numeric(maxSNV) || length(maxSNV) != 1L || is.na(maxSNV) || maxSNV < minSNV)
    cli::cli_abort("{.arg maxSNV} must be a number >= {.arg minSNV} (or Inf)")

  idx <- .vcf_filter_index(vcf_arrow)

  cli::cli_alert_info("Applying linked SNV filter")

  # select variants
  b <- .snv_blocks(vcf_arrow@variants, block_size)
  # count SNVs per block, and rank SNVs within each block (by POS)
  snvs_in_block <- tabulate(b$run)[b$run]
  rank <- seq_along(b$run) - which(b$first)[b$run] + 1L
  # keep blocks with at least minSNV SNVs; blocks with more than maxSNV are
  # removed, or truncated to their first maxSNV SNVs
  keep <- if (over_max == "remove") {
    b$row_id[snvs_in_block >= minSNV & snvs_in_block <= maxSNV]
  } else {
    b$row_id[snvs_in_block >= minSNV & rank <= maxSNV]
  }
  n_over <- sum(b$first & snvs_in_block > maxSNV)

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} (linked SNVs)"
  )
  if (is.finite(maxSNV))
    cli::cli_alert_info(
      "{n_over} block{?s} with more than {maxSNV} SNVs \\
       {if (over_max == 'remove') 'removed' else 'truncated'}"
    )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
