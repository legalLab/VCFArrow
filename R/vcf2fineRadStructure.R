#' @title vcf2fineradstructure
#'
#' @description
#' Converts a VCFArrow object to FineRadStructure format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of file to output, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param block_size -> size in bp of the blocks that define RAD loci, as in
#'   vcf_filter_multiSNV(), default 10000 (integer)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external fineRADstructure formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The first line holds the sample names; each further line is one RAD
#' locus, defined as in vcf_filter_multiSNV(): the SNVs on one CHROM within
#' the same block of block_size bases (counted from the CHROM's first SNV).
#' Each cell holds a sample's two haplotypes over the locus's SNVs in
#' RADpainter's syntax, e.g. 'TGAT/CGAT', or 'TGAT' if both are the same; a
#' sample missing at any SNV of the locus has an empty cell.
#' The genotypes are taken as phased: the first allele of each SNV (a1)
#' forms one haplotype and the second (a2) the other.  This holds for
#' phased VCFs, and for loci reconstructed as phased genotypes such as the
#' output of vcf_filter_multiSNV(); for unphased heterozygotes at several
#' SNVs of a locus the haplotypes are arbitrary.
#' Use the same block_size as in vcf_filter_multiSNV().
#' All genotypes are held in memory while writing (2 bytes per sample and
#' SNV, see vcf_memory_estimate(format = "individual")).
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_linked <- vcf_filter_multiSNV(vcf)
#' vcf2fineradstructure(vcf_linked, out_file = tempfile(fileext = ".txt"))
#'
#' @export
#'

vcf2fineradstructure <- function(vcf_arrow, out_file, keep_groups = NULL,
                                 block_size = 10000) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  # RAD loci: SNVs on one CHROM within the same block of block_size bases,
  # as in vcf_filter_multiSNV(); in CHROM / POS order
  b <- .snv_blocks(setup$variants, block_size)
  cols <- setup$var_pos[b$row_id]
  starts <- which(b$first)

  # a locus's SNVs can be in different chunk files, so all genotypes are
  # collected first
  acc <- .accumulate_individuals(setup, "fineRADstructure")

  cli::cli_alert_info(
    "Writing fineRADstructure file: {length(starts)} {cli::qty(length(starts))}loc{?us/i} \\
     from {length(cols)} SNV{?s}..."
  )
  write_fineradstructure_cpp(acc$a1, acc$a2, setup$variants$REF, setup$variants$ALT,
                             setup$samples, as.integer(cols), as.integer(starts), out_file)

  cli::cli_alert_success("fineRADstructure file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
