#' @title vcf2plink_ped
#'
#' @description
#' Converts a VCFArrow object to a PLINK .ped format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> output file prefix; extensions are added automatically, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param sex -> vector of sexes of samples, default NULL (character)
#' @param pheno -> vector of phenotypes of samples, default NULL (character)
#' @param chrom_code -> how CHROM is written: "auto" (keep numeric CHROM,
#'   otherwise 0), "index" (1..n by first appearance), "zero" (all 0) or
#'   "keep" (verbatim), default "auto" (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external PLINK .ped formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Sex and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz", package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2plink_ped(vcf, out_file = file.path(tempdir(), "plink_out"))
#'
#' @export
#'
#
# =============================================================================
# vcf2plink_ped()  – accumulate individuals, write .ped + .map
# =============================================================================
#
# Produces two files:
#   <out_file>.ped  – text pedigree + genotype table (individual-major)
#   <out_file>.map  – variant metadata (CHROM, ID, cM, POS)
#
# Alleles written as REF/ALT nucleotide letters; missing: "0".
# Memory: O(n_samples x n_var), same pattern as vcf2structure()/vcf2arlequin().
#
# The .map shares columns 1-4 with vcf2plink_bed()'s .bim and is built by the
# same .plink_variant_fields() helper, which is what enforces integer
# chromosome codes.  A contig name in column 1 does not produce a helpful
# message on this input path — ADMIXTURE reports only "PLINK Input file error".

vcf2plink_ped <- function(vcf_arrow, out_file, keep_groups = NULL,
                          sex = NULL, pheno = NULL,
                          chrom_code = c("auto", "index", "zero", "keep")) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  chrom_code <- match.arg(chrom_code)
  
  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  
  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  fields <- .plink_fam_fields(setup, sex, pheno)
  vfields <- .plink_variant_fields(setup, chrom_code)
  acc <- .accumulate_individuals(setup, "PLINK .ped")

  cli::cli_alert_info("Writing PLINK file...")
  write_plink_ped_cpp(
    acc$a1, acc$a2,
    setup$variants$REF, setup$variants$ALT,
    setup$samples,
    fields$fid, fields$pat, fields$mat, fields$sex, fields$pheno,
    paste0(out_file, ".ped")
  )

  # ── .map ─────────────────────────────────────────────────────────────────
  # Four fields, no header: CHROM, ID, genetic distance (dummy 0), POS.
  utils::write.table(
    data.frame(CHROM = vfields$chrom,
               ID = vfields$id,
               cM = vfields$cm,
               POS = vfields$pos,
               stringsAsFactors = FALSE),
    file = paste0(out_file, ".map"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )
  
  # ── .chrommap ────────────────────────────────────────────────────────────
  .write_chrom_map(vfields, out_file)

  cli::cli_alert_success(
    "PLINK text fileset written to \\
     {.file {paste0(out_file, '.ped')}}, {.file {paste0(out_file, '.map')}}"
  )

  invisible(vcf_arrow)
}
