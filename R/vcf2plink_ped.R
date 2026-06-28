#' @title vcf2fasta
#'
#' @description
#' Converts a VCFArrow object to a FASTA format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'fasta_infile.fas' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external NEXUS formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2fasta(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "fasta_infile.fas")
#' vcf2fasta(vcf_arrow, my_groups, out_file = "fasta_infile.fas")
#' vcf2fasta(vcf_arrow)
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

vcf2plink_ped <- function(vcf_arrow, keep_groups = NULL,
                          out_file = "plink_out",
                          sex = NULL, pheno = NULL) {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  fields <- .plink_fam_fields(setup, sex, pheno)
  acc <- .accumulate_individuals(setup, "PLINK .ped")

  cli::cli_alert_info("Writing PLINK file...")
  write_plink_ped_cpp(
    acc$a1, acc$a2,
    setup$variants$REF, setup$variants$ALT,
    setup$samples,
    fields$fid, fields$pat, fields$mat, fields$sex, fields$pheno,
    paste0(out_file, ".ped")
  )

  # ── genetic position ───────────────────────────────────────────────────────
  #
  # PLINKS's .bim column 3 is the genetic position in Morgans.
  # Without a recombination map, dividing physical position (bp) by 1 000 000
  # yields position in Mb, which is a standard proxy for Morgans and is
  # accepted by all major EIGENSTRAT-compatible tools (ADMIXTOOLS, smartpca).
  #
  # setup$variants$POS is the physical position for every retained, filtered
  # variant, already arranged in .row_id order by .vcf_export_setup().

  rel_pos_vec <- setup$variants$POS / 1e6

  # ── .map ─────────────────────────────────────────────────────────────────
  # Four fields, no header: CHROM, ID, genetic distance (dummy 0), POS.
  utils::write.table(
    data.frame(CHROM = setup$variants$CHROM,
               ID = setup$loci,
               cM = rel_pos_vec,
               POS = setup$variants$POS,
               stringsAsFactors = FALSE),
    file = paste0(out_file, ".map"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  cli::cli_alert_success(
    "PLINK text fileset written to \\
     {.file {paste0(out_file, '.ped')}}, {.file {paste0(out_file, '.map')}}"
  )

  invisible(vcf_arrow)
}
