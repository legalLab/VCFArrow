#' @title vcf2plink_bed
#'
#' @description
#' Converts a VCFArrow object to a PLINK .bed format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'plink_out' (character)
#' @param sex -> vector of sexes of samples, default NULL (character)
#' @param pheno -> vector of phenotypes of samples, default NULL (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external PLINK .bed formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Sex and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
#'
#' @examples
#' vcf2plink_bed(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "plink_out", sex = sex, pheno = pheno)
#' vcf2plink_bed(vcf_arrow, my_groups, out_file = "plink_out")
#' vcf2plink_bed(vcf_arrow)
#'
#' @export
#'
#
# =============================================================================
# vcf2plink_bed()  – chunk-by-chunk .bed writer + .bim/.fam metadata
# =============================================================================
#
# Produces three files:
#   <out_file>.bed  – binary genotype table (SNP-major), written chunk-by-chunk
#   <out_file>.bim  – variant metadata (CHROM, ID, cM, POS, allele1, allele2)
#   <out_file>.fam  – sample metadata (FID, IID, PAT, MAT, SEX, PHENOTYPE)
#
# allele1 = REF (clear bits in .bed), allele2 = ALT (set bits) — see the
# mapping note in write_plink_bed_chunk_cpp()'s header comment.
# Memory: O(chunk), identical pattern to vcf2smartsnp().

vcf2plink_bed <- function(vcf_arrow, keep_groups = NULL,
                          out_file = "plink_out",
                          sex = NULL, pheno = NULL) {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  fields <- .plink_fam_fields(setup, sex, pheno)

  # ── .bed ─────────────────────────────────────────────────────────────────
  bed_file <- paste0(out_file, ".bed")
  write_plink_bed_header_cpp(bed_file)

  cli::cli_alert_info("Building PLINK: {setup$n_var} variant{?s} x {setup$n_samples} sample{?s} \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)")
  cli::cli_alert_info("Writing PLINK files...")
  cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc)) write_plink_bed_chunk_cpp(rc$a1, rc$a2, bed_file)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

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

  # ── .bim ─────────────────────────────────────────────────────────────────
  # Six fields, no header: CHROM, ID, genetic distance (dummy 0), POS,
  # allele1 (REF, clear bits), allele2 (ALT, set bits).
  utils::write.table(
    data.frame(CHROM = setup$variants$CHROM,
               ID = setup$loci,
               cM = rel_pos_vec,
               POS = setup$variants$POS,
               allele1 = setup$variants$REF,
               allele2 = setup$variants$ALT,
               stringsAsFactors = FALSE),
    file = paste0(out_file, ".bim"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  # ── .fam ─────────────────────────────────────────────────────────────────
  # Six fields, no header: FID, IID, PAT, MAT, SEX, PHENOTYPE.
  utils::write.table(
    data.frame(FID = fields$fid,
               IID = setup$samples,
               PAT = fields$pat,
               MAT = fields$mat,
               SEX = fields$sex,
               PHENO = fields$pheno,
               stringsAsFactors = FALSE),
    file = paste0(out_file, ".fam"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  cli::cli_alert_success(
    "PLINK binary fileset written to \\
     {.file {bed_file}}, {.file {paste0(out_file, '.bim')}}, \\
     {.file {paste0(out_file, '.fam')}}"
  )

  invisible(vcf_arrow)
}
