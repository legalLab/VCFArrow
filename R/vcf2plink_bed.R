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
#
# Columns 1-4 of the .bim are built by .plink_variant_fields(), shared with
# vcf2plink_ped()'s .map writer.  That helper is also what enforces integer
# chromosome codes; writing CHROM verbatim is what makes every ADMIXTURE
# run abort with "Invalid chromosome code! Use integers."

vcf2plink_bed <- function(vcf_arrow, keep_groups = NULL,
                          out_file = "plink_out",
                          sex = NULL, pheno = NULL,
                          chrom_code = c("auto", "index", "zero", "keep")) {

  chrom_code <- match.arg(chrom_code)
  
  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  
  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  fields <- .plink_fam_fields(setup, sex, pheno)
  vfields <- .plink_variant_fields(setup, chrom_code)

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

  # ── .bim ─────────────────────────────────────────────────────────────────
  # Six fields, no header: CHROM, ID, genetic distance (dummy 0), POS,
  # allele1 (REF, clear bits), allele2 (ALT, set bits).
  utils::write.table(
    data.frame(CHROM = vfields$chrom,
               ID = vfields$id,
               cM = vfields$cm,
               POS = vfields$pos,
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
  
  # ── .chrommap ────────────────────────────────────────────────────────────
  .write_chrom_map(vfields, out_file)

  cli::cli_alert_success(
    "PLINK binary fileset written to \\
     {.file {bed_file}}, {.file {paste0(out_file, '.bim')}}, \\
     {.file {paste0(out_file, '.fam')}}"
  )

  invisible(vcf_arrow)
}
