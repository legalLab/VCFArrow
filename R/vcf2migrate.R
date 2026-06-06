#' @title vcf2migrate
#'
#' @description
#' Converts a VCFArrow object to a MIGRATE-N format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'migrateN_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external MIGRATE-N formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The function implements two different SNP formats,
#' new (un)linked SNPs (S) and new (un)linked SNPs (N).
#' Alleles are sequence data when SNPs are mapped to a reference.
#' The S option generates blocks of a specific number of SNPs and treats them as linked;
#' this is appropriate if SNPs are extracted without a reference.
#' The N format extracts all SNPs within a chromosome, within the block size and treats them as linked;
#' this is appropriate if SNPs are mapped against a reference.
#' The size of the linked block is determined by the block_size parameter.
#' See https://peterbeerli.com/programs/migrate/distribution_4.x/migratedoc4.x.pdf for format detail.
#'
#' @examples
#' vcf2migrate(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "migrateN_infile.txt")
#' vcf2migrate(vcf_arrow, my_groups, out_file = "migrateN_infile.txt")
#' vcf2migrate(vcf_arrow)
#'

vcf2migrate <- function(vcf_arrow, keep_groups = NULL,
                        out_file = "migrateN_infile.txt") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_pops(setup, "Migrate-N")
  loci <- setup$loci

  cli::cli_alert_info("Writing Migrate-N file ({setup$n_var} variants)...")

  write_migrate_cpp(acc$ref, acc$alt, acc$nobs,
                    setup$group_names, loci, out_file)

  cli::cli_alert_success("Migrate-N file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
