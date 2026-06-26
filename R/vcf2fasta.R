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

vcf2fasta <- function(vcf_arrow, keep_groups = NULL,
                      out_file = "fasta_infile.fas") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "FASTA")

  cli::cli_alert_info("Writing FASTA file...")

  write_fasta_cpp(acc$a1, acc$a2,
                  setup$variants$REF, setup$variants$ALT,
                  setup$samples, out_file)

  cli::cli_alert_success("FASTA file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
