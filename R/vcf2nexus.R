#' @title vcf2nexus
#'
#' @description
#' Converts a VCFArrow object to a NEXUS format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of file to output, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external NEXUS formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Only biallelic SNVs are exported (multiallelic sites and indels are left
#' out), and every sample must be assigned to a group (see set_vcf_groups()),
#' also when keep_groups is NULL.
#'
#' Each sample is one sequence of its genotypes at the SNVs: the REF or
#' ALT nucleotide for homozygotes, the IUPAC ambiguity code of the two
#' nucleotides for heterozygotes (e.g. R for A/G), and '?' if missing.
#' The matrix is declared as DNA data (FORMAT DATATYPE=DNA MISSING=?).
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2nexus(vcf, out_file = tempfile(fileext = ".nex"))
#'
#' @export
#'

vcf2nexus <- function(vcf_arrow, out_file, keep_groups = NULL) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Nexus")

  cli::cli_alert_info("Writing Nexus file...")

  write_nexus_cpp(acc$a1, acc$a2,
                  setup$variants$REF, setup$variants$ALT,
                  setup$samples, 0L, out_file)   # 0 = IUPAC DNA

  cli::cli_alert_success("Nexus file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
