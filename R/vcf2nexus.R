#' @title vcf2nexus
#'
#' @description
#' Converts a VCFArrow object to a NEXUS format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'nexus_infile.nex' (character)
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
#' vcf2nexus(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "nexus_infile.nex")
#' vcf2nexus(vcf_arrow, my_groups, out_file = "nexus_infile.nex")
#' vcf2nexus(vcf_arrow)
#'

vcf2nexus <- function(vcf_arrow, keep_groups = NULL,
                      out_file = "nexus_infile.nex") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Nexus")

  cli::cli_alert_info("Writing Nexus file...")

  write_nexus_cpp(acc$a1, acc$a2,
                  setup$variants$REF, setup$variants$ALT,
                  setup$samples, 0L, out_file)   # 0 = IUPAC DNA

  cli::cli_alert_success("Nexus file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
