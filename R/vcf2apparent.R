#' @title vcf2apparent
#'
#' @description
#' Converts a VCFArrow object to an Apparent format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param key -> relationship type (All, Pa, Mo, Fa, Off), default All (character)
#' @param out_file -> name of file to output, default 'apparent_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Possible relationships defined by the parameter 'kee' are All, Pa, Mo, Fa, Off.
#'
#' @examples
#' vcf2apparent(vcf_arrow = my_vcf, keep_groups = my_groups, key = my_key, out_file = "apparent_infile.txt")
#' vcf2apparent(vcf_arrow, my_groups, my_key, out_file = "apparent_infile.txt")
#' vcf2apparent(vcf_arrow)
#'
#' @export
#'

vcf2apparent <- function(vcf_arrow, keep_groups = NULL,
                         key = "All",
                         out_file = "apparent_infile.txt") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Apparent")

  cli::cli_alert_info("Writing Apparent file...")

  write_apparent_cpp(acc$a1, acc$a2,
                     setup$variants$REF, setup$variants$ALT,
                     setup$samples, key, out_file)

  cli::cli_alert_success("Apparent file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
