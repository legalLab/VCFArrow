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
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Possible relationships defined by the parameter 'kee' are All, Pa, Mo, Fa, Off.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz", package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2apparent(vcf, out_file = tempfile(fileext = ".txt"))
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
