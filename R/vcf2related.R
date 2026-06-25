#' @title vcf2related
#'
#' @description
#' Converts a VCFArrow object to Related format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'related_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external Related formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2related(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "related_infile.txt")
#' vcf2related(vcf_arrow, my_groups, out_file = "related_infile.txt")
#' vcf2related(vcf_arrow)
#'
#' @export
#'

vcf2related <- function(vcf_arrow, keep_groups = NULL,
                        out_file = "related_infile.txt") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Related")

  cli::cli_alert_info("Writing Related file...")

  write_related_cpp(acc$a1, acc$a2,
                    setup$variants$REF, setup$variants$ALT,
                    setup$samples, out_file)

  cli::cli_alert_success("Related file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
