#' @title vcf2bayesass
#'
#' @description
#' Converts a VCFArrow object to BayesAss format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'bayesass3_infile.immanc' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external BayesAss formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2bayesass(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "bayesass3_infile.immanc")
#' vcf2bayesass(vcf_arrow, my_groups, out_file = "bayesass3_infile.immanc")
#' vcf2bayesass(vcf_arrow)
#'
#' @export
#'

vcf2bayesass <- function(vcf_arrow, keep_groups = NULL,
                         out_file = "bayesass3_infile.immanc") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "BayesAss")

  cli::cli_alert_info("Writing BayesAss file...")

  write_bayesass_cpp(acc$a1, acc$a2,
                     setup$variants$REF, setup$variants$ALT,
                     setup$samples, setup$group_names, setup$group_sizes,
                     setup$loci, out_file)

  cli::cli_alert_success("BayesAss file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
