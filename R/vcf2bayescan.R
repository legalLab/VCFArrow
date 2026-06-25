#' @title vcf2bayescan
#'
#' @description
#' Converts a VCFArrow object to Bayescan format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'bayescan_infile.txt' (character)
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
#' vcf2bayescan(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "bayescan_infile.txt")
#' vcf2bayescan(vcf_arrow, my_groups, out_file = "bayescan_infile.txt")
#' vcf2bayescan(vcf_arrow)
#'

vcf2bayescan <- function(vcf_arrow, keep_groups = NULL,
                         out_file = "bayescan_infile.txt") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_pops_lowmem(setup, "BayesScan")

  cli::cli_alert_info("Writing BayesScan file...")

  write_bayescan_cpp(acc$ref, acc$alt, acc$nobs,
                     setup$group_names, out_file)

  cli::cli_alert_success("BayesScan file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
