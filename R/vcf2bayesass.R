#' @title vcf2bayesass
#'
#' @description
#' Converts a VCFArrow object to BayesAss format infile
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
#' This function converts a VCFArrow object to an external BayesAss formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz", package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2bayesass(vcf, out_file = tempfile(fileext = ".immanc"))
#'
#' @export
#'

vcf2bayesass <- function(vcf_arrow, out_file, keep_groups = NULL) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

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
