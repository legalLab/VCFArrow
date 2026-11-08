#' @title vcf2structure
#'
#' @description
#' Converts a VCFArrow object to a Structure or FastStructure infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of file to output, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param method -> flag for Structure/FastStructure formats, default 'S' (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external Structure
#' or FastStructure formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The flag parameter controls whether Structure (flag = 'S') or
#' FastStructure (flag = 'F') formatted output is written out.
#' Each sample is written on two lines (one per allele copy): its name, its
#' group number and, for FastStructure, four columns of 0 (FastStructure
#' ignores the first six columns), then one allele per SNP: 0 (REF), 1 (ALT)
#' or -9 (missing).
#' FastStructure adds the extension '.str' to the name given with --input,
#' so give out_file the extension .str (e.g. "data.str", read with
#' --input=data --format=str).
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2structure(vcf, out_file = tempfile(fileext = ".str"))
#' vcf2structure(vcf, out_file = tempfile(fileext = ".str"),
#'               method = "F")
#'
#' @export
#'

vcf2structure <- function(vcf_arrow, out_file, keep_groups = NULL,
                          method = c("S", "F")) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  method <- match.arg(method)

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Structure")

  cli::cli_alert_info("Writing Structure file...")

  write_structure_cpp(acc$a1, acc$a2,
                      setup$samples, setup$group_ids,
                      if (method == "S") 0L else 1L, out_file)

  cli::cli_alert_success("Structure file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
