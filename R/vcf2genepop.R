#' @title vcf2genepop
#'
#' @description
#' Converts a VCFArrow object to a Genepop format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of file to output, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param pop_names -> write the group name after Pop on the line that starts
#'   each population ('pop WA') instead of Pop alone, default FALSE (Boolean)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external Genepop formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Each population starts with a line holding the word Pop alone, as the
#' GENEPOP documentation shows it; GENEPOP itself names each population
#' after its last individual.  Some readers, e.g. adegenet::read.genepop(),
#' accept only this form.  With pop_names = TRUE, the line also holds the
#' group name ('pop WA'), which other programs may read as the population
#' name.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2genepop(vcf, out_file = tempfile(fileext = ".gen"))
#'
#' @export
#'

vcf2genepop <- function(vcf_arrow, out_file, keep_groups = NULL,
                        pop_names = FALSE) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Genepop")

  cli::cli_alert_info("Writing Genepop file...")

  write_genepop_cpp(acc$a1, acc$a2,
                    setup$variants$REF, setup$variants$ALT,
                    setup$samples, setup$group_names, setup$group_sizes,
                    setup$loci, isTRUE(pop_names), out_file)

  cli::cli_alert_success("Genepop file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
