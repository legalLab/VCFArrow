#' @title vcf2arlequin
#'
#' @description
#' Converts a VCFArrow object to an Arlequin infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'arlequin_infile.arp' (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external Arlequin formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz", package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2arlequin(vcf, out_file = tempfile(fileext = ".arp"))
#'
#' @export
#'

vcf2arlequin <- function(vcf_arrow, keep_groups = NULL,
                         out_file = "arlequin_infile.arp") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc   <- .accumulate_individuals(setup, "Arlequin")

  cli::cli_alert_info("Writing Arlequin file...")

  write_arlequin_cpp(acc$a1, acc$a2,
                     setup$variants$REF, setup$variants$ALT,
                     setup$samples, setup$group_names, setup$group_sizes,
                     out_file)

  cli::cli_alert_success("Arlequin file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
