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

vcf2snmf <- function(vcf_arrow, keep_groups = NULL,
                     out_file = "snmf_infile.geno") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  # write_snmf_cpp opens in append mode, so the file must exist and be empty
  # before the first chunk call
  file.create(out_file)

  cli::cli_alert_info("Writing sNMF .geno: {setup$n_var} variant{?s} x {setup$n_samples} sample{?s}")
  cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc)) {
      write_snmf_cpp(rc$a1, rc$a2, out_file)
    }
    cli::cli_progress_update()
  }

  cli::cli_progress_done()
  cli::cli_alert_success("sNMF .geno file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
