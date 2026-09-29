#' @title vcf2treemix
#'
#' @description
#' Converts a VCFArrow object to Treemix format infile
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
#' This function converts a VCFArrow object to an external Treemix formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2treemix(vcf, out_file = tempfile(fileext = ".txt"))
#'
#' @export
#'

vcf2treemix <- function(vcf_arrow, out_file, keep_groups = NULL) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  write_treemix_header_cpp(setup$group_names, out_file)

  cli::cli_alert_info("Building Treemix: {setup$n_var} variant{?s} x {setup$n_pops} pop{?s} \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)")
  cli::cli_alert_info("Writing Treemix file...")
  cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id","sample","a1","a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc)) {
      pc <- .pop_counts_from_chunk(rc, setup)
      write_treemix_chunk_cpp(pc$ref, pc$alt, out_file)
    }

    cli::cli_progress_update()
  }

  cli::cli_progress_done()
  cli::cli_alert_success("Treemix file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
