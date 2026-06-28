#' @title vcf2fineradstructure
#'
#' @description
#' Converts a VCFArrow object to FineRadStructure format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'fineradstructure_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external Treemix formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2fineradstructure(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "fineradstructure_infile.txt")
#' vcf2fineradstructure(vcf_arrow, my_groups, out_file = "fineradstructure_infile.txt")
#' vcf2fineradstructure(vcf_arrow)
#'
#' @export
#'

vcf2fineradstructure <- function(vcf_arrow, keep_groups = NULL,
                                 out_file = "fineradstructure_infile.txt") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  write_fineradstructure_header_cpp(setup$samples, out_file)

  cli::cli_alert_info("Building fineRADstructure: {setup$n_var} variants x {setup$n_samples} samples \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)")
  cli::cli_alert_info("Writing fineRADstructure file...")
  cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id","sample","a1","a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc))
      write_fineradstructure_chunk_cpp(rc$a1, rc$a2,
                                       setup$variants$REF[rc$col_idx],
                                       setup$variants$ALT[rc$col_idx],
                                       out_file)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  cli::cli_alert_success("fineRADstructure file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
