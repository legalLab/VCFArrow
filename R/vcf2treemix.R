#' @title vcf2treemix
#'
#' @description
#' Converts a VCFArrow object to Treemix format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'treemix_infile.txt' (character)
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
#' vcf2treemix(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "treemix_infile.txt")
#' vcf2treemix(vcf_arrow, my_groups, out_file = "treemix_infile.txt")
#' vcf2treemix(vcf_arrow)
#'
#' @export
#'

vcf2treemix <- function(vcf_arrow, keep_groups = NULL,
                        out_file = "treemix_infile.txt") {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  write_treemix_header_cpp(setup$group_names, out_file)

  cli::cli_alert_info("Writing Treemix: {setup$n_var} variants x {setup$n_pops} pops")
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
