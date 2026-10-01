#' @title vcf2smartsnp
#'
#' @description
#' Converts a VCFArrow object to smartsnp tabular format infile
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
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2smartsnp(vcf, out_file = tempfile(fileext = ".txt"))
#'
#' @export
#'

vcf2smartsnp <- function(vcf_arrow, out_file, keep_groups = NULL) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  write_smartsnp_header_cpp(setup$samples, out_file)

  cli::cli_alert_info("Building SmartSNP: {setup$n_var} variant{?s} x {setup$n_samples} sample{?s} \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)")
  cli::cli_alert_info("Writing SmartSNP file...")
  .write_chunks_ordered(.file_tasks(setup$feather_files), .smartsnp_chunk,
                        shared = .reshape_shared(setup), out_file = out_file,
                        label = "Writing chunk")
  cli::cli_alert_success("SmartSNP file written to {.file {out_file}}")

  invisible(vcf_arrow)
}

# Worker side of the chunk loop (see .write_chunks_ordered()): append one chunk to task$part
.smartsnp_chunk <- function(task, shared) {
  rc <- .read_reshape_chunk(task$fpath, shared)
  if (!is.null(rc)) write_smartsnp_chunk_cpp(rc$a1, rc$a2, task$part)
  invisible(NULL)
}
