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
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`
#'   and its sample file.
#'
#' @details
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The genotype file has one line per SNP and one column per sample,
#' separated by spaces: the number of ALT alleles (0, 1, 2), 9 if missing.
#' It has no header line, as smartsnp requires ("without row names or
#' column headings").
#' The samples and their groups, in the order of the columns, are written to
#' a separate tab-separated file with a header line (sample, group), named
#' after out_file without its extension plus '_samples.smartsnp' (e.g.
#' 'data.txt' gives 'data_samples.smartsnp').  Its group column is what
#' smartsnp's analyses take as sample_group, e.g.
#' `s <- read.delim("data_samples.smartsnp")` then
#' `smartsnp::smart_pca("data.txt", sample_group = s$group)` or
#' `smartsnp::smart_permanova("data.txt", sample_group = s$group)`.
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

  # samples and groups, in column order, in their own file: smartsnp reads
  # the genotype file without a header
  samples_file <- paste0(tools::file_path_sans_ext(out_file), "_samples.smartsnp")
  utils::write.table(
    data.frame(sample = setup$samples, group = setup$samples_groups,
               stringsAsFactors = FALSE),
    file = samples_file, sep = "\t", quote = FALSE, row.names = FALSE
  )

  # the chunks are appended to out_file, so it must exist and be empty
  file.create(out_file)

  cli::cli_alert_info("Building SmartSNP: {setup$n_var} variant{?s} x {setup$n_samples} sample{?s} \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)")
  cli::cli_alert_info("Writing SmartSNP file...")
  .write_chunks_ordered(.file_tasks(setup$feather_files), .smartsnp_chunk,
                        shared = .reshape_shared(setup), out_file = out_file,
                        label = "Writing chunk")
  cli::cli_alert_success(
    "SmartSNP files written to {.file {out_file}} and {.file {samples_file}}"
  )

  invisible(vcf_arrow)
}

# Worker side of the chunk loop (see .write_chunks_ordered()): append one chunk to task$part
.smartsnp_chunk <- function(task, shared) {
  rc <- .read_reshape_chunk(task$fpath, shared)
  if (!is.null(rc)) write_smartsnp_chunk_cpp(rc$a1, rc$a2, task$part)
  invisible(NULL)
}
