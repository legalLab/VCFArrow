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
#' @param gzip -> compress the file with gzip, as TreeMix expects, default TRUE
#'   if out_file ends in ".gz" (Boolean)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external Treemix formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The first line holds the group names; each further line one SNP, with the
#' REF and ALT allele counts of each group as "ref,alt", separated by spaces.
#' TreeMix reads gzip-compressed input (treemix -i file.gz): name out_file
#' with the extension .gz to compress it.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2treemix(vcf, out_file = tempfile(fileext = ".gz"))
#'
#' @export
#'

vcf2treemix <- function(vcf_arrow, out_file, keep_groups = NULL,
                        gzip = grepl("\\.gz$", out_file)) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  # with gzip, the text is written to a temporary file and then compressed
  # into out_file
  final_file <- out_file
  if (gzip) {
    out_file <- tempfile("vcfarrow_treemix_")
    on.exit(unlink(out_file), add = TRUE)
  }

  write_treemix_header_cpp(setup$group_names, out_file)

  cli::cli_alert_info("Building Treemix: {setup$n_var} variant{?s} x {setup$n_pops} pop{?s} \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)")
  cli::cli_alert_info("Writing Treemix file...")
  .write_chunks_ordered(.file_tasks(setup$feather_files), .treemix_chunk,
                        shared = c(.reshape_shared(setup), list(P = setup$P)),
                        out_file = out_file, label = "Writing chunk")
  if (gzip) .gzip_file(out_file, final_file)
  cli::cli_alert_success("Treemix file written to {.file {final_file}}")

  invisible(vcf_arrow)
}

# Worker side of the chunk loop (see .write_chunks_ordered()): append one
# chunk's per-population allele counts to task$part
.treemix_chunk <- function(task, shared) {
  rc <- .read_reshape_chunk(task$fpath, shared)
  if (!is.null(rc)) {
    pc <- .pop_counts_from_chunk(rc, shared)
    write_treemix_chunk_cpp(pc$ref, pc$alt, task$part)
  }
  invisible(NULL)
}

# Compress file `from` into gzip file `to`, in blocks (so memory stays small)
.gzip_file <- function(from, to) {
  con_in <- file(from, "rb")
  on.exit(close(con_in), add = TRUE)
  con_out <- gzfile(to, "wb")
  on.exit(close(con_out), add = TRUE)
  while (length(b <- readBin(con_in, "raw", 16L * 1024L^2)) > 0L) writeBin(b, con_out)
  invisible(to)
}
