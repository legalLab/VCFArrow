#' @title write_vcf
#'
#' @description
#' Write a VCFArrow object to an external VCF file
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of the VCF file to be written to, default 'output.vcf' (character)
#' @param gzip -> a flag to GZIP VCF when writing, default FALSE (Boolean)
#'
#' @return NULL
#'
#' @details
#' This function writes a VCFArrow object to an external VCF file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' It writes both uncompressed and gz compressed files. Compressing increases
#' writing time be about 50%.
#' For large files, it is recommended to output an uncompressed VCF file,
#' and then compress with GZIP or PIGZ.
#'
#' @examples
#' write_vcf(vcf = my_vcf, out_file = "output.vcf", gzip = FALSE)
#' write_vcf(my_vcf, "output.vcf", FALSE)
#' write_vcf(my_vcf)
#'

write_vcf <- function(vcf_arrow, out_file = "output.vcf", gzip = FALSE) {

  # write header
  header <- vcf_arrow@header
  header[length(header)] <- paste(
    c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
      vcf_arrow@samples),
    collapse = "\t"
  )
  writeLines(header, out_file)

  variants <- vcf_arrow@variants
  n_samples <- length(vcf_arrow@samples)
  samples <- vcf_arrow@samples

  # locate and sort feather files by chunk number
  # natural numeric sort avoids chunk_10 < chunk_2 with plain sort()
  feather_files <- list.files(vcf_arrow@path, pattern = "\\.arrow$",
                              full.names = TRUE)
  if (length(feather_files) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")

  chunk_nums <- as.integer(stringr::str_extract(basename(feather_files), "\\d+"))
  feather_files <- feather_files[order(chunk_nums)]

  # join FORMAT lookup (per variant) keyed by .row_id
  fmt_lookup <- vcf_arrow@format

  # chunk size message
  cli::cli_alert_info("VCF is being written in {length(feather_files)} chunks")

  # set up progress bar
  cli::cli_progress_bar("Writing VCF chunk", total = length(feather_files))

  for (fpath in feather_files) {

    chunk <- arrow::read_feather(fpath)  # columns: .row_id, sample, a1, a2, phased, fmt, ...

    # remove filtered samples from chunks
    # filtering here mirrors the .reshape_chunk() logic used in all exporters
    chunk <- chunk[chunk$sample %in% samples, , drop = FALSE]

    # sort so rows are: variant 1 sample 1, variant 1 sample 2, ..., variant 2 sample 1, ...
    # in practice read_vcf() writes them this way already, but sort defensively
    chunk <- chunk[order(chunk$.row_id, match(chunk$sample, samples)), ]

    row_ids <- unique(chunk$.row_id) # integer vector, one entry per variant in chunk
    n_chunk <- length(row_ids)

    # fmt_vec: flat row-major character vector, length n_chunk * n_samples
    # matrix() with byrow=TRUE then as.vector() gives column-major; use t() to stay row-major
    fmt_mat <- matrix(chunk$fmt, nrow = n_chunk, ncol = n_samples, byrow = TRUE)
    fmt_vec <- as.vector(t(fmt_mat)) # row-major: [var1_s1, var1_s2, ..., var2_s1, ...]

    # variant metadata for this chunk (indexed by .row_id)
    vi <- match(row_ids, variants$.row_id)

    # progress update
    cli::cli_progress_update()

    write_vcf_chunk_cpp(
      output_file = out_file,
      chrom = variants$CHROM[vi],
      pos = variants$POS[vi],
      id = variants$ID[vi],
      ref = variants$REF[vi],
      alt = variants$ALT[vi],
      qual = variants$QUAL[vi],
      filter_col = variants$FILTER[vi],
      info = vcf_arrow@info[vi],
      format_col = fmt_lookup$FORMAT[match(row_ids, fmt_lookup$.row_id)],
      fmt_vec = fmt_vec,
      n_samples = n_samples,
      gzip = gzip
    )
  }

  # end of progress bar
  cli::cli_progress_done()
  cli::cli_alert_success("VCFArrow object successfully written to {.file {out_file}}")

  invisible(out_file)
}
