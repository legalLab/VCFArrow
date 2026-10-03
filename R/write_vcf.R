#' @title write_vcf
#'
#' @description
#' Write a VCFArrow object to an external VCF file
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of the VCF file to be written to, no default (character)
#' @param gzip -> a flag to GZIP VCF when writing, default FALSE (Boolean)
#'
#' @return Invisibly returns the path of the written file.
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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' write_vcf(vcf, out_file = tempfile(fileext = ".vcf"))
#' write_vcf(vcf, out_file = tempfile(fileext = ".vcf.gz"),
#'           gzip = TRUE)
#'
#' @export
#'

write_vcf <- function(vcf_arrow, out_file, gzip = FALSE) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  # write header (with gzip, as its own gzip member: the chunks are appended
  # as further members, which together form a valid multi-member gzip file)
  header <- vcf_arrow@header
  header[length(header)] <- paste(
    c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
      vcf_arrow@samples),
    collapse = "\t"
  )
  con <- if (gzip) gzfile(out_file, "w") else file(out_file, "w")
  writeLines(header, con)
  close(con)

  variants <- vcf_arrow@variants

  # locate and sort feather files by chunk number
  # natural numeric sort avoids chunk_10 < chunk_2 with plain sort()
  feather_files <- .get_sorted_feather_files(vcf_arrow@path)

  # Each chunk's task carries the metadata of the variants in its .row_id
  # range (CHROM ... INFO, and FORMAT from the per-variant lookup), so the
  # chunks can be written in parallel without sending all metadata to every
  # worker.
  fmt_lookup <- vcf_arrow@format
  meta <- data.frame(
    .row_id = variants$.row_id,
    CHROM = variants$CHROM, POS = variants$POS, ID = variants$ID,
    REF = variants$REF, ALT = variants$ALT, QUAL = variants$QUAL,
    FILTER = variants$FILTER, INFO = vcf_arrow@info,
    FORMAT = fmt_lookup$FORMAT[.match_row_id(variants$.row_id,
                                             .row_id_pos(fmt_lookup$.row_id))],
    stringsAsFactors = FALSE
  )
  tasks <- lapply(feather_files, function(f) {
    id <- as.vector(arrow::read_feather(f, col_select = ".row_id",
                                        as_data_frame = FALSE)$.row_id)
    in_range <- if (length(id)) meta$.row_id >= min(id) & meta$.row_id <= max(id)
                else logical(nrow(meta))
    list(fpath = f, meta = meta[in_range, , drop = FALSE])
  })
  meta <- NULL

  # chunk size message
  cli::cli_alert_info("VCF is being written in {length(feather_files)} chunk{?s}")

  .write_chunks_ordered(tasks, .write_vcf_chunk,
                        shared = list(samples = vcf_arrow@samples, gzip = gzip),
                        out_file = out_file, label = "Writing VCF chunk")

  cli::cli_alert_success("VCFArrow object successfully written to {.file {out_file}}")

  invisible(out_file)
}

# Worker side of write_vcf() (see .write_chunks_ordered()): append the VCF
# lines of one chunk to task$part.  task$meta: metadata of the variants in the
# chunk's .row_id range.

.write_vcf_chunk <- function(task, shared) {
  samples <- shared$samples
  meta <- task$meta
  var_pos <- .row_id_pos(meta$.row_id)

  # only the per-sample FORMAT strings are written; the other genotype
  # columns are derived from them at read time
  chunk <- arrow::read_feather(task$fpath, col_select = c(".row_id", "sample", "fmt"),
                               as_data_frame = FALSE)

  # remove filtered variants and samples from chunks
  # filtering here mirrors the .reshape_chunk() logic used in all exporters
  sample_order <- .sample_index(chunk, samples)
  row_id <- as.vector(chunk$.row_id)
  keep <- !is.na(.match_row_id(row_id, var_pos)) & !is.na(sample_order)
  # nothing to write if no variant of this chunk is retained (e.g. all were
  # filtered out); this also avoids integer(0)[TRUE], which is NA, below
  if (!any(keep)) return(invisible(NULL))

  # sort so rows are: variant 1 sample 1, variant 1 sample 2, ..., variant 2 sample 1, ...
  # in practice read_vcf() writes them this way already, but sort defensively
  # (rows: positions in the chunk of the kept rows, in that order)
  rows <- which(keep)
  rows <- rows[order(row_id[rows], sample_order[rows])]
  row_id <- row_id[rows]

  # one entry per variant in chunk (row_id is sorted)
  row_ids <- row_id[c(TRUE, row_id[-1L] != row_id[-length(row_id)])]

  # per-sample FORMAT strings in that order (flat, row-major:
  # [var1_s1, var1_s2, ..., var2_s1, ...]), selected in Arrow and handed to
  # the C++ writer as an Arrow array, so they never become R strings
  fmt <- chunk$fmt$Take(arrow::Array$create(rows - 1L))
  fmt <- if (fmt$num_chunks == 1L) fmt$chunk(0L)
         else do.call(arrow::concat_arrays, fmt$chunks)
  fmt_c <- arrow_c_alloc_cpp()
  fmt$export_to_c(fmt_c$array, fmt_c$schema)
  fmt <- chunk <- NULL

  # variant metadata for this chunk (indexed by .row_id)
  vi <- .match_row_id(row_ids, var_pos)

  write_vcf_chunk_cpp(
    output_file = task$part,
    chrom = meta$CHROM[vi],
    pos = meta$POS[vi],
    id = meta$ID[vi],
    ref = meta$REF[vi],
    alt = meta$ALT[vi],
    qual = meta$QUAL[vi],
    filter_col = meta$FILTER[vi],
    info = meta$INFO[vi],
    format_col = meta$FORMAT[vi],
    fmt_array = fmt_c$array,
    fmt_schema = fmt_c$schema,
    n_samples = length(samples),
    gzip = shared$gzip
  )
  invisible(NULL)
}
