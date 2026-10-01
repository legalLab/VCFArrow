#' @title read_vcf
#'
#' @description
#' Read VCF file and store its content within a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_file -> VCF file
#' @param chunk_size -> number of variants to read in at a time, default 50000 (integer)
#'
#' @return VCFArrow object
#'
#' @details
#' This function read a VCF file into an S4 class object in chunks,
#' returning a VCFArrow object.
#' It accepts both uncompressed and gz compressed files.
#' The GT field is stored as an Apache Arrow in Long format.
#' Various metrics are precalculated for fast and easy filtering.
#' GT slot content stored in a TEMP directory for lazy loading.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f)
#' vcf
#'
#' @export
#'

read_vcf <- function(vcf_file, chunk_size = 50000) {

  # helper function
  parse_format_header <- function(header) {

    fmt_lines <- grep("^##FORMAT=", header, value = TRUE)

    fmt_ids <- sub(".*ID=([^,]+).*", "\\1", fmt_lines)
    fmt_types <- sub(".*Type=([^,]+).*", "\\1", fmt_lines)

    return(data.frame(ID = fmt_ids, Type = fmt_types, stringsAsFactors = FALSE))
  }

  # generate a unique path name within the session tmp folder
  tmp_dir <- tempfile(pattern = "arrow_vcf_")

  # physically create the directory at that path
  dir.create(tmp_dir)

  # open the VCF (plain or gzip-compressed) with the C++ chunk reader
  con <- vcf_open_cpp(normalizePath(vcf_file, mustWork = TRUE))
  on.exit(vcf_close_cpp(con))

  # read header
  hdr <- vcf_read_header_cpp(con)
  header <- hdr$header
  if (hdr$status == 1L) cli::cli_abort("Unexpected EOF before header")
  if (hdr$status == 2L) cli::cli_abort("Malformed VCF: missing #CHROM line")

  # extract samples
  header_fields <- strsplit(header[length(header)], "\t", fixed = TRUE)[[1]]
  samples <- header_fields[-(1:9)]

  # placeholder for group information
  groups <- rep.int(NA_character_, length(samples))

  # storage
  variant_buffer <- list()
  info_buffer <- list()
  format_buffer <- list()
  chunk_id <- 1

  # R metadata that write_feather() attaches to a data.frame, so chunks
  # written from the C++-built Arrow batches read back exactly as before
  gt_metadata <- arrow::arrow_table(data.frame(
    .row_id = integer(0), sample = character(0), a1 = integer(0), a2 = integer(0),
    phased = logical(0), fmt = character(0), DP = numeric(0), GQ = numeric(0),
    ADR = numeric(0), stringsAsFactors = FALSE
  ))$metadata

  # chunk size message
  cli::cli_alert_info("VCF is being read in chunks of {chunk_size} variant{?s}")

  # set up progress bar
  cli::cli_progress_bar("Reading in VCF chunk", total = NA)

  # Chunks are parsed into Arrow record batches and written to feather files;
  # the per-variant fields are collected here.  With workers (see
  # vcf_set_workers()), this process reads raw chunks of lines and the
  # workers parse and write them, one chunk per worker per round.
  n_workers <- .vcf_workers()
  shared <- list(samples = samples, metadata = gt_metadata)
  chunk_path <- function(id) file.path(tmp_dir, paste0("chunk_", id, ".arrow"))

  repeat {
    if (n_workers == 1L) {
      parsed <- vcf_read_chunk_cpp(con, chunk_size, samples,
                                   as.integer((chunk_id - 1) * chunk_size))
      if (parsed$n == 0L) break
      done <- list(.write_gt_batch(parsed, chunk_path(chunk_id), gt_metadata))
    } else {
      tasks <- list()
      for (w in seq_len(n_workers)) {
        rc <- vcf_read_raw_cpp(con, chunk_size)
        if (rc$n == 0L) break
        id <- chunk_id + w - 1
        tasks[[w]] <- list(raw = rc$raw, row_offset = as.integer((id - 1) * chunk_size),
                           out = chunk_path(id))
      }
      if (length(tasks) == 0L) break
      done <- .map_chunks(tasks, .read_vcf_chunk, shared = shared, progress = FALSE)
    }

    for (parsed in done) {
      # progress update
      cli::cli_progress_update()

      # FORMAT is per-variant → for memory efficiency keep in a separate FORMAT lookup
      format_df <- data.frame(
        FORMAT = parsed$format,
        .row_id = as.integer((chunk_id - 1) * chunk_size) + seq_len(parsed$n),
        stringsAsFactors = FALSE
      )
      format_buffer[[chunk_id]] <- format_df

      variant_buffer[[chunk_id]] <- as.matrix(parsed$variants)
      info_buffer[[chunk_id]] <- as.character(parsed$info)

      chunk_id <- chunk_id + 1
    }
  }

  # end of progress bar
  cli::cli_progress_done()

  # combine FORMAT field
  format_df_all <- do.call(rbind, format_buffer)

  # combine metadata
  variants_df <- as.data.frame(do.call(rbind, variant_buffer), stringsAsFactors = FALSE)
  colnames(variants_df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER")
  variants_df$POS <- suppressWarnings(as.integer(variants_df$POS))
  info_vec <- unlist(info_buffer, use.names = FALSE)

  # DiscoSNP-RAD paralog metrics from INFO: RANK (Rk) and REPEAT (RPT);
  # NA where absent
  variants_df$Rk <- .info_numeric(info_vec, "Rk")
  variants_df$RPT <- .info_numeric(info_vec, "RPT")

  # extract variant information for filtering
  variants_df <- variants_df |>
    dplyr::mutate(
      n_alt = stringr::str_count(ALT, ",") + 1,
      is_biallelic = n_alt == 1,
      # indel: REF or any ALT allele not exactly one character long
      # (an allele of 2+ characters, or an empty allele before a comma)
      is_indel = (nchar(REF) != 1 |
                    grepl("[^,]{2,}", ALT) | grepl("(^|,),", ALT))
    )

  variants_df$.row_id <- as.integer(seq_len(nrow(variants_df)))

  # open Arrow Dataset (lazy, unified view)
  gt_arrow <- arrow::open_dataset(tmp_dir, format = "feather")
  # option to store as Apache parquet
  #gt_arrow <- arrow::open_dataset(tmp_dir, format = "parquet")

  # CHROM/POS/REF/ALT of every .row_id, kept unfiltered (row i = .row_id i);
  # shares its columns with variants_df until variants are filtered
  loci <- variants_df[c("CHROM", "POS", "REF", "ALT")]
  rownames(loci) <- NULL

  vcf_arrow <- .new_vcfarrow(
    header,
    info_vec,
    format_df_all,
    variants_df,
    gt_arrow,
    samples,
    groups,
    tmp_dir,
    loci = loci
  )

  cli::cli_alert_success("VCF successfully read into a VCFArrow object")

  return(vcf_arrow)
}

# Write the genotype record batch of a parsed chunk to a feather file and
# return the per-variant fields.  gt_metadata: the R metadata write_feather()
# attaches to a data.frame, so chunks read back exactly as before.

.write_gt_batch <- function(parsed, out, gt_metadata) {
  gt_long <- arrow::Table$create(
    arrow::RecordBatch$import_from_c(parsed$array, parsed$schema)
  )
  gt_long$metadata <- gt_metadata
  arrow::write_feather(gt_long, out)
  # option to store as Apache parquet
  #arrow::write_parquet(gt_long, sub("\\.arrow$", ".parquet", out))
  gt_long <- NULL
  .release_arrow_memory()
  parsed[c("n", "variants", "info", "format")]
}

# Worker side of parallel reading: parse one raw chunk (see
# vcf_read_raw_cpp()) and write its genotypes.

.read_vcf_chunk <- function(task, shared) {
  parsed <- vcf_parse_raw_cpp(task$raw, shared$samples, task$row_offset)
  .write_gt_batch(parsed, task$out, shared$metadata)
}

# Numeric value of INFO key `key` for each variant (NA where the key is absent
# or its value is not a number).  The key is matched only at the start of an
# INFO entry, so e.g. "RPT" does not match inside "MRPT=3".

.info_numeric <- function(info, key) {
  value <- stringr::str_match(info, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]
  suppressWarnings(as.numeric(value))
}
