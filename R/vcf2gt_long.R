#' @title vcf2gt_long
#'
#' @description
#' Converts a VCFArrow object to tidy long format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of file to output, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param format -> one of three output formats (arrow, parquet, CSV) (character)
#' @param col_select -> optional selection of columns to save, default ALL
#'
#' @return Invisibly returns the path of the written file.
#'
#' @details
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The tidy data can be saved in either Arrow, Parquet or CSV formats.
#' File extension is added automatically if missing.
#' Optionally, specific columns can be saved, by default all columns are saved.
#' The gt long slot contains pre-calculated metrics in addition to just genotypes.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2gt_long(vcf, out_file = tempfile(), format = "csv")
#'
#' @export
#'

vcf2gt_long <- function(vcf_arrow, out_file, keep_groups = NULL,
                        format = c("feather", "parquet", "csv"),
                        col_select = NULL) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  format <- match.arg(format)
  setup <-.vcf_export_setup(vcf_arrow, keep_groups)

  # Resolve output path (add extension if missing)
  ext_map <- c(feather = ".feather", parquet = ".parquet", csv = ".csv")
  ext <- ext_map[[format]]
  if (!endsWith(out_file, ext)) {
    out_file <- paste0(out_file, ext)
  }

  cli::cli_alert_info(
    "Writing gt_long: {setup$n_var} variant{?s} x {setup$n_samples} sample{?s}"
  )

  # Chunks are filtered, sorted and appended to out_file one at a time, so
  # memory stays O(chunk).  Sorting each chunk (variant-major, then sample
  # group order) sorts the whole table because feather chunks hold
  # consecutive, non-overlapping .row_id ranges (checked below).
  # With workers (see vcf_set_workers()) and CSV output, the chunks are
  # filtered, sorted and formatted as CSV text in parallel, in temporary part
  # files that this process appends to out_file in order.  Feather and
  # Parquet output stay serial: the rows must pass through this process's
  # single writer anyway, and writing and re-reading part files cost more
  # than the parallel sort saves (measured: 1.6x slower at 4-8 workers).
  read_cols <- if (is.null(col_select)) NULL else union(c(".row_id", "sample"), col_select)
  parallel <- format == "csv" && .vcf_workers() > 1L && length(setup$feather_files) > 1L
  part_dir <- NULL
  if (parallel) {
    part_dir <- tempfile("vcfarrow_parts_")
    dir.create(part_dir)
  }
  sink <- NULL
  writer <- NULL
  last_id <- 0L
  on.exit({
    if (!is.null(writer) && format != "csv") {
      if (format == "feather") writer$close() else writer$Close()
    }
    if (!is.null(sink)) sink$close()
    if (!is.null(part_dir)) unlink(part_dir, recursive = TRUE)
  }, add = TRUE)

  write_part <- function(tbl) {
    if (format == "csv") {
      df <- as.data.frame(tbl)
      if (is.null(writer)) {
        utils::write.csv(df, out_file, row.names = FALSE)
        writer <<- "csv"
      } else {
        # same settings as write.csv(), appended without a header
        utils::write.table(df, out_file, append = TRUE, sep = ",", dec = ".",
                           qmethod = "double", row.names = FALSE, col.names = FALSE)
      }
      return(invisible(NULL))
    }
    if (is.null(writer)) {
      sink <<- arrow::FileOutputStream$create(out_file)
      writer <<- if (format == "feather") {
        arrow::RecordBatchFileWriter$create(sink, tbl$schema)
      } else {
        arrow::ParquetFileWriter$create(
          tbl$schema, sink,
          properties = arrow::ParquetWriterProperties$create(names(tbl))
        )
      }
    }
    if (format == "feather") writer$write_table(tbl) else writer$WriteTable(tbl, chunk_size = tbl$num_rows)
  }

  # columns of the output, for the header / an empty result
  empty <- arrow::read_feather(setup$feather_files[[1L]], col_select = read_cols,
                               as_data_frame = FALSE)$Slice(0, 0)
  if (!is.null(col_select)) empty <- empty[, col_select]
  if (parallel) write_part(empty)  # CSV header

  tasks <- lapply(seq_along(setup$feather_files), function(i) list(
    fpath = setup$feather_files[[i]],
    part = if (parallel) file.path(part_dir, paste0("part_", i)) else NULL
  ))
  shared <- list(samples = setup$samples, var_pos = setup$var_pos,
                 read_cols = read_cols, col_select = col_select,
                 parallel = parallel)

  it <- .chunk_iterator(tasks, .gt_long_chunk, shared, label = "Writing chunk",
                        local = !parallel)
  on.exit(it$done(), add = TRUE)
  while (!is.null(w <- it$next_wave())) {
    for (r in w$res) {
      if (is.null(r)) next
      if (r$first <= last_id)
        cli::cli_abort("Internal error: genotype chunks are not in .row_id order.")
      last_id <- r$last
      if (!parallel) {
        write_part(r$tbl)
      } else {
        if (!file.append(out_file, r$part))
          cli::cli_abort("Could not append to {.file {out_file}}.")
        unlink(r$part)
      }
      r <- NULL
      .release_arrow_memory()
    }
  }

  # no genotypes retained: write an empty table with the same columns
  if (is.null(writer)) write_part(empty)

  cli::cli_alert_success("gt_long table written to {.file {out_file}}")

  invisible(out_file)
}

# Filter one chunk to valid variants and requested samples and sort it
# (variant-major, then sample group order).  Returns NULL if no rows remain,
# else list(first, last: its first and last .row_id, and tbl: the sorted
# Arrow Table, or part: the CSV file it was written to when shared$parallel).
.gt_long_chunk <- function(task, shared) {
  chunk <- arrow::read_feather(task$fpath, col_select = shared$read_cols,
                               as_data_frame = FALSE)
  row_id <- as.vector(chunk$.row_id)
  sample_order <- .sample_index(chunk, shared$samples)
  keep <- which(!is.na(.match_row_id(row_id, shared$var_pos)) & !is.na(sample_order))
  if (length(keep) == 0L) return(NULL)

  o <- keep[order(row_id[keep], sample_order[keep])]
  out <- list(first = row_id[o[1L]], last = row_id[o[length(o)]])
  chunk <- chunk$Take(arrow::Array$create(o - 1L))
  if (!is.null(shared$col_select)) chunk <- chunk[, shared$col_select]

  if (!shared$parallel) {
    out$tbl <- chunk
  } else {
    # same settings as write.csv(), without a header
    utils::write.table(as.data.frame(chunk), task$part, sep = ",", dec = ".",
                       qmethod = "double", row.names = FALSE, col.names = FALSE)
    out$part <- task$part
    chunk <- NULL
    .release_arrow_memory()
  }
  out
}
