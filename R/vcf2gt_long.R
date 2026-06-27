#' @title vcf2gt_long
#'
#' @description
#' Converts a VCFArrow object to tidy long format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'gt_long' (character)
#' @param format -> one of three output formats (arrow, parquet, CSV) (character)
#' @param col_select -> optional selection of columns to save, default ALL
#'
#' @return NULL
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
#' vcf2gt_long(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "gt_long", format = "csv")
#' vcf2gt_long(vcf_arrow, my_groups, format = "csv")
#' vcf2gt_long(vcf_arrow, format = "csv")
#'
#' @export
#'

vcf2gt_long <- function(vcf_arrow, keep_groups = NULL,
                        out_file = "gt_long",
                        format = c("feather", "parquet", "csv"),
                        col_select = NULL) {

  format <- match.arg(format)
  setup <-.vcf_export_setup(vcf_arrow, keep_groups)

  # Resolve output path (add extension if missing)
  ext_map <- c(feather = ".feather", parquet = ".parquet", csv = ".csv")
  ext <- ext_map[[format]]
  if (!endsWith(out_file, ext)) {
    out_file <- paste0(out_file, ext)
  }

  cli::cli_alert_info(
    "Building gt_long: {setup$n_var} variants x {setup$n_samples} samples \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)"
  )
  cli::cli_progress_bar("Reading chunk", total = length(setup$feather_files))

  chunks <- vector("list", length(setup$feather_files))

  for (i in seq_along(setup$feather_files)) {
    fpath <- setup$feather_files[[i]]

    chunk <- if (is.null(col_select)) {
      arrow::read_feather(fpath)
    } else {
      arrow::read_feather(fpath, col_select = col_select)
    }

    # Filter to valid variants and requested samples
    chunk <- chunk[chunk$.row_id %in% setup$valid_row_ids &
                     chunk$sample %in% setup$samples, ]

    chunks[[i]] <- chunk
    cli::cli_progress_update()
  }

  cli::cli_progress_done()
  cli::cli_alert_info("Combining and writing GT long table...")

  combined <- do.call(rbind, chunks)

  # Sort: variant-major, then by sample group order
  sample_order <- match(combined$sample, setup$samples)
  combined <- combined[order(combined$.row_id, sample_order), ]
  rownames(combined) <- NULL

  switch(format,
         feather = arrow::write_feather(combined, out_file),
         parquet = arrow::write_parquet(combined, out_file),
         csv = utils::write.csv(combined, out_file, row.names = FALSE)
  )

  cli::cli_alert_success("gt_long table written to {.file {out_file}}")

  invisible(out_file)
}
