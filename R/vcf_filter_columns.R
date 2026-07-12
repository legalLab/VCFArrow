#' @title vcf_filter_columns
#'
#' @description
#' Unified API for filtering VCFArrow objects by row IDs.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep -> rows to keep (numeric or logical)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report removed samples and final % missing data, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes all rows in the 'keep' parameter from a VCFArrow object,
#' returning a new VCFArrow object.
#' Optionally will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' Optionally will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_filter_rank(vcf_arrow = my_vcf, keep = rows_to_keep)
#' vcf_filter_rank(my_vcf, rows_to_keep)
#'
#' @export
#'

.vcf_filter_columns <- function(vcf_arrow, keep, f_invar = TRUE, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  s <- vcf_arrow@samples

  keep_ids <- if (is.logical(keep)) {
    if (length(keep) != length(s))
      cli::cli_abort("Logical filter length must match number of samples")
    s[keep]
  } else if (is.numeric(keep)) {
    s[keep]
  } else if (is.character(keep)) {
    intersect(s, keep)
  } else {
    cli::cli_abort("{.arg keep} must be character, logical, or numeric")
  }

  removed <- setdiff(s, keep_ids)
  idx <- match(keep_ids, s)

  if (length(removed) > 0L) {

    # ── Compact feather files ─────────────────────────────────────────────────
    # Read each chunk, retain only kept samples AND currently-live variants,
    # write to a new temp directory. This produces the smallest possible
    # backing files, benefiting every downstream operation.
    ffiles <- .get_sorted_feather_files(vcf_arrow@path)

    # Fast logical-vector lookup for row_id membership (same pattern used in
    # all filter scan loops throughout this package).
    valid_row_ids <- vcf_arrow@variants$.row_id
    max_id <- if (length(valid_row_ids) > 0L) max(valid_row_ids) else 0L
    lv <- logical(max_id + 1L)
    lv[valid_row_ids] <- TRUE

    tmp_dir <- tempfile("arrow_vcf_samp_")
    dir.create(tmp_dir)
    out_idx <- 0L

    cli::cli_alert_info(
      "Compacting GT: {length(s)} -> {length(keep_ids)} sample{?s} \\
       across {length(ffiles)} chunk{?s}"
    )
    cli::cli_progress_bar("Compacting chunk", total = length(ffiles))

    for (fpath in ffiles) {
      chunk <- arrow::read_feather(fpath)
      chunk <- chunk[lv[chunk$.row_id] & chunk$sample %in% keep_ids, , drop = FALSE]
      if (nrow(chunk) > 0L) {
        out_idx <- out_idx + 1L
        arrow::write_feather(
          chunk,
          file.path(tmp_dir, paste0("chunk_", out_idx, ".arrow"))
        )
      }
      chunk <- NULL; gc(verbose = FALSE, full = FALSE)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    if (out_idx == 0L)
      cli::cli_abort("No genotype data remained after sample filtering.")

    gt_arrow <- suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather"))

    # Build a new VCFArrow through the canonical constructor so registration,
    # finalizer, and invariant_removed tracking are all handled correctly.
    vcf_arrow <- .new_vcfarrow(
      header = vcf_arrow@header,
      info = vcf_arrow@info,
      format = vcf_arrow@format,
      variants = vcf_arrow@variants,
      gt = gt_arrow,
      samples = keep_ids,
      groups = vcf_arrow@groups[idx],
      path = tmp_dir,
      invariant_removed = vcf_arrow@invariant_removed
    )

  } else {
    # No samples removed — pure metadata update, no I/O needed.
    vcf_arrow@samples <- keep_ids
    vcf_arrow@groups <- vcf_arrow@groups[idx]
  }

  # only filter invariants if samples were removed
  if (f_invar && length(removed) > 0L) vcf_arrow <- vcf_filter_invariant(vcf_arrow)

  if (verbose) {
    if (length(removed) > 0L)
      cli::cli_alert_info("Removed samples: {removed}")
    cli::cli_alert_info(
      "Variants retained: {nrow(vcf_arrow@variants)} | \\
       Samples retained: {length(vcf_arrow@samples)}"
    )
  }

  return(vcf_arrow)
}
