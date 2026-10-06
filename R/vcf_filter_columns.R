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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' .vcf_filter_columns(vcf, keep = vcf@groups != "OG")
#'
#' @keywords internal
#' @noRd
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

    tmp_dir <- tempfile("arrow_vcf_samp_")
    dir.create(tmp_dir)

    cli::cli_alert_info(
      "Compacting GT: {length(s)} -> {length(keep_ids)} sample{?s} \\
       across {length(ffiles)} chunk{?s}"
    )
    # Each chunk is written to part_<i>.arrow (in parallel with workers);
    # non-empty parts are then renamed chunk_1, chunk_2, ... in order.
    tasks <- lapply(seq_along(ffiles), function(i) list(
      fpath = ffiles[[i]], out = file.path(tmp_dir, paste0("part_", i, ".arrow"))
    ))
    shared <- list(pos = .row_id_pos(vcf_arrow@variants$.row_id), samples = keep_ids)
    written <- unlist(.map_chunks(tasks, .compact_chunk, shared = shared,
                                  label = "Compacting chunk"))
    parts <- vapply(tasks[written], `[[`, character(1), "out")
    out_idx <- length(parts)
    file.rename(parts, file.path(tmp_dir, paste0("chunk_", seq_len(out_idx), ".arrow")))

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
      invariant_removed = vcf_arrow@invariant_removed,
      loci = vcf_arrow@loci
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

# Write the rows of live variants (shared$pos) and kept samples
# (shared$samples) of one chunk to task$out.  The chunk stays an Arrow Table:
# converting it (including its per-sample FORMAT strings) to R and back is the
# dominant cost.  Returns FALSE (and writes nothing) if no rows remain.

.compact_chunk <- function(task, shared) {
  chunk <- arrow::read_feather(task$fpath, as_data_frame = FALSE)
  in_samples <- as.vector(arrow::call_function(
    "is_in", chunk$sample,
    options = list(value_set = arrow::Array$create(shared$samples), skip_nulls = FALSE)
  ))
  keep <- !is.na(.match_row_id(as.vector(chunk$.row_id), shared$pos)) & in_samples
  chunk <- chunk$Filter(arrow::Array$create(keep))
  written <- chunk$num_rows > 0L
  if (written) arrow::write_feather(chunk, task$out)
  chunk <- NULL
  .release_arrow_memory()
  written
}
