#' @title vcf_bind_sparse
#'
#' @description
#' Bind two or more VCFArrow objects into a new VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param ... -> a collection of VCFArrow objects
#'
#' @return VCFArrow object
#'
#' @details
#' This function binds two or more VCFArrow objects, returning new VCFArrow object.
#' The VCFArrow objects need not have the same SNPs, and must have unique individuals.
#'
#' @examples
#' vcf_bind_sparse(my_vcf1, my_vcf2, other_vcf, ...)
#'
#' @export
#'

vcf_bind_sparse <- function(...) {

  vcfs <- list(...)

  # input validation
  if (length(vcfs) < 2L)
    cli::cli_abort("Provide at least two VCFArrow objects to bind.")

  for (i in seq_along(vcfs))
    if (!inherits(vcfs[[i]], "VCFArrow"))
      cli::cli_abort("Argument {i} is not a VCFArrow object.")

  # sample name collision check
  all_samples <- unlist(lapply(vcfs, function(v) v@samples))
  all_groups <- unlist(lapply(vcfs, function(v) v@groups))
  dup <- all_samples[duplicated(all_samples)]
  if (length(dup))
    cli::cli_abort(c(
      "Duplicate sample names found across VCFArrow objects.",
      "x" = "Duplicate(s): {.val {unique(dup)}}"
    ))

  # variant intersection
  common_ids <- sort(
    Reduce(intersect, lapply(vcfs, function(v) as.integer(v@variants$.row_id)))
  )
  n_common <- length(common_ids)

  if (n_common == 0L)
    cli::cli_abort("No common variants found between VCFArrow objects.")

  cli::cli_alert_info(
    "Binding {length(vcfs)} VCFArrow object{?s}: \\
     {n_common} common variant{?s}, {length(all_samples)} total sample{?s}."
  )

  # ── Helper: read one VCFArrow's feather files, filter to target row_ids ───────
  # Returns a plain data.frame — no Arrow lazy objects.
  .collect_gt <- function(vcf_obj, target_ids) {
    ffiles <- list.files(vcf_obj@path, pattern = "\\.arrow$", full.names = TRUE)
    chunk_idx <- as.integer(stringr::str_extract(basename(ffiles), "\\d+"))
    ffiles <- ffiles[order(chunk_idx)]
    samp_keep <- vcf_obj@samples

    parts <- vector("list", length(ffiles))
    for (i in seq_along(ffiles)) {
      # Read only the four columns needed for genotype assembly
      df <- arrow::read_feather(
        ffiles[[i]],
        col_select = c(".row_id", "sample", "a1", "a2",
                       "phased", "fmt", "DP", "GQ", "ADR")
      )
      df <- df[df$.row_id %in% target_ids & df$sample %in% samp_keep, , drop = FALSE]
      if (nrow(df) > 0L) parts[[i]] <- df
    }
    non_null <- Filter(Negate(is.null), parts)
    if (length(non_null) == 0L) return(NULL)
    do.call(rbind, non_null)
  }

  # collect gt data from every VCFArrow object
  # NOTE: cli_progress_update() must run in the *same* environment as
  # cli_progress_bar(). Use a plain for loop to keep everything in one environment.
  gt_parts <- vector("list", length(vcfs))
  cli::cli_progress_bar("Reading object", total = length(vcfs))
  for (i in seq_along(vcfs)) {
    gt_parts[[i]] <- .collect_gt(vcfs[[i]], common_ids)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # remove any NULLs (VCFArrow had no data for common_ids; should not happen)
  gt_parts <- Filter(Negate(is.null), gt_parts)
  if (length(gt_parts) == 0L)
    cli::cli_abort("No genotype data were recovered for the common variants.")

  # combine data.frames
  merged_gt <- do.call(rbind, gt_parts)
  rownames(merged_gt) <- NULL

  # sort: variant-major, then by global sample position
  sample_pos <- setNames(seq_along(all_samples), all_samples)
  row_ord <- order(merged_gt$.row_id, sample_pos[merged_gt$sample])
  merged_gt <- merged_gt[row_ord, ]
  rownames(merged_gt) <- NULL

  # remap .row_ids to 1 .. n_common
  # unname() is essential here: id_map is a NAMED vector (setNames(...)), and
  # indexing a named vector by character subscript preserves those names on
  # the result.
  id_map <- setNames(seq_len(n_common), as.character(common_ids))
  merged_gt$.row_id <- unname(id_map[as.character(merged_gt$.row_id)])

  # write output feather chunks
  tmp_dir <- tempfile("arrow_vcf_bind_")
  dir.create(tmp_dir)
  chunk_size <- 50000L
  n_chunks <- ceiling(n_common / chunk_size)

  cli::cli_progress_bar("Writing chunk", total = n_chunks)
  for (ci in seq_len(n_chunks)) {
    id_lo <- (ci - 1L) * chunk_size + 1L
    id_hi <- min(ci * chunk_size, n_common)
    mask <- merged_gt$.row_id >= id_lo & merged_gt$.row_id <= id_hi
    arrow::write_feather(
      merged_gt[mask, ],
      file.path(tmp_dir, paste0("chunk_", ci, ".arrow"))
    )
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  rm(merged_gt)  # free memory before building metadata

  # variant metadata
  # source: first object (most-filtered, authoritative variant set).
  vars_src <- as.data.frame(vcfs[[1]]@variants)
  in_cmn <- vars_src$.row_id %in% common_ids
  vars_out <- vars_src[in_cmn, ]
  vars_out <- vars_out[order(vars_out$.row_id), ]
  vars_out$.row_id  <- unname(id_map[as.character(vars_out$.row_id)])
  rownames(vars_out) <- NULL

  # info vector
  # @info is parallel to @variants (same length, same order)
  info_src <- vcfs[[1]]@info
  if (length(info_src) != nrow(vars_src))
    cli::cli_abort(c(
      "@info length ({length(info_src)}) does not match @variants length \\
       ({nrow(vars_src)}) for the first VCFArrow object.",
      "i" = "@info appears to be misaligned with @variants — likely from \\
             filtering performed before the @info alignment fix in \\
             vcf_filter_rows() was applied.",
      "i" = "One-time repair: vcf@info <- vcf@info[vcf@variants$.row_id], \\
             then retry vcf_bind_sparse()."
    ))

  info_vec <- info_src[in_cmn][order(vars_src$.row_id[in_cmn])]

  # FORMAT lookup
  fmt_src <- as.data.frame(vcfs[[1]]@format)
  if (nrow(fmt_src) > 0L && ".row_id" %in% colnames(fmt_src)) {
    fmt_out <- fmt_src[fmt_src$.row_id %in% common_ids, , drop = FALSE]
    fmt_out$.row_id <- unname(id_map[as.character(fmt_out$.row_id)])
    fmt_out <- fmt_out[order(fmt_out$.row_id), ]
    rownames(fmt_out) <- NULL
  } else {
    fmt_out <- fmt_src
  }

  # assemble new VCFArrow
  # suppressWarnings: Arrow emits "Invalid metadata$r" when it re-parses the
  # R-specific IPC schema annotations that write_feather() embeds in new files
  # benign - no data affected
  gt_arrow <- suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather"))

  new_vcfarrow <- .new_vcfarrow(
    header = vcfs[[1]]@header,
    info = info_vec,
    format = fmt_out,
    variants = vars_out,
    gt = gt_arrow,
    samples = all_samples,
    groups = all_groups,
    path = tmp_dir
  )

  return(new_vcfarrow)
}
