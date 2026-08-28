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

vcf_bind_sparse <- function(...,
                            mode = c("intersect", "union"),
                            absent_as = c("missing", "hom_ref")) {

  # Capture BEFORE match.arg() forces the promise, so we can tell whether the
  # caller explicitly supplied absent_as, vs. it falling through to the
  # default. This only matters for mode = "union" — see below.
  absent_as_supplied <- !missing(absent_as)

  mode <- match.arg(mode)

  if (mode == "union") {
    if (!absent_as_supplied) {
      cli::cli_warn(c(
        "{.arg absent_as} not specified for {.code mode = \"union\"} — \\
         defaulting to {.val missing} (./.) for genotypes absent from a \\
         source object.",
        "i" = "Set {.code absent_as = \"hom_ref\"} explicitly if absent \\
               positions represent homozygous-reference calls (e.g. \\
               independently called VCFs against the same reference)."
      ))
    }
    absent_as <- match.arg(absent_as)
  }
  # mode == "intersect": absent_as is irrelevant here and is intentionally
  # NEVER touched, validated, or evaluated — whatever was (or wasn't) passed
  # for it, including an invalid value, has no effect and causes no error.

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

  # ── target id set ────────────────────────────────────────────────────────
  if (mode == "intersect") {
    target_ids <- sort(
      Reduce(intersect, lapply(vcfs, function(v) as.integer(v@variants$.row_id)))
    )
  } else {
    own_universe <- function(v) {
      removed_ids <- if (is.null(v@invariant_removed) || nrow(v@invariant_removed) == 0L) {
        integer(0)
      } else {
        as.integer(v@invariant_removed$.row_id)
      }
      union(as.integer(v@variants$.row_id), removed_ids)
    }
    target_ids <- sort(Reduce(union, lapply(vcfs, own_universe)))
  }

  n_common <- length(target_ids)
  if (n_common == 0L)
    cli::cli_abort("No variants found to bind (target set is empty).")

  cli::cli_alert_info(
    "Binding {length(vcfs)} VCFArrow object{?s} ({mode}): \\
     {n_common} variant{?s}, {length(all_samples)} total sample{?s}."
  )
  if (mode == "union")
    cli::cli_alert_info("Absent genotypes will be filled as: {absent_as}")

  # ── gt collection ────────────────────────────────────────────────────────
  gt_parts <- vector("list", length(vcfs))
  cli::cli_progress_bar("Reading object", total = length(vcfs))
  for (i in seq_along(vcfs)) {
    gt_parts[[i]] <- if (mode == "intersect") {
      .collect_gt(vcfs[[i]], target_ids)
    } else {
      .collect_gt_union(vcfs[[i]], target_ids, absent_as)
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  gt_parts <- Filter(Negate(is.null), gt_parts)
  if (length(gt_parts) == 0L)
    cli::cli_abort("No genotype data was recovered for the target variants.")

  merged_gt <- do.call(rbind, gt_parts)
  rownames(merged_gt) <- NULL

  sample_pos <- setNames(seq_along(all_samples), all_samples)
  row_ord <- order(merged_gt$.row_id, sample_pos[merged_gt$sample])
  merged_gt <- merged_gt[row_ord, ]
  rownames(merged_gt) <- NULL

  id_map <- setNames(seq_len(n_common), as.character(target_ids))
  merged_gt$.row_id <- unname(id_map[as.character(merged_gt$.row_id)])

  # ── write output feather chunks ───────────────────────────────────────────
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
  rm(merged_gt)

  # ── variant metadata + info ──────────────────────────────────────────────
  if (mode == "intersect") {

    vars_src <- as.data.frame(vcfs[[1]]@variants)
    in_cmn <- vars_src$.row_id %in% target_ids
    vars_out <- vars_src[in_cmn, ]
    vars_out <- vars_out[order(vars_out$.row_id), ]
    vars_out$.row_id <- unname(id_map[as.character(vars_out$.row_id)])
    rownames(vars_out) <- NULL

    info_src <- vcfs[[1]]@info
    if (length(info_src) != nrow(vars_src))
      cli::cli_abort(c(
        "@info length ({length(info_src)}) does not match @variants length \\
         ({nrow(vars_src)}) for the first VCFArrow object.",
        "i" = "One-time repair: vcf@info <- vcf@info[vcf@variants$.row_id]."
      ))
    info_vec <- info_src[in_cmn][order(vars_src$.row_id[in_cmn])]

  } else {
    # union: gather LIVE + RECOVERABLE metadata from every object. rbind()
    # concatenates in vcfs[[1]], vcfs[[2]], ... order, and !duplicated() keeps
    # the FIRST occurrence — so the first object's version wins for any id
    # present in more than one, consistent with intersect mode's convention.
    gather_meta <- function(v) {
      live <- v@variants
      live$.info_str <- v@info  # positionally parallel, per the @info fix
      removed <- v@invariant_removed
      if (!is.null(removed) && nrow(removed) > 0L) rbind(live, removed) else live
    }
    combined <- do.call(rbind, lapply(vcfs, gather_meta))
    combined <- combined[!duplicated(combined$.row_id), , drop = FALSE]
    combined <- combined[combined$.row_id %in% target_ids, , drop = FALSE]
    combined <- combined[order(combined$.row_id), ]

    vars_out <- combined[, setdiff(names(combined), ".info_str"), drop = FALSE]
    vars_out$.row_id <- unname(id_map[as.character(vars_out$.row_id)])
    rownames(vars_out) <- NULL

    info_vec <- combined$.info_str
  }

  # ── FORMAT lookup ──────────────────────────────────────────────────────────
  # @format carries its own .row_id column and is matched by value, never
  # filtered down — always fully recoverable regardless of mode.
  fmt_src <- as.data.frame(vcfs[[1]]@format)
  if (nrow(fmt_src) > 0L && ".row_id" %in% colnames(fmt_src)) {
    fmt_out <- fmt_src[fmt_src$.row_id %in% target_ids, , drop = FALSE]
    fmt_out$.row_id <- unname(id_map[as.character(fmt_out$.row_id)])
    fmt_out <- fmt_out[order(fmt_out$.row_id), ]
    rownames(fmt_out) <- NULL
  } else {
    fmt_out <- fmt_src
  }

  # ── assemble new VCFArrow ───────────────────────────────────────────────────
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
    # invariant_removed intentionally omitted — defaults to empty.
    # Fresh start: everything recoverable from any input object has just
    # been materialized into vars_out, so there is nothing left to track.
  )

  return(new_vcfarrow)
}

##################
# ── Shared helper: read one object's data, filling gaps for union mode ────────
.collect_gt_union <- function(vcf_obj, target_ids, absent_as) {

  ffiles <- list.files(vcf_obj@path, pattern = "\\.arrow$", full.names = TRUE)
  chunk_idx <- as.integer(stringr::str_extract(basename(ffiles), "\\d+"))
  ffiles <- ffiles[order(chunk_idx)]
  samp_keep <- vcf_obj@samples

  parts <- vector("list", length(ffiles))
  for (i in seq_along(ffiles)) {
    df <- arrow::read_feather(
      ffiles[[i]],
      col_select = c(".row_id", "sample", "a1", "a2",
                     "phased", "fmt", "DP", "GQ", "ADR")
    )
    df <- df[df$.row_id %in% target_ids & df$sample %in% samp_keep, , drop = FALSE]
    if (nrow(df) > 0L) parts[[i]] <- df
  }

  non_null <- Filter(Negate(is.null), parts)
  real <- if (length(non_null) == 0L) {
    data.frame(.row_id = integer(0),
               sample = character(0),
               a1 = integer(0),
               a2 = integer(0),
               phased = logical(0),
               fmt = character(0),
               DP = integer(0),
               GQ = integer(0),
               ADR = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    do.call(rbind, non_null)
  }

  # ── Complete the grid ──────────────────────────────────────────────────────
  # Every (target_id, sample) pair for this object's samples should appear
  # exactly once. This single mechanism covers BOTH cases: a variant this
  # object's source VCF never called at all (every sample missing for that
  # id), and the defensive edge case of a partially-compacted object missing
  # only some samples at an id it otherwise has. Either way, the gap gets
  # filled according to absent_as.
  expected_n <- length(target_ids) * length(samp_keep)
  if (nrow(real) < expected_n) {
    have_key <- paste(real$.row_id, real$sample, sep = "\r")
    full_grid <- expand.grid(.row_id = target_ids,
                             sample = samp_keep,
                             KEEP.OUT.ATTRS = FALSE)
    want_key <- paste(full_grid$.row_id, full_grid$sample, sep = "\r")
    gap_mask <- !(want_key %in% have_key)

    if (any(gap_mask)) {
      fill_a <- if (absent_as == "hom_ref") 0L else NA_integer_
      synth <- data.frame(
        .row_id = full_grid$.row_id[gap_mask],
        sample = full_grid$sample[gap_mask],
        a1 = fill_a,
        a2 = fill_a,
        phased = NA,
        fmt = NA_character_,
        DP = NA_integer_,
        GQ = NA_integer_,
        ADR = NA_real_,
        stringsAsFactors = FALSE
      )
      real <- rbind(real, synth)
    }
  }

  real
}

# ── Helper: read one VCFArrow's feather files, filter to target row_ids ───────
# for mode = "intersect"
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
