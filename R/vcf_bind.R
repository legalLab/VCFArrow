#' @title vcf_bind
#'
#' @description
#' Bind two or more VCFArrow objects into a new VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param ... -> a collection of VCFArrow objects
#' @param mode -> "intersect" keeps variants present in all objects, "union"
#'   keeps variants present in any object, default "intersect" (character)
#' @param absent_as -> for mode = "union" or outgroup = TRUE, how genotypes
#'   absent from a source object are coded, "missing" or "hom_ref",
#'   default "missing" (character)
#' @param outgroup -> add the other objects (e.g. outgroups) to the first
#'   object (e.g. a filtered ingroup), keeping exactly the variants of the
#'   first object; mode and recover_loci are then not used, default FALSE
#'   (Boolean)
#' @param recover_loci -> for mode = "union", restore the genotypes of
#'   variants an object removed by filtering, where its genotype files still
#'   hold them, instead of coding them according to absent_as, default FALSE
#'   (Boolean)
#' @param .check -> deprecated and ignored; objects no longer need to share
#'   the same variants
#'
#' @return VCFArrow object
#'
#' @details
#' This function binds two or more VCFArrow objects with different sets of
#' individuals, returning a new VCFArrow object. Sample names must be unique
#' across objects.
#' The objects need not have the same variants: mode = "intersect" keeps the
#' variants present in all objects, mode = "union" keeps the variants present
#' in any object and codes the genotypes an object lacks according to
#' absent_as.
#' For objects with the same variants both modes give the same result.
#' In mode = "union", variants an object removed stay removed for it: a
#' variant that one object removed by filtering but that another object
#' contains is included, with the genotypes of the object that removed it
#' coded according to absent_as. Binding in steps and binding all objects in
#' one call then give the same result. With recover_loci = TRUE, the real
#' genotypes of such variants are restored instead, read from the genotype
#' files of the object that removed them (which keep every variant the object
#' had when its files were written). This suits joining separately filtered
#' datasets; binding in steps can then differ from binding in one call.
#'
#' outgroup = TRUE is meant for adding outgroups to an ingroup that was
#' filtered on its own (so that properties of the outgroup, such as its
#' missingness, cannot remove ingroup variants): the result keeps exactly the
#' variants of the first object, and variants removed from it stay removed.
#' The other objects contribute their genotypes for these variants, read from
#' their genotype files wherever present (including variants they removed,
#' for example as invariant); genotypes they do not have are coded according
#' to absent_as.
#'
#' Variants are matched by CHROM, POS, REF and ALT, so the objects may come
#' from different VCF files (e.g. separate sequencing runs); REF and ALT must
#' be written the same way in all of them (including the order of alternative
#' alleles). This is assured when all datasets were called against the same
#' reference genome assembly, which gives the same CHROM, POS and REF for the
#' same locus. It is generally not possible for de novo (reference-free) SNP
#' calling, for example with DiscoSnp, where locus names and positions are
#' generated anew in each calling run, so the same locus will usually not
#' match between separately called datasets; such datasets should be called
#' together in a single run instead. Multi-allelic sites can list their ALT
#' alleles in a different order in different files, so it is best to
#' normalize the VCFs first, for example by splitting multi-allelic sites with
#' `bcftools norm -m -any`. Variants keep the .row_id of the first object where its source
#' VCF has them, other variants get new ids; the variant metadata of the first
#' object and, in mode = "intersect", its record of removed invariant variants
#' are kept, so bound objects can be bound again.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' ingroup <- vcf_extract_groups(vcf, c("GS", "BS", "WA"))
#' outgroup <- vcf_extract_groups(vcf, "OG")
#' vcf_bind(ingroup, outgroup, mode = "intersect")
#' vcf_bind(ingroup, outgroup, mode = "union", absent_as = "missing")
#' vcf_bind(vcf_filter_missingness(ingroup, 0.2), outgroup, outgroup = TRUE,
#'          absent_as = "missing")
#'
#' @export
#'

vcf_bind <- function(...,
                     mode = c("intersect", "union"),
                     absent_as = c("missing", "hom_ref"),
                     outgroup = FALSE,
                     recover_loci = FALSE,
                     .check = NULL) {

  if (!is.null(.check))
    cli::cli_warn(c(
      "{.arg .check} is deprecated and ignored.",
      "i" = "Objects no longer need to share the same variants; see {.arg mode}."
    ))

  # Capture BEFORE match.arg() forces the promise, so we can tell whether the
  # caller explicitly supplied absent_as, vs. it falling through to the
  # default. This only matters for mode = "union" and outgroup = TRUE — see
  # below.
  absent_as_supplied <- !missing(absent_as)

  if (!isTRUE(outgroup) && !isFALSE(outgroup))
    cli::cli_abort("{.arg outgroup} must be TRUE or FALSE.")
  if (!isTRUE(recover_loci) && !isFALSE(recover_loci))
    cli::cli_abort("{.arg recover_loci} must be TRUE or FALSE.")
  if (outgroup && !missing(mode))
    cli::cli_warn("{.arg mode} is not used when {.code outgroup = TRUE}.")

  # which variants to keep: those of all objects ("intersect"), of any object
  # ("union"), or of the first object ("outgroup")
  mode <- if (outgroup) "outgroup" else match.arg(mode)
  if (recover_loci && mode != "union")
    cli::cli_warn("{.arg recover_loci} is only used with {.code mode = \"union\"}.")

  if (mode != "intersect") {
    if (!absent_as_supplied) {
      setting <- if (outgroup) "outgroup = TRUE" else "mode = \"union\""
      cli::cli_warn(c(
        "{.arg absent_as} not specified for {.code {setting}} - \\
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

  # ── target variants, matched by locus (CHROM, POS, REF, ALT) ──────────────
  # Output .row_ids: the first object's .row_id for loci its source VCF has
  # (see @loci), so binding objects from one read_vcf() keeps their .row_ids;
  # loci only found in other objects get new ids after those.
  first <- vcfs[[1]]
  first_vars <- as.data.frame(first@variants)

  if (mode %in% c("intersect", "outgroup")) {
    keep <- rep(TRUE, nrow(first_vars))
    if (mode == "intersect")
      for (v in vcfs[-1])
        keep <- keep & vctrs::vec_in(.locus_key(first_vars), .locus_key(v@variants))
    ord <- order(first_vars$.row_id[keep])
    vars_out <- first_vars[keep, , drop = FALSE][ord, , drop = FALSE]

    info_src <- first@info
    if (length(info_src) != nrow(first_vars))
      cli::cli_abort(c(
        "@info length ({length(info_src)}) does not match @variants length \\
         ({nrow(first_vars)}) for the first VCFArrow object.",
        "i" = "One-time repair: vcf@info <- vcf@info[vcf@variants$.row_id]."
      ))
    info_vec <- info_src[keep][ord]
    format_out <- first@format
    loci_out <- first@loci

    # variants the first object removed as invariant stay recoverable
    invariant_removed <- first@invariant_removed

  } else {
    # union: gather LIVE + RECOVERABLE metadata from every object. rbind()
    # concatenates in vcfs[[1]], vcfs[[2]], ... order, and the first
    # occurrence of each locus is kept — so the first object's version wins
    # for any locus present in more than one, consistent with intersect mode.
    combined <- do.call(rbind, lapply(seq_along(vcfs), function(k) {
      v <- vcfs[[k]]
      live <- as.data.frame(v@variants)
      live$.info_str <- v@info  # positionally parallel, per the @info fix
      removed <- v@invariant_removed
      m <- if (!is.null(removed) && nrow(removed) > 0L) rbind(live, removed) else live
      m$.src <- rep(k, nrow(m))
      m
    }))
    combined <- combined[vctrs::vec_unique_loc(.locus_key(combined)), , drop = FALSE]

    # output ids; loci unknown to the first object's source are appended
    out_id <- vctrs::vec_match(.locus_key(combined), .locus_key(first@loci))
    new_loci <- is.na(out_id)
    out_id[new_loci] <- nrow(first@loci) + seq_len(sum(new_loci))
    src_id <- combined$.row_id
    combined$.row_id <- as.integer(out_id)
    ord <- order(combined$.row_id)
    combined <- combined[ord, , drop = FALSE]
    src_id <- src_id[ord]
    new_loci <- new_loci[ord]

    vars_out <- combined[, setdiff(names(combined), c(".info_str", ".src")), drop = FALSE]
    rownames(vars_out) <- NULL
    info_vec <- combined$.info_str

    # FORMAT and locus lookups: the first object's, plus the appended loci
    # (their FORMAT taken from the object they came from)
    format_out <- first@format
    loci_out <- first@loci
    if (any(new_loci)) {
      add <- combined[new_loci, , drop = FALSE]
      add_src_id <- src_id[new_loci]
      add_fmt <- character(nrow(add))
      for (k in unique(add$.src)) {
        r <- add$.src == k
        f <- vcfs[[k]]@format
        add_fmt[r] <- f$FORMAT[.match_row_id(add_src_id[r], .row_id_pos(f$.row_id))]
      }
      format_out <- rbind(format_out,
                          data.frame(FORMAT = add_fmt, .row_id = add$.row_id,
                                     stringsAsFactors = FALSE))
      loci_out <- rbind(loci_out, .locus_key(add))
    }

    # every object's removed invariant variants are live again in vars_out,
    # so there is nothing left to track
    invariant_removed <- NULL
  }

  n_common <- nrow(vars_out)
  if (n_common == 0L)
    cli::cli_abort("No variants found to bind (target set is empty).")

  cli::cli_alert_info(
    "Binding {length(vcfs)} VCFArrow object{?s} ({mode}): \\
     {n_common} variant{?s}, {length(all_samples)} total sample{?s}."
  )
  if (mode != "intersect")
    cli::cli_alert_info("Absent genotypes will be filled as: {absent_as}")

  # ── genotypes ───────────────────────────────────────────────────────────
  # Each object's own .row_id for each target locus (NA: it has no genotypes
  # to contribute, filled as absent_as in union / outgroup mode):
  #   intersect: its variants (every object has all target loci)
  #   union: its own variants only (removed ones stay removed), unless
  #     recover_loci: any locus of its source VCF, read from its genotype
  #     files wherever they still hold it
  #   outgroup: any locus of its source VCF (so loci it removed, e.g. as
  #     invariant, still get its real genotypes)
  target_key <- .locus_key(vars_out)
  obj_ids <- lapply(vcfs, function(v) {
    if (mode == "intersect" || (mode == "union" && !recover_loci)) {
      own <- as.data.frame(v@variants)
      if (mode == "union" && nrow(v@invariant_removed) > 0L)
        own <- rbind(own[c(.locus_cols, ".row_id")],
                     as.data.frame(v@invariant_removed)[c(.locus_cols, ".row_id")])
      own$.row_id[vctrs::vec_match(target_key, .locus_key(own))]
    } else {
      vctrs::vec_match(target_key, .locus_key(v@loci))  # row i of @loci = .row_id i
    }
  })
  tmp_dir <- .bind_write_gt(vcfs, obj_ids, vars_out$.row_id,
                            if (mode == "intersect") "intersect" else "union",
                            absent_as)

  # ── assemble new VCFArrow ───────────────────────────────────────────────────
  gt_arrow <- suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather"))

  .new_vcfarrow(
    header = first@header,
    info = info_vec,
    format = format_out,
    variants = vars_out,
    gt = gt_arrow,
    samples = all_samples,
    groups = all_groups,
    path = tmp_dir,
    invariant_removed = invariant_removed,
    loci = loci_out
  )
}

#' @title vcf_bind_sparse
#'
#' @description
#' Deprecated: use vcf_bind(), which now binds objects with different
#' variants (see its mode and absent_as arguments). vcf_bind_sparse(...) is
#' vcf_bind(..., recover_loci = TRUE).
#'
#' @param ... -> a collection of VCFArrow objects
#' @param mode -> see vcf_bind()
#' @param absent_as -> see vcf_bind()
#'
#' @return VCFArrow object
#'
#' @keywords internal
#' @export
#'

vcf_bind_sparse <- function(...,
                            mode = c("intersect", "union"),
                            absent_as = c("missing", "hom_ref")) {
  cli::cli_warn(c(
    "{.fn vcf_bind_sparse} is deprecated; use {.fn vcf_bind} instead.",
    "i" = "{.fn vcf_bind} takes the same {.arg mode} and {.arg absent_as} arguments."
  ), .frequency = "regularly", .frequency_id = "VCFArrow_vcf_bind_sparse")
  # forward absent_as only if supplied, so vcf_bind() still warns about an
  # unspecified absent_as in union mode
  # recover_loci = TRUE: the behaviour vcf_bind_sparse() always had
  if (missing(absent_as)) vcf_bind(..., mode = mode, recover_loci = TRUE)
  else vcf_bind(..., mode = mode, absent_as = absent_as, recover_loci = TRUE)
}

# Write the genotypes of `vcfs` for the target variants to a new chunk
# directory, samples in binding order, and return its path.  obj_ids[[k]]
# gives, for each target variant (in output order), the .row_id of that
# variant in object k, or NA if k contributes no genotypes for it; output
# .row_ids are out_row_ids.  mode = "union" fills (variant, sample) pairs
# absent from an object according to absent_as; "intersect" leaves them
# absent.
#
# Output chunk ci holds window lo..hi of the target variants.  For each
# window, each object's rows are read only from its feather files whose
# .row_id range overlaps the object's ids in the window, kept as Arrow
# Tables, so memory is O(window) rather than O(whole dataset).

.bind_write_gt <- function(vcfs, obj_ids, out_row_ids, mode, absent_as = NULL) {
  n_common <- length(out_row_ids)
  tmp_dir <- tempfile("arrow_vcf_bind_")
  dir.create(tmp_dir)
  chunk_size <- 50000L
  n_chunks <- ceiling(n_common / chunk_size)

  # Windows are independent (each writes its own chunk file), so they run in
  # parallel with workers.  `state` holds the per-object sources (file
  # ranges, id maps, file cache), built on first use: serially one state
  # serves all windows, so each file is read once; each worker builds its own.
  shared <- list(
    objs = lapply(seq_along(vcfs), function(k) list(
      files = .get_sorted_feather_files(vcfs[[k]]@path),
      samples = vcfs[[k]]@samples, obj_ids = obj_ids[[k]])),
    sample_offset = cumsum(c(0L, lengths(lapply(vcfs, function(v) v@samples)))),
    out_row_ids = out_row_ids, mode = mode, absent_as = absent_as,
    tmp_dir = tmp_dir, state = new.env(parent = emptyenv())
  )
  tasks <- lapply(seq_len(n_chunks), function(ci) list(
    ci = ci, win_lo = (ci - 1L) * chunk_size + 1L, win_hi = min(ci * chunk_size, n_common)
  ))
  n_rows <- sum(unlist(.map_chunks(tasks, .bind_window_task, shared = shared,
                                   label = "Writing chunk")))

  if (n_rows == 0) {
    unlink(tmp_dir, recursive = TRUE)
    cli::cli_abort("No genotype data was recovered for the target variants.")
  }
  tmp_dir
}

# Write output chunk task$ci (target window task$win_lo..task$win_hi) of
# .bind_write_gt(); returns its number of rows.
.bind_window_task <- function(task, shared) {
  st <- shared$state
  if (is.null(st$srcs)) {
    st$srcs <- lapply(shared$objs, function(o) .bind_source(o$files, o$samples, o$obj_ids))
    st$canon <- .bind_schema(st$srcs)
  }
  srcs <- st$srcs
  canon <- st$canon
  win_lo <- task$win_lo
  win_hi <- task$win_hi

  parts <- list(); new_id <- list(); gpos <- list()
  for (k in seq_along(srcs)) {
    w <- .bind_window_rows(srcs[[k]], win_lo, win_hi, canon)
    if (shared$mode == "union")
      w <- .bind_fill_gaps(w, win_lo, win_hi, srcs[[k]]$samples,
                           shared$absent_as, canon)
    if (!is.null(w)) {
      parts[[length(parts) + 1L]] <- w$tbl
      new_id[[length(new_id) + 1L]] <- w$new_id
      gpos[[length(gpos) + 1L]] <- w$s + shared$sample_offset[k]
    }
  }

  # sort: variant-major, then samples in binding order; write output .row_ids
  out <- if (length(parts)) do.call(arrow::concat_tables, parts) else NULL
  n_rows <- 0
  if (!is.null(out) && out$num_rows > 0L) {
    new_id <- unlist(new_id, use.names = FALSE)
    o <- order(new_id, unlist(gpos, use.names = FALSE))
    out <- out$Take(arrow::Array$create(o - 1L))
    out[[".row_id"]] <- as.integer(shared$out_row_ids[new_id[o]])
    n_rows <- out$num_rows
  } else {
    out <- .bind_empty_table(canon)
  }
  arrow::write_feather(out, file.path(shared$tmp_dir, paste0("chunk_", task$ci, ".arrow")))
  parts <- out <- NULL
  .release_arrow_memory()
  n_rows
}

##################
# ── Streaming genotype helpers ────────────────────────────────────────────────

# Variants are identified across objects by locus
.locus_cols <- c("CHROM", "POS", "REF", "ALT")

.locus_key <- function(df) {
  key <- as.data.frame(df)[.locus_cols]
  rownames(key) <- NULL
  key
}

# Per-object source: its feather files, each file's .row_id range, the
# object's .row_id for each target position (obj_ids) and the reverse map
# (id_pos: .row_id -> target position), and a cache (environment, so it
# persists across windows) of files already read.
.bind_source <- function(files, samples, obj_ids) {
  rng <- vapply(files, function(f) {
    id <- as.vector(arrow::read_feather(f, col_select = ".row_id",
                                        as_data_frame = FALSE)$.row_id)
    if (length(id)) as.numeric(range(id)) else c(Inf, -Inf)
  }, numeric(2), USE.NAMES = FALSE)
  has <- which(!is.na(obj_ids))
  id_pos <- integer(if (length(has)) max(obj_ids[has]) else 0L)
  id_pos[obj_ids[has]] <- has
  src <- new.env(parent = emptyenv())
  src$files <- files
  src$lo <- rng[1L, ]
  src$hi <- rng[2L, ]
  src$samples <- samples
  src$obj_ids <- obj_ids
  src$id_pos <- id_pos
  src$cache <- list()
  src
}

.bind_gt_cols <- c(".row_id", "sample", "a1", "a2", "phased", "fmt", "DP", "GQ", "ADR")

# Common output schema (column types of the first object, without R metadata)
.bind_schema <- function(srcs) {
  sch <- arrow::read_feather(srcs[[1L]]$files[1L], col_select = .bind_gt_cols,
                             as_data_frame = FALSE)$schema
  types <- lapply(.bind_gt_cols, function(n) sch$GetFieldByName(n)$type)
  do.call(arrow::schema, stats::setNames(types, .bind_gt_cols))
}

.bind_empty_table <- function(canon) {
  arrow::arrow_table(
    .row_id = integer(0), sample = character(0), a1 = integer(0), a2 = integer(0),
    phased = logical(0), fmt = character(0), DP = numeric(0), GQ = numeric(0),
    ADR = numeric(0), schema = canon
  )
}

# One object's rows for the target window win_lo..win_hi, restricted to its
# ids for target variants (src$obj_ids) and to its own samples.
# Returns list($tbl, $new_id: target position, $s: sample position within the
# object) or NULL.
.bind_window_rows <- function(src, win_lo, win_hi, canon) {
  ids <- src$obj_ids[win_lo:win_hi]
  if (all(is.na(ids))) return(NULL)
  id_lo <- min(ids, na.rm = TRUE)
  id_hi <- max(ids, na.rm = TRUE)
  need <- src$files[src$lo <= id_hi & src$hi >= id_lo]
  # keep only files needed now: with objects from one source (ids in target
  # order), windows move forward through the files and each is read once;
  # otherwise a file may be read again for a later window
  src$cache <- src$cache[intersect(names(src$cache), need)]

  parts <- lapply(need, function(f) {
    if (is.null(src$cache[[f]])) {
      tbl <- arrow::read_feather(f, col_select = .bind_gt_cols, as_data_frame = FALSE)
      new_id <- .match_row_id(as.vector(tbl$.row_id), src$id_pos)
      s <- .sample_index(tbl, src$samples)
      keep <- which(!is.na(new_id) & !is.na(s))
      src$cache[[f]] <- list(tbl = tbl$Take(arrow::Array$create(keep - 1L))$cast(canon),
                             new_id = new_id[keep], s = s[keep])
    }
    ch <- src$cache[[f]]
    in_win <- which(ch$new_id >= win_lo & ch$new_id <= win_hi)
    if (!length(in_win)) return(NULL)
    list(tbl = ch$tbl$Take(arrow::Array$create(in_win - 1L)),
         new_id = ch$new_id[in_win], s = ch$s[in_win])
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  list(tbl = do.call(arrow::concat_tables, lapply(parts, `[[`, "tbl")),
       new_id = unlist(lapply(parts, `[[`, "new_id"), use.names = FALSE),
       s = unlist(lapply(parts, `[[`, "s"), use.names = FALSE))
}

# mode = "union": every (variant, sample) pair of the window must appear once
# for this object's samples.  This covers both a variant the object's source
# VCF never called and the defensive case of a partially-compacted object
# missing only some samples at an id it otherwise has.  Gaps are filled
# according to absent_as.
.bind_fill_gaps <- function(w, win_lo, win_hi, samples, absent_as, canon) {
  n_s <- length(samples)
  n_cell <- (win_hi - win_lo + 1L) * n_s
  present <- logical(n_cell)
  if (!is.null(w)) present[(w$new_id - win_lo) * n_s + w$s] <- TRUE
  miss <- which(!present)
  if (!length(miss)) return(w)

  fill_a <- if (absent_as == "hom_ref") 0L else NA_integer_
  # per-sample FORMAT text, as write_vcf() writes it: "0/0", or missing
  # (written as ".")
  fill_fmt <- if (absent_as == "hom_ref") "0/0" else NA_character_
  new_id <- win_lo + (miss - 1L) %/% n_s
  s <- (miss - 1L) %% n_s + 1L
  n <- length(miss)
  synth <- arrow::arrow_table(
    .row_id = new_id, sample = samples[s], a1 = rep(fill_a, n), a2 = rep(fill_a, n),
    phased = rep(NA, n), fmt = rep(fill_fmt, n), DP = rep(NA_real_, n),
    GQ = rep(NA_real_, n), ADR = rep(NA_real_, n),
    schema = canon
  )
  if (is.null(w)) return(list(tbl = synth, new_id = new_id, s = s))
  list(tbl = arrow::concat_tables(w$tbl, synth),
       new_id = c(w$new_id, new_id), s = c(w$s, s))
}
