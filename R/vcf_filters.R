#' @noRd

# vcf_filters.R
#
# ── Shared chunk-reading helpers ───────────────────────────────────────────────

.get_sorted_feather_files <- function(path) {
  files <- list.files(path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(files) == 0L) cli::cli_abort("No .arrow files found in {path}")
  nums <- as.integer(stringr::str_extract(basename(files), "\\d+"))
  files[order(nums)]
}

# Per-chunk work of .scan_vcf_gt(): per-sample totals and missing counts
# (missing = either allele NA, consistent with .pop_counts_from_chunk()
# elsewhere in this package, which treats a1/a2 symmetrically), plus an
# optional random subsample of (sample, DP) pairs.

.scan_chunk <- function(task, shared) {
  cols <- if (shared$collect_dp) c("a1", "a2", "DP") else c("a1", "a2")
  ch <- .read_live_chunk(task$fpath, cols, shared)
  cnt <- sample_counts_cpp(ch$s, ch$a1, ch$a2, numeric(0), length(shared$samples))
  out <- list(total = cnt$total, missing = cnt$total - cnt$called, dp = NULL)

  if (shared$collect_dp) {
    dp_ok <- !is.na(ch$DP)
    if (any(dp_ok)) {
      sub <- tibble::tibble(sample = shared$samples[ch$s[dp_ok]], DP = ch$DP[dp_ok])
      if (shared$keep_p < 1) {
        u <- .with_seed(task$seed, stats::runif(nrow(sub)))
        sub <- sub[u < shared$keep_p, , drop = FALSE]
      }
      if (nrow(sub) > 0L) out$dp <- sub
    }
  }
  out
}

# Evaluate expr with the RNG seeded by `seed`, then restore the caller's RNG
# state (so a chunk's own seed never disturbs the session's random stream).

.with_seed <- function(seed, expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(
    if (had_seed) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  )
  set.seed(seed)
  expr
}

# ── SNV blocks (vcf_filter_oneSNV(), vcf_filter_multiSNV()) ────────────────────
#
# Variants sorted by CHROM and POS, each assigned to a block of `block_size`
# bases counted from the first position of its CHROM:
#   block = (POS - min(POS of the CHROM)) %/% block_size + 1
# Done with vectorized operations on the sorted variants rather than
# dplyr::group_by(CHROM): CHROM can have nearly one value per variant (e.g.
# DiscoSnp path names), and per-group work over millions of groups was slow.
#
# The sort is stable and in the C locale (as dplyr::arrange()).  Variants
# without a position (POS NA, e.g. not a number or beyond R's integer range)
# cannot be placed in a block: they are left out, so they are removed by the
# calling filters, and do not affect the blocks of the other variants.
#
# Returns list($row_id: .row_id in sorted order, $run: block id (runs of
# equal CHROM and block, numbered 1, 2, ... in sorted order), $first: TRUE for
# the first variant of each block).

.snv_blocks <- function(variants, block_size) {
  no_pos <- is.na(variants$POS)
  if (any(no_pos)) {
    cli::cli_alert_warning(
      "{sum(no_pos)} variant{?s} without a valid POS removed (cannot be placed in a block)."
    )
    variants <- variants[!no_pos, , drop = FALSE]
  }
  o <- order(variants$CHROM, variants$POS, method = "radix")
  chrom <- variants$CHROM[o]
  pos <- variants$POS[o]
  n <- length(o)
  if (n == 0L)
    return(list(row_id = variants$.row_id[o], run = integer(0), first = logical(0)))

  # NA-safe "differs from the previous element"
  differs <- function(x) {
    if (length(x) < 2L) return(logical(0))
    a <- x[-1L]; b <- x[-length(x)]
    ne <- a != b
    ne[is.na(ne)] <- xor(is.na(a), is.na(b))[is.na(ne)]
    ne
  }

  # CHROM runs; the minimum POS of each is its first (sorted) POS
  chrom_start <- c(TRUE, differs(chrom))
  chrom_id <- cumsum(chrom_start)
  min_pos <- pos[chrom_start]

  block <- (pos - min_pos[chrom_id]) %/% block_size + 1

  first <- chrom_start | c(TRUE, differs(block))
  list(row_id = variants$.row_id[o], run = cumsum(first), first = first)
}

# ── Row-id lookup tables ───────────────────────────────────────────────────────
#
# Chunk loops must not look up .row_id with %in%, match() or a named vector:
# each call re-hashes the whole table (all variants), so the per-chunk cost
# grows with dataset size and total runtime grows quadratically.  Instead,
# build an integer vector indexed directly by .row_id once, then subset it.

# pos[id] = position of id in row_ids (first occurrence, as match()), 0 if absent.
.row_id_pos <- function(row_ids) {
  ok <- !is.na(row_ids) & !duplicated(row_ids)
  pos <- integer(if (any(ok)) max(row_ids[ok]) else 0L)
  pos[row_ids[ok]] <- which(ok)
  pos
}

# Drop-in replacement for match(ids, row_ids), given pos = .row_id_pos(row_ids).
.match_row_id <- function(ids, pos) {
  if (is.integer(ids)) return(match_row_id_cpp(ids, pos))
  p <- pos[ids]
  p[!is.na(p) & p == 0L] <- NA_integer_
  p
}

# ── Arrow memory release ───────────────────────────────────────────────────────
#
# Loops that keep chunks as Arrow Tables (rather than converting them to R)
# allocate in Arrow's memory pool, which R's garbage collector does not see:
# dropped Tables are only freed when R happens to collect, so they pile up.
# Collect once the pool exceeds `limit` bytes, instead of after every chunk.

.release_arrow_memory <- function(limit = 512 * 1024^2) {
  if (arrow::default_memory_pool()$bytes_allocated > limit)
    gc(verbose = FALSE, full = FALSE)
  invisible(NULL)
}

# ── Chunk reading and per-variant accumulation ─────────────────────────────────
#
# Position of each row's sample in `samples` (NA if absent), as match().
# For an Arrow Table this is computed in Arrow, which avoids converting the
# chunk's sample strings (one per genotype) to R strings — the largest
# remaining per-chunk cost.

.sample_index <- function(chunk, samples) {
  if (inherits(chunk, "ArrowTabular")) {
    idx <- arrow::call_function(
      "index_in", chunk$sample,
      options = list(value_set = arrow::Array$create(samples), skip_nulls = FALSE)
    )
    as.vector(idx) + 1L
  } else {
    match(chunk$sample, samples)
  }
}
#
# Read one feather chunk and keep only rows of live variants (idx$pos) and
# retained samples (idx$samples).  Returns a list of the requested columns as
# plain vectors plus
#   $pos  position of each row's variant in @variants (1..n_var)
#   $s    position of each row's sample in idx$samples
# Subsetting plain vectors is much cheaper than subsetting the tibble, and
# sample names are matched once per chunk.

.read_live_chunk <- function(fpath, cols, idx) {
  chunk <- arrow::read_feather(fpath, col_select = unique(c(".row_id", "sample", cols)),
                               as_data_frame = FALSE)
  pos <- .match_row_id(as.vector(chunk$.row_id), idx$pos)
  s <- .sample_index(chunk, idx$samples)
  if (anyNA(pos) || anyNA(s)) {
    keep <- !is.na(pos) & !is.na(s)
    out <- lapply(cols, function(col) as.vector(chunk[[col]])[keep])
    pos <- pos[keep]
    s <- s[keep]
  } else {  # every row kept (the common case): no subsetting copies
    out <- lapply(cols, function(col) as.vector(chunk[[col]]))
  }
  names(out) <- cols
  out$pos <- pos
  out$s <- s
  out
}

# Counts (.range_count) or sums of w (.range_sum) per position, covering only
# the range of positions present in this chunk.  Returns list($at, $value);
# accumulate with  acc[r$at] <- acc[r$at] + r$value  in the caller's own frame,
# which updates acc in place.  Tabulating over all n_var positions instead
# would make every chunk cost O(n_var).  Returns NULL for empty input.

.range_count <- function(pos) {
  if (length(pos) == 0L) return(NULL)
  lo <- min(pos); n <- max(pos) - lo + 1L
  list(at = lo - 1L + seq_len(n), value = tabulate(pos - lo + 1L, n))
}

.range_sum <- function(pos, w) {
  if (length(pos) == 0L) return(NULL)
  if (is.unsorted(pos)) {
    o <- order(pos); pos <- pos[o]; w <- w[o]
  }
  lo <- pos[1L]; n <- pos[length(pos)] - lo + 1L
  last <- c(pos[-1L] != pos[-length(pos)], TRUE)  # last row of each position
  sums <- diff(c(0, cumsum(as.numeric(w))[last]))
  value <- numeric(n)
  value[pos[last] - lo + 1L] <- sums
  list(at = lo - 1L + seq_len(n), value = value)
}

# Sum per-chunk range results (list(at, value), see above) stored under `key`
# in each element of `parts` into `acc`.  The additions happen in this
# function's own frame, so acc is updated in place.

.merge_ranges <- function(parts, key, acc) {
  for (p in parts) {
    r <- p[[key]]
    if (!is.null(r)) acc[r$at] <- acc[r$at] + r$value
  }
  acc
}

# Build the two fast-lookup structures used by every filter's chunk loop.
# Returns a list with $lv (logical vector for row_id membership) and
# $pos (integer vector: row_id → 1-based position, see .row_id_pos()).

.vcf_filter_index <- function(vcf_arrow) {
  valid <- vcf_arrow@variants$.row_id
  max_id <- max(valid)
  lv <- logical(max_id + 1L)
  lv[valid] <- TRUE
  list(
    lv = lv,
    pos = .row_id_pos(valid),
    n_var = length(valid),
    samples = vcf_arrow@samples
  )
}

# ── Shared chunk-reading helper ────────────────────────────────────────────────
#
# Reads every feather file directly, filtering each chunk to the VCFArrow's
# CURRENT variant set (vcf_arrow@variants$.row_id) and CURRENT sample list
# (vcf_arrow@samples) — the same filtering every exporter in this package
# applies, and the same fix that was needed in write_vcf() to avoid emitting
# already-filtered-out variants.
#
# Accumulates missingness counts in O(n_samples) memory.  Optionally also
# collects a bounded random subsample of (sample, DP) pairs for the coverage
# violin plot, controlled by max_points_per_sample so memory stays bounded
# regardless of dataset size.

.scan_vcf_gt <- function(vcf_arrow, collect_dp = FALSE,
                         max_points_per_sample = 5000L) {

  samples <- vcf_arrow@samples
  valid_row_ids <- vcf_arrow@variants$.row_id
  n_var <- length(valid_row_ids)
  n_samples <- length(samples)

  feather_files <- list.files(vcf_arrow@path, pattern = "\\.arrow$",
                              full.names = TRUE)
  if (length(feather_files) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")
  chunk_nums <- as.integer(stringr::str_extract(basename(feather_files), "\\d+"))
  feather_files <- feather_files[order(chunk_nums)]

  total_loci <- stats::setNames(integer(n_samples), samples)
  missing_n <- stats::setNames(integer(n_samples), samples)

  # Probability of keeping a row for the DP subsample: each sample has
  # exactly n_var rows total in the dense gt_long format, so a flat
  # per-row keep-probability of max_points_per_sample / n_var yields
  # ~max_points_per_sample retained rows per sample, independent of which
  # chunk they fall in.
  keep_p <- if (collect_dp) min(1, max_points_per_sample / max(n_var, 1L)) else 0

  # The DP subsample of each chunk uses its own seed, drawn here from the
  # session's RNG, so results are reproducible with set.seed() and identical
  # whether chunks run serially or on workers.
  shared <- list(pos = .row_id_pos(valid_row_ids), samples = samples,
                 collect_dp = collect_dp, keep_p = keep_p)
  tasks <- lapply(seq_along(feather_files), function(i) list(
    fpath = feather_files[[i]],
    seed = if (collect_dp && keep_p < 1) sample.int(.Machine$integer.max, 1L) else NA_integer_
  ))

  parts <- .map_chunks(tasks, .scan_chunk, shared = shared, label = "Scanning chunk")
  for (p in parts) {
    total_loci <- total_loci + p$total
    missing_n <- missing_n + p$missing
  }

  list(
    total_loci = total_loci,
    missing_n = missing_n,
    dp_df = if (collect_dp) do.call(rbind, Filter(Negate(is.null), lapply(parts, `[[`, "dp")))
    else NULL
  )
}
