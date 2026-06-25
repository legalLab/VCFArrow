#' @noRd

# vcfarrow_lifecycle.R
#
# VCFArrow object lifecycle: registration, finalizers, temp-directory cleanup.
#
# Architecture
# ─────────────────────────────────────────────────────────────────────────────
# Two registries live in the package namespace:
#
#   .vcfarrow_registry  – paths of LIVE objects (ref-count > 0)
#   .vcfarrow_pending   – paths whose ref-count hit 0 but unlink() has not
#                         yet succeeded (Arrow file handles may still be open)
#
# Lifecycle of a temp directory
#   read_vcf / vcf_bind_sparse  →  .register_vcfarrow()   [count = 1]
#   each filter that inherits the same path  →  .register_vcfarrow()  [count + 1]
#   GC of a VCFArrow object fires finalizer  →  .deregister_vcfarrow()
#     if count → 0: move path to .vcfarrow_pending, attempt immediate unlink
#   vcf_gc() or session exit  →  retry unlink for everything in .vcfarrow_pending
#
# Why finalizer_env is a SLOT (not attr())
#   S4 slots of type "environment" have reference semantics: they are not
#   copied when the S4 object is copied.  The reference count on the
#   environment therefore precisely tracks all live VCFArrow copies.
#   attr() on S4 objects has undefined copy behaviour and must not be used.
#
# Why .onLoad has ONE definition
#   R silently keeps only the last .onLoad in a file.  All initialisation
#   (registry setup + stale-dir cleanup) must be in a single function.

# ── Arrow allocator: switch from jemalloc to system malloc ────────────────────
#
# Must be called BEFORE any Arrow objects are created in the session.
# The best place is .onLoad(); doing it lazily (on first use) is too late
# because Arrow initialises its pool at package load time.
#
# Effect: freed Arrow C++ memory is returned to the OS immediately, so cgroup
# memory tracks R's actual working set rather than Arrow's held pool.
#
# Arrow exposes this via the ARROW_DEFAULT_MEMORY_POOL environment variable.
# Setting it before Arrow loads its C++ library is the reliable path; calling
# arrow::set_io_thread_count() or similar after the fact does not help.

.set_arrow_allocator <- function() {
  # Only override if the user has not set their own preference.
  if (Sys.getenv("ARROW_DEFAULT_MEMORY_POOL") == "") {
    Sys.setenv(ARROW_DEFAULT_MEMORY_POOL = "system")
  }
}


# ── Package-level registries (overwritten by .onLoad; defined here for R CMD check) ──
.vcfarrow_registry <- new.env(parent = emptyenv())
.vcfarrow_pending  <- new.env(parent = emptyenv())


# ── Single .onLoad ────────────────────────────────────────────────────────────

.onLoad <- function(libname, pkgname) {
  ns <- environment()   # the package namespace
  assign(".vcfarrow_registry", new.env(parent = emptyenv()), envir = ns)
  assign(".vcfarrow_pending",  new.env(parent = emptyenv()), envir = ns)

  # Remove stale directories left by a previous crashed/hard-killed session.
  # Normal shutdown fires onexit finalizers, but a hard crash skips them.
  .cleanup_old_vcfarrow()

  # switch from jemalloc to system malloc so system memory gets deallocated
  .set_arrow_allocator()

  # Register a session-exit hook so cleanup runs even if vcf_gc() is never called.
  reg.finalizer(ns, function(e) {
    .flush_pending(verbose = FALSE)
  }, onexit = TRUE)
}


# ── Low-level helpers ─────────────────────────────────────────────────────────

# Attempt to delete path; return TRUE on success.
.try_unlink <- function(path) {
  if (!dir.exists(path)) return(TRUE)  # already gone
  res <- tryCatch(
    unlink(path, recursive = TRUE, force = TRUE),
    error = function(e) 1L
  )
  res == 0L && !dir.exists(path)  # 0 = success from unlink()
}

# Process everything in .vcfarrow_pending.
.flush_pending <- function(verbose = TRUE) {
  paths <- ls(.vcfarrow_pending)
  n_ok  <- 0L
  for (path in paths) {
    if (.try_unlink(path)) {
      rm(list = path, envir = .vcfarrow_pending)
      n_ok <- n_ok + 1L
    }
  }
  if (verbose && n_ok > 0L)
    message("VCFArrow: cleaned up ", n_ok, " pending temp director",
            if (n_ok == 1L) "y." else "ies.")
  invisible(n_ok)
}


# ── Registry helpers ──────────────────────────────────────────────────────────

.register_vcfarrow <- function(path) {
  if (exists(path, envir = .vcfarrow_registry, inherits = FALSE)) {
    .vcfarrow_registry[[path]] <- .vcfarrow_registry[[path]] + 1L
  } else {
    .vcfarrow_registry[[path]] <- 1L
  }
}

.deregister_vcfarrow <- function(path) {
  if (!exists(path, envir = .vcfarrow_registry, inherits = FALSE)) return()

  .vcfarrow_registry[[path]] <- .vcfarrow_registry[[path]] - 1L

  if (.vcfarrow_registry[[path]] <= 0L) {
    rm(list = path, envir = .vcfarrow_registry)

    # Move to pending queue; attempt immediate deletion.
    # If Arrow still holds file handles, unlink() may fail silently on Linux
    # or return non-zero on Windows; .flush_pending() will retry later.
    .vcfarrow_pending[[path]] <- TRUE
    .try_unlink(path)
    if (!dir.exists(path) &&
        exists(path, envir = .vcfarrow_pending, inherits = FALSE)) {
      rm(list = path, envir = .vcfarrow_pending)
    }
  }
}


# ── Finalizer factory ─────────────────────────────────────────────────────────

# Returns a fresh lightweight environment with a finalizer that calls
# .deregister_vcfarrow(path) when collected.  Store the return value in the
# VCFArrow object's `finalizer_env` SLOT (not as an attr()).
.make_finalizer_env <- function(path) {
  ptr <- new.env(parent = emptyenv())
  reg.finalizer(ptr, function(e) .deregister_vcfarrow(path), onexit = TRUE)
  ptr
}


# ── Constructor ───────────────────────────────────────────────────────────────
#
# The VCFArrow S4 class MUST include:
#   finalizer_env = "environment"
# in its setClass() representation.

.new_vcfarrow <- function(header, info, format, variants, gt,
                          samples, groups, path) {
  .register_vcfarrow(path)
  new("VCFArrow",
      header = header,
      info = info,
      format = format,
      variants = variants,
      gt = gt,
      samples = samples,
      groups = groups,
      path = path,
      finalizer_env = .make_finalizer_env(path)
  )
}


# ── Startup cleanup ───────────────────────────────────────────────────────────

.cleanup_old_vcfarrow <- function() {
  tmp   <- tempdir()
  dirs  <- list.dirs(tmp, recursive = FALSE, full.names = TRUE)
  stale <- dirs[grepl("arrow_vcf_", basename(dirs))]
  for (d in stale) try(.try_unlink(d), silent = TRUE)
}


# ── User-facing cleanup ───────────────────────────────────────────────────────

#' Garbage-collect VCFArrow temp directories
#'
#' Triggers R's garbage collector (three full passes to handle the
#' finalizer → pending-queue → unlink chain), then flushes any directories
#' whose reference count has already reached zero.
#'
#' @param force Logical.  If TRUE, also force-deregister and delete directories
#'   for objects that are still nominally live in the registry.  Use this when
#'   you have called \code{rm()} on all VCFArrow objects but the temp
#'   directories have not been cleaned up (common in RStudio, which can hold
#'   display references that delay GC).
#' @param verbose Logical.  Print a status message.
#'
#' @section Typical workflow:
#' \preformatted{
#'   rm(vcf1, vcf2, vcf3)
#'   vcf_gc()              # usually sufficient
#'   vcf_gc(force = TRUE)  # if directories are still present after rm()
#' }
vcf_gc <- function(force = FALSE, verbose = TRUE) {

  # Three full GC passes:
  #   pass 1 – collect VCFArrow objects → fires finalizer_env finalizers
  #             → .deregister_vcfarrow() → moves paths to .vcfarrow_pending
  #   pass 2 – collect anything freed in pass 1 (including Arrow Dataset SEXPs
  #             which were referenced by VCFArrow@gt)
  #   pass 3 – collect anything freed in pass 2 (Arrow's internal C++ objects)
  for (i in seq_len(3L)) gc(verbose = FALSE, full = TRUE)

  # Process pending queue (paths whose ref-count already hit 0)
  .flush_pending(verbose = FALSE)

  # Optional: force-deregister live registry entries
  if (force) {
    live <- ls(.vcfarrow_registry)
    if (length(live) > 0L) {
      for (path in live) {
        .try_unlink(path)
        if (exists(path, envir = .vcfarrow_registry, inherits = FALSE))
          rm(list = path, envir = .vcfarrow_registry)
        if (exists(path, envir = .vcfarrow_pending, inherits = FALSE))
          rm(list = path, envir = .vcfarrow_pending)
      }
      if (verbose)
        message("VCFArrow: force-cleaned ", length(live), " live director",
                if (length(live) == 1L) "y." else "ies.")
    }
  }

  if (verbose) {
    n_live    <- length(ls(.vcfarrow_registry))
    n_pending <- length(ls(.vcfarrow_pending))

    if (n_live == 0L && n_pending == 0L) {
      message("VCFArrow: all temp directories cleaned up.")
    } else {
      if (n_live > 0L) {
        live_paths <- ls(.vcfarrow_registry)
        message(
          "VCFArrow: ", n_live, " temp director",
          if (n_live == 1L) "y is" else "ies are",
          " still live (VCFArrow objects not yet GC'd):\n",
          paste0("  ", live_paths, collapse = "\n"), "\n",
          "  \u2192 rm() all VCFArrow objects, then call vcf_gc() again.\n",
          "  \u2192 Or call vcf_gc(force = TRUE) to delete regardless."
        )
      }
      if (n_pending > 0L) {
        message(
          "VCFArrow: ", n_pending, " director",
          if (n_pending == 1L) "y" else "ies",
          " pending (Arrow file handles may still be open).\n",
          "  \u2192 Call vcf_gc() again in a moment."
        )
      }
    }
  }

  invisible(NULL)
}


# ── Memory estimation ─────────────────────────────────────────────────────────
#
# Call before running an export to check whether enough RAM is available.
# Reports both the accumulation-matrix footprint and the Arrow pool overhead.

#' Estimate RAM required for a VCFArrow export operation
#'
#' @param vcf_arrow A VCFArrow object.
#' @param keep_groups Groups to export (NULL = all).
#' @param format     One of "individual" (Structure, Arlequin, FASTA, …) or
#'                   "pop" (BayesScan, Treemix, Migrate-N C, …) or
#'                   "chunk" (SmartSNP, fineRADstructure, sNMF, EIGENSTRAT, …).
#' @param chunk_size Feather chunk size used at read_vcf() time.
#' @param lowmem     If TRUE, estimate uses raw-byte matrices (vcf2*() lowmem
#'                   variants); otherwise integer matrices.
vcf_memory_estimate <- function(vcf_arrow,
                                keep_groups = NULL,
                                format      = c("individual", "pop", "chunk"),
                                chunk_size  = 100000L,
                                lowmem      = FALSE) {
  format <- match.arg(format)

  all_samples <- vcf_arrow@samples
  all_groups  <- vcf_arrow@groups
  if (is.null(keep_groups)) keep_groups <- unique(all_groups)
  n_samples <- sum(all_groups %in% keep_groups)

  n_var <- vcf_arrow@variants |>
    dplyr::filter(is_biallelic, !is_indel) |>
    nrow()

  n_pops       <- length(keep_groups)
  bytes_per_el <- if (lowmem) 1L else 4L

  mb <- function(x) paste0(round(x / 1024^2, 1), " MiB")

  # Arrow C++ pool: one chunk fully materialised per read_feather() call.
  # ~5 columns (row_id, sample, a1, a2, phased) × 4 bytes each.
  chunk_arrow_mb <- chunk_size * n_samples * 5L * 4L
  # Accumulation matrix (a1 + a2)
  mat_bytes <- 2L * n_samples * n_var * bytes_per_el   # individual
  pop_bytes <- 2L * n_pops   * n_var * 4L              # pop-level (always int)

  total <- switch(format,
                  individual = mat_bytes + chunk_arrow_mb,
                  pop        = pop_bytes + chunk_arrow_mb,
                  chunk      = chunk_arrow_mb * 2L   # double-buffer for reshape
  )

  cli::cli_h2("VCFArrow memory estimate")
  cli::cli_ul(c(
    "Variants (filtered): {n_var}",
    "Samples retained:    {n_samples}",
    "Populations:         {n_pops}",
    "Chunk size:          {format(chunk_size, big.mark=',')} variants",
    "Low-memory mode:     {lowmem}"
  ))
  cli::cli_h3("Per-operation footprint")
  if (format %in% c("individual", "chunk"))
    cli::cli_bullets(c(
      "*" = "Arrow pool per chunk read: {mb(chunk_arrow_mb)}",
      "*" = "Accumulation matrices (a1+a2): {mb(mat_bytes)}  [{bytes_per_el} B/cell]",
      "*" = "Estimated peak RAM: {mb(total)}"
    ))
  else
    cli::cli_bullets(c(
      "*" = "Arrow pool per chunk read: {mb(chunk_arrow_mb)}",
      "*" = "Pop-count matrices (ref+alt+nobs): {mb(pop_bytes)}",
      "*" = "Estimated peak RAM: {mb(total)}"
    ))

  invisible(list(
    n_var             = n_var,
    n_samples         = n_samples,
    n_pops            = n_pops,
    chunk_arrow_bytes = chunk_arrow_mb,
    matrix_bytes      = if (format == "individual") mat_bytes else pop_bytes,
    peak_bytes        = total
  ))
}


# ── Recommended chunk_size for available RAM ──────────────────────────────────

#' Suggest a chunk_size for read_vcf() given available RAM
#'
#' @param available_gb RAM available in gigabytes.
#' @param n_samples    Number of samples.
#' @param n_columns    Number of columns per gt row (default 5: row_id, sample,
#'                     a1, a2, phased).
vcf_suggest_chunk_size <- function(available_gb, n_samples, n_columns = 5L) {
  # Leave half for OS + R overhead + accumulation matrices
  arrow_budget <- available_gb * 1024^3 / 2
  bytes_per_variant <- n_samples * n_columns * 4L
  chunk <- floor(arrow_budget / bytes_per_variant)
  chunk <- max(1000L, min(chunk, 500000L))
  cli::cli_alert_info(
    "For {available_gb} GB RAM and {n_samples} samples, \\
     recommended chunk_size = {format(chunk, big.mark=',')}
     (use: read_vcf(vcf_file, chunk_size = {format(chunk, big.mark=',')}))"
  )
  invisible(chunk)
}
