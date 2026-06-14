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
