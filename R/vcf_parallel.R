#' @title vcf_set_workers
#'
#' @description
#' Set the number of parallel worker processes used by VCFArrow
#'
#' @author Tomas Hrbek September 2026
#'
#' @param workers -> number of worker processes, default 1 (no parallelism) (integer)
#'
#' @return Invisibly, the previous number of workers.
#'
#' @details
#' Functions that process the on-disk genotype chunks independently
#' distribute the chunks over `workers` background R processes: read_vcf(),
#' the genotype-based filters, sample and group extraction, vcf_stats() and
#' the assess functions, vcf_bind(), write_vcf(), the vcf2*() exporters and
#' vcf2gt_long() with CSV output (with feather or Parquet output it runs
#' serially). Results do not depend on the number of workers.  The workers are started once and reused
#' until `vcf_set_workers(1)` is called or the session ends.
#' Each worker holds one chunk in memory at a time, so memory use grows with
#' the number of workers (see vcf_memory_estimate()).  Useful values are up
#' to the number of physical CPU cores; disk speed limits the gain beyond that.
#'
#' @examples
#' \dontrun{
#' vcf_set_workers(4)
#' vcf <- read_vcf(f)
#' vcf_set_workers(1)
#' }
#'
#' @export
#'

vcf_set_workers <- function(workers = 1L) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L)
    cli::cli_abort("{.arg workers} must be a positive integer.")

  old <- .vcf_workers()
  .stop_cluster()
  if (workers > 1L) {
    cl <- parallel::makePSOCKcluster(workers)
    # workers load VCFArrow from the library this session loaded it from
    # (first on their library path, so they run the same installed version)
    # and keep Arrow to one thread each, so n workers do not start
    # n x n_cpu Arrow threads
    lib <- c(dirname(getNamespaceInfo("VCFArrow", "path")), .libPaths())
    init <- .init_worker
    # sent without the VCFArrow namespace as its environment: unserializing a
    # namespace function would load VCFArrow before the library path is set
    environment(init) <- baseenv()
    tryCatch(parallel::clusterCall(cl, init, unique(lib)), error = function(e) {
      parallel::stopCluster(cl)
      cli::cli_abort(c("Could not start VCFArrow worker processes.", "x" = conditionMessage(e)))
    })
    .vcfarrow_parallel$cl <- cl
    cli::cli_alert_info("VCFArrow will use {workers} worker process{?es}.")
  }
  invisible(old)
}

# Package-level parallel state: the worker cluster (NULL = serial) and, in a
# worker process, the data shared by the chunks of the current call.
.vcfarrow_parallel <- new.env(parent = emptyenv())

.vcf_workers <- function() {
  cl <- .vcfarrow_parallel$cl
  if (is.null(cl)) 1L else length(cl)
}

.stop_cluster <- function() {
  cl <- .vcfarrow_parallel$cl
  .vcfarrow_parallel$cl <- NULL
  if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
  invisible(NULL)
}

.init_worker <- function(lib_paths) {
  .libPaths(lib_paths)
  if ("VCFArrow" %in% loadedNamespaces())
    stop("VCFArrow was loaded in the worker before its library path was set")
  loadNamespace("VCFArrow")
  arrow::set_cpu_count(1L)
  arrow::set_io_thread_count(2L)
  invisible(NULL)
}

.onUnload <- function(libpath) .stop_cluster()

.set_worker_shared <- function(shared) {
  .vcfarrow_parallel$shared <- shared
  invisible(NULL)
}

.run_worker_task <- function(x, FUN) FUN(x, .vcfarrow_parallel$shared)


# ── Apply FUN to every chunk ──────────────────────────────────────────────────
#
# .chunk_iterator(xs, FUN, shared, label) runs FUN(x, shared) over xs in
# waves — one x per worker with workers (see vcf_set_workers()), one x at a
# time otherwise — and returns an iterator: it$next_wave() gives the next
# wave's results as list(idx, res) (NULL when done), it$done() cleans up.
# Consuming results wave by wave lets callers merge them into accumulators in
# their own frame (updated in place) without holding all results at once.
# `shared` (e.g. row-id lookup tables) is sent to each worker once;
# local = TRUE runs everything in this process (e.g. when FUN returns Arrow
# objects, which cannot be sent back from workers).  The
# progress bar (if `label` is given) belongs to .envir, the calling
# function's frame, so cli closes it when that function exits.
#
# FUN must be a function defined at the top level of the package namespace,
# and xs / shared must hold only what FUN needs: a closure defined inside
# another function would carry (and serialize) its whole enclosing
# environment, such as the VCFArrow object.

.chunk_iterator <- function(xs, FUN, shared = list(), label = NULL,
                            local = FALSE, .envir = parent.frame()) {
  it <- new.env(parent = emptyenv())
  n <- length(xs)
  cl <- .vcfarrow_parallel$cl
  parallel <- !local && !is.null(cl) && n > 1L
  wave <- if (parallel) length(cl) else 1L
  pos <- 0L
  if (parallel) parallel::clusterCall(cl, .set_worker_shared, shared)
  if (!is.null(label)) cli::cli_progress_bar(label, total = n, .envir = .envir)

  it$next_wave <- function() {
    if (pos >= n) return(NULL)
    idx <- (pos + 1L):min(pos + wave, n)
    res <- if (parallel) {
      parallel::clusterApply(cl, xs[idx], .run_worker_task, FUN = FUN)
    } else {
      lapply(xs[idx], FUN, shared)
    }
    pos <<- max(idx)
    if (!is.null(label)) cli::cli_progress_update(inc = length(idx), .envir = .envir)
    list(idx = idx, res = res)
  }
  it$done <- function() {
    if (!is.null(label)) cli::cli_progress_done(.envir = .envir)
    if (parallel) try(parallel::clusterCall(cl, .set_worker_shared, NULL), silent = TRUE)
    invisible(NULL)
  }
  it
}

# All results of FUN over xs, in order (see .chunk_iterator()).
.map_chunks <- function(xs, FUN, shared = list(), label = "Processing chunk",
                        progress = TRUE) {
  out <- vector("list", length(xs))
  it <- .chunk_iterator(xs, FUN, shared, if (progress) label)
  on.exit(it$done())
  while (!is.null(w <- it$next_wave())) out[w$idx] <- w$res
  out
}

# One task per chunk file, for .write_chunks_ordered()
.file_tasks <- function(files) lapply(files, function(f) list(fpath = f))

# Write chunks to one output file, in order.  FUN(task, shared) appends a
# chunk to task$part: serially task$part is out_file itself; with workers each
# chunk goes to its own temporary part file, and the parts are appended to
# out_file in order after each wave.  tasks: list of lists, to which $part is
# added.

.write_chunks_ordered <- function(tasks, FUN, shared, out_file,
                                  label = "Writing chunk") {
  parallel <- .vcf_workers() > 1L && length(tasks) > 1L
  part_dir <- NULL
  if (parallel) {
    part_dir <- tempfile("vcfarrow_parts_")
    dir.create(part_dir)
    on.exit(unlink(part_dir, recursive = TRUE), add = TRUE)
  }
  for (i in seq_along(tasks))
    tasks[[i]]$part <- if (parallel) file.path(part_dir, paste0("part_", i)) else out_file

  it <- .chunk_iterator(tasks, FUN, shared, label)
  on.exit(it$done(), add = TRUE)
  while (!is.null(w <- it$next_wave())) {
    if (!parallel) next
    for (i in w$idx) {
      part <- tasks[[i]]$part
      if (file.exists(part)) {
        if (!file.append(out_file, part))
          cli::cli_abort("Could not append to {.file {out_file}}.")
        unlink(part)
      }
    }
  }
  invisible(out_file)
}
