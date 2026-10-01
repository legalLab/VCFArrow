#' @title vcf_filter_adr
#'
#' @description
#' Correct or remove genotypes with > % ALT/REF ratio from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param mode -> switch between 'correct' or 'remove' mode (character)
#' @param threshold -> decimal missing threshold, default 0.1 (numeric)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function either changes genotypes or makes genotypes missing in
#' a VCFArrow object if they have above/below threshold normalized ADR ratio,
#' returning a new VCFArrow object.
#' The ADR is calculated as ADR = ALT / (REF + ALT).
#' If ADR > threshold, the genotype becomes ALT homozygous.
#' If ADR < threshold, the genotype becomes REF homozygous.
#' By default will remove any loci that may have become invariant as the
#' result of the change/removal of genotypes
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_adr(vcf, mode = "correct")
#'
#' @export
#'

vcf_filter_adr <- function(vcf_arrow, mode = c("correct", "remove"),
                           threshold = 0.1, f_invar = TRUE) {

  mode <- match.arg(mode)

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  if (threshold <= 0 || threshold >= 0.5)
    cli::cli_abort("{.arg threshold} must be strictly between 0 and 0.5")

  ffiles <- .get_sorted_feather_files(vcf_arrow@path)
  tmp_dir <- tempfile("arrow_vcf_adr_")
  dir.create(tmp_dir)

  cli::cli_alert_info("Applying ADR filter")

  cli::cli_alert_info(
    "Applying ADR {mode} (threshold = {threshold}) across {length(ffiles)} chunk{?s}..."
  )
  tasks <- lapply(seq_along(ffiles), function(i) list(
    fpath = ffiles[[i]], out = file.path(tmp_dir, paste0("chunk_", i, ".arrow"))
  ))
  .map_chunks(tasks, .adr_chunk, shared = list(mode = mode, threshold = threshold),
              label = "Rewriting chunk")

  gt_arrow <- suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather"))

  new_vcfarrow <- .new_vcfarrow(
    header = vcf_arrow@header,
    info = vcf_arrow@info,
    format = vcf_arrow@format,
    variants = vcf_arrow@variants,
    gt = gt_arrow,
    samples = vcf_arrow@samples,
    groups = vcf_arrow@groups,
    path = tmp_dir,
    invariant_removed = vcf_arrow@invariant_removed,
    loci = vcf_arrow@loci
  )

  if (f_invar) new_vcfarrow <- vcf_filter_invariant(new_vcfarrow)

  return(new_vcfarrow)
}

# Apply the ADR correction/removal to one chunk and write it to task$out.
# The chunk stays an Arrow Table and only ADR, a1 and a2 are converted to R:
# converting the other columns (notably the per-sample FORMAT strings) to R
# and back is the dominant cost.

.adr_chunk <- function(task, shared) {
  threshold <- shared$threshold
  chunk <- arrow::read_feather(task$fpath, as_data_frame = FALSE)
  adr <- as.vector(chunk$ADR)
  a1 <- as.vector(chunk$a1)
  a2 <- as.vector(chunk$a2)

  # adr_flag: TRUE only where ADR is actually observed AND outside bounds.
  # NA ADR → FALSE (genotype passes through untouched) — this is the fix.
  adr_flag <- !is.na(adr) &
    (adr < threshold | adr > (1 - threshold))
  adr_dir  <- ifelse(adr <= threshold, 0L,
                     ifelse(adr >= (1 - threshold), 1L, NA_integer_))

  if (shared$mode == "correct") {
    a1[adr_flag] <- adr_dir[adr_flag]
    a2[adr_flag] <- adr_dir[adr_flag]
  } else {  # mode == "remove"
    a1[adr_flag] <- NA_integer_
    a2[adr_flag] <- NA_integer_
  }
  chunk[["a1"]] <- a1
  chunk[["a2"]] <- a2

  # Keep the per-sample FORMAT strings (written by write_vcf()) consistent:
  # replace the GT field of every changed genotype, keeping its separator.
  # Only those strings are converted to R.
  if (any(adr_flag)) {
    fmt <- chunk$fmt
    fmt <- if (fmt$num_chunks == 1L) fmt$chunk(0L) else do.call(arrow::concat_arrays, fmt$chunks)
    old <- as.vector(fmt$Take(arrow::Array$create(which(adr_flag) - 1L)))
    sep <- ifelse(substr(old, 2L, 2L) == "|", "|", "/")
    allele <- ifelse(is.na(a1[adr_flag]), ".", as.character(a1[adr_flag]))
    new <- paste0(allele, sep, allele, sub("^[^:]*", "", old))
    chunk[["fmt"]] <- arrow::call_function(
      "replace_with_mask", fmt, arrow::Array$create(adr_flag), arrow::Array$create(new)
    )
  }

  arrow::write_feather(chunk, task$out)
  chunk <- NULL
  .release_arrow_memory()
  invisible(NULL)
}
