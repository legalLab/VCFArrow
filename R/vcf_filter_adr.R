#' @title vcf_filter_adr
#'
#' @description
#' Correct or remove genotypes with allele imbalance greater than a threshold 
#' from a VCFArrow object
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
#' This function either changes heterozygous genotypes or makes them missing
#' in a VCFArrow object if they have above/below threshold normalized ADR
#' ratio, returning a new VCFArrow object.
#' The ADR is calculated as ADR = ALT / (REF + ALT), from the AD field.
#' Only called heterozygous genotypes are affected; homozygous and missing
#' genotypes are left unchanged, even where AD gives them an ADR.
#' In 'correct' mode, a heterozygote with ADR < threshold becomes REF
#' homozygous, and one with 1 - ADR < threshold becomes ALT homozygous.
#' In 'remove' mode, these heterozygotes become missing.
#' The number of genotypes corrected or removed is reported.
#' By default will remove any loci that may have become invariant as the
#' result of the change/removal of genotypes.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_adr(vcf, mode = "correct")
#'
#' @export
#'

vcf_filter_adr <- function(vcf_arrow, threshold = 0.1, 
                           mode = c("correct", "remove"), f_invar = TRUE) {

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
  counts <- .map_chunks(tasks, .adr_chunk, shared = list(mode = mode, threshold = threshold),
                        label = "Rewriting chunk")
  counts <- Reduce(`+`, counts)
  if (mode == "correct") {
    cli::cli_alert_info(
      "Corrected {counts[['changed']]} / {counts[['called']]} genotype{?s}"
    )
  } else {
    cli::cli_alert_info(
      "Removed {counts[['changed']]} / {counts[['called']]} genotype{?s} (set to missing)"
    )
  }

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
# Returns the number of genotypes changed (a1 or a2 differs from before) and
# of genotypes called before the change.
# The chunk stays an Arrow Table and only ADR, a1 and a2 are converted to R:
# converting the other columns (notably the per-sample FORMAT strings) to R
# and back is the dominant cost.

.adr_chunk <- function(task, shared) {
  threshold <- shared$threshold
  chunk <- arrow::read_feather(task$fpath, as_data_frame = FALSE)
  adr <- as.vector(chunk$ADR)
  a1 <- as.vector(chunk$a1)
  a2 <- as.vector(chunk$a2)
  a1_old <- a1
  a2_old <- a2

  # adr_flag: TRUE only for called heterozygotes whose ADR is observed and
  # outside bounds. Missing genotypes can have an ADR (AD is written for
  # uncalled genotypes too) and homozygotes have ADR near 0 or 1, so both are
  # excluded; NA ADR → FALSE (genotype passes through untouched).
  adr_flag <- !is.na(a1) & !is.na(a2) & a1 != a2 & !is.na(adr) &
    (adr < threshold | adr > (1 - threshold))
  adr_dir <- ifelse(adr <= threshold, 0L,
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
  differs <- function(x, y) is.na(x) != is.na(y) | (!is.na(x) & !is.na(y) & x != y)
  counts <- c(changed = sum(differs(a1_old, a1) | differs(a2_old, a2)),
              called = sum(!is.na(a1_old) & !is.na(a2_old)))

  chunk <- .set_fmt_gt(chunk, adr_flag, a1)

  arrow::write_feather(chunk, task$out)
  chunk <- NULL
  .release_arrow_memory()
  counts
}
