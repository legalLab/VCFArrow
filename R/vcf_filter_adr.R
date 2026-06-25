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
#' vcf_filter_adr(vcf_arrow = my_vcf, mode = "correct", threshold = my_threshold, f_invar = TRUE)
#' vcf_filter_adr(my_vcf, "correct", my_threshold, TRUE)
#' vcf_filter_adr(my_vcf, "correct")
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

  cli::cli_alert_info(
    "Applying ADR {mode} (threshold = {threshold}) across {length(ffiles)} chunk(s)..."
  )
  cli::cli_progress_bar("Rewriting chunk", total = length(ffiles))

  for (i in seq_along(ffiles)) {
    chunk <- arrow::read_feather(ffiles[[i]])   # full columns, full rows

    # adr_flag: TRUE only where ADR is actually observed AND outside bounds.
    # NA ADR → FALSE (genotype passes through untouched) — this is the fix.
    adr_flag <- !is.na(chunk$ADR) &
      (chunk$ADR < threshold | chunk$ADR > (1 - threshold))
    adr_dir  <- ifelse(chunk$ADR <= threshold, 0L,
                       ifelse(chunk$ADR >= (1 - threshold), 1L, NA_integer_))

    if (mode == "correct") {
      chunk$a1[adr_flag] <- adr_dir[adr_flag]
      chunk$a2[adr_flag] <- adr_dir[adr_flag]
    } else {   # mode == "remove"
      chunk$a1[adr_flag] <- NA_integer_
      chunk$a2[adr_flag] <- NA_integer_
    }

    arrow::write_feather(chunk, file.path(tmp_dir, paste0("chunk_", i, ".arrow")))
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  gt_arrow <- suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather"))

  new_vcfarrow <- .new_vcfarrow(
    vcf_arrow@header,
    vcf_arrow@info,
    vcf_arrow@format,
    vcf_arrow@variants,
    gt_arrow,
    vcf_arrow@samples,
    vcf_arrow@groups,
    tmp_dir
  )

  if (f_invar) new_vcfarrow <- vcf_filter_invariant(new_vcfarrow)

  return(new_vcfarrow)
}
