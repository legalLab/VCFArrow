#' @title vcf_filter_missing
#'
#' @description
#' Remove samples with > % missing data from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal missing threshold, default 0.5 (numeric)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report filtering stats, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes samples from a VCFArrow object if they have
#' above threshold missing loci, returning a new VCFArrow object.
#' By default will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' By default will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_filter_missing(vcf_arrow = my_vcf, threshold = my_threshold)
#' vcf_filter_missing(my_vcf, my_threshold)
#' vcf_filter_missing(my_vcf)
#'
#' @export
#'

vcf_filter_missing <- function(vcf_arrow, threshold = 0.5,
                               f_invar = TRUE, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  samples <- idx$samples
  miss_n <- stats::setNames(integer(length(samples)), samples)
  total_n <- stats::setNames(integer(length(samples)), samples)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying sample missingness filter")

  cli::cli_progress_bar("Scanning chunk", total = length(ffiles))
  for (fpath in ffiles) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1"))
    chunk <- chunk[idx$lv[chunk$.row_id] & chunk$sample %in% samples, , drop = FALSE]
    if (nrow(chunk) > 0L) {
      tt <- table(chunk$sample)
      total_n[names(tt)] <- total_n[names(tt)] + as.integer(tt)
      na_samp <- chunk$sample[is.na(chunk$a1)]
      if (length(na_samp) > 0L) {
        nt <- table(na_samp)
        miss_n[names(nt)] <- miss_n[names(nt)] + as.integer(nt)
      }
    }
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  p_miss <- ifelse(total_n > 0L, miss_n / total_n, 1)
  keep <- samples[p_miss < threshold]

  if (length(keep) == 0L)
    cli::cli_abort("All samples removed by vcf_filter_missing(threshold = {threshold})")

  vcf_arrow <- .vcf_filter_columns(vcf_arrow, keep, f_invar, verbose)

  return(vcf_arrow)
}
