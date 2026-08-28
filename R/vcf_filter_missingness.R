#' @title vcf_filter_missingness
#'
#' @description
#' Remove loci above missingness threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal missingness threshold, default 0.1 (numeric)
#' @param verbose -> flag to report total % missing data after filtering, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci above missingness threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' Missingness is a locus focused metric, i.e. missing data per locus.
#'
#' @examples
#' vcf_filter_missingness(vcf_arrow = my_vcf, threshold = p_miss)
#' vcf_filter_missingness(my_vcf, p_miss, verbose = FALSE)
#' vcf_filter_missingness(my_vcf, p_miss)
#'
#' @export
#'

vcf_filter_missingness <- function(vcf_arrow, threshold = 0.1, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  miss_n <- integer(idx$n_var)  # NA count per variant position
  total_n <- integer(idx$n_var)  # row count per variant position
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying locus missingness filter")

  cli::cli_progress_bar("Scanning chunk", total = length(ffiles))
  for (fpath in ffiles) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1"))
    chunk <- chunk[idx$lv[chunk$.row_id] & chunk$sample %in% idx$samples, , drop = FALSE]
    if (nrow(chunk) > 0L) {
      pos <- idx$col_idx[as.character(chunk$.row_id)]
      tt <- tabulate(pos, nbins = idx$n_var)
      total_n <- total_n + tt
      na_pos <- pos[is.na(chunk$a1)]
      if (length(na_pos) > 0L)
        miss_n <- miss_n + tabulate(na_pos, nbins = idx$n_var)
    }
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  p_miss <- ifelse(total_n > 0L, miss_n / total_n, 1)
  keep_pos <- which(p_miss <= threshold)
  keep <- vcf_arrow@variants$.row_id[keep_pos]

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  if (verbose)
    cli::cli_alert_info(
      "Retained {length(keep)} / {idx$n_var} variant{?s} \\
       (per-variant missingness <= {threshold})"
    )

  return(vcf_arrow)
}
