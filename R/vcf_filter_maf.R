#' @title vcf_filter_maf
#'
#' @description
#' Remove loci below MAF threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal MAF threshold, default 0.05 (numeric)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci below MAF threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_maf(vcf_arrow = my_vcf, threshold = 0.05)
#' vcf_filter_maf(my_vcf, 0.05)
#' vcf_filter_maf(my_vcf)
#'
#' @export
#'

vcf_filter_maf <- function(vcf_arrow, threshold = 0.05) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  alt_sum <- integer(idx$n_var)
  n_called <- integer(idx$n_var)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying MAF filter")

  cli::cli_progress_bar("Scanning chunk", total = length(ffiles))
  for (fpath in ffiles) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))
    chunk <- chunk[idx$lv[chunk$.row_id] & chunk$sample %in% idx$samples, , drop = FALSE]
    if (nrow(chunk) > 0L) {
      called <- !is.na(chunk$a1) & !is.na(chunk$a2)
      if (any(called)) {
        sub <- chunk[called, , drop = FALSE]
        pos <- idx$col_idx[as.character(sub$.row_id)]
        n_called <- n_called + tabulate(pos, nbins = idx$n_var)
        # tapply gives sum of (a1+a2) per position
        rs <- tapply(sub$a1 + sub$a2, pos, sum)
        alt_sum[as.integer(names(rs))] <- alt_sum[as.integer(names(rs))] + as.integer(rs)
      }
    }
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  af <- ifelse(n_called > 0L, alt_sum / (2L * n_called), NA_real_)
  maf <- pmin(af, 1 - af, na.rm = FALSE)
  pass <- !is.na(maf) & maf >= threshold
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variants (MAF >= {threshold})"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
