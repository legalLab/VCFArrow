#' @title vcf_filter_hets
#'
#' @description
#' Remove loci above a heterozigosity threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> heterozigosity threshold, default 0.5 (numeric)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci above a heterozigosity threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#' High heterozigosities are indicative of potential paralogs.
#'
#' @examples
#' vcf_filter_hets(vcf_arrow = my_vcf, threshold = 0.5)
#' vcf_filter_hets(my_vcf, 0.5)
#' vcf_filter_hets(my_vcf)
#'
#' @export
#'

vcf_filter_hets <- function(vcf_arrow, threshold = 0.5) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  n_het <- integer(idx$n_var)
  n_called <- integer(idx$n_var)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_progress_bar("Scanning chunk (heterozygosity)", total = length(ffiles))
  for (fpath in ffiles) {
    chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))
    chunk <- chunk[idx$lv[chunk$.row_id] & chunk$sample %in% idx$samples, , drop = FALSE]
    if (nrow(chunk) > 0L) {
      called <- !is.na(chunk$a1) & !is.na(chunk$a2)
      if (any(called)) {
        sub <- chunk[called, , drop = FALSE]
        pos <- idx$col_idx[as.character(sub$.row_id)]
        n_called <- n_called + tabulate(pos, nbins = idx$n_var)
        het_pos <- pos[sub$a1 != sub$a2]
        if (length(het_pos) > 0L) {
          n_het <- n_het + tabulate(het_pos, nbins = idx$n_var)
        }
      }
    }
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  het_rate <- ifelse(n_called > 0L, n_het / n_called, NA_real_)
  # n_called == 0 → excluded (matches the original's implicit-omission behaviour)
  pass <- n_called > 0L & het_rate < threshold
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variants \\
     (heterozygosity rate < {threshold})"
  )
  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
