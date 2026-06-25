#' @title vcf_filter_coverage
#'
#' @description
#' Remove genotypes below read DP threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> DP threshold, default 10 (integer)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes genotypes below a DP threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' vcf_filter_coverage(vcf_arrow = my_vcf, threshold = 10)
#' vcf_filter_coverage(my_vcf, 10)
#' vcf_filter_coverage(my_vcf)
#'
#' @export
#'

vcf_filter_coverage <- function(vcf_arrow, threshold = 10) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  min_gs <- rep(3L,  idx$n_var)  # min allele_sum per variant; init > max possible (2)
  max_gs <- rep(-1L, idx$n_var)  # max allele_sum per variant; init < min possible (0)
  dp_pass <- integer(idx$n_var)   # count of DP-passing called genotypes
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_progress_bar("Scanning chunk (coverage)", total = length(ffiles))
  for (fpath in ffiles) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2", "DP"))
    chunk <- chunk[idx$lv[chunk$.row_id] & chunk$sample %in% idx$samples, , drop = FALSE]
    if (nrow(chunk) > 0L) {
      # Apply coverage threshold: keep only called genotypes with DP >= threshold
      ok <- !is.na(chunk$a1) & !is.na(chunk$a2) &
        !is.na(chunk$DP) & chunk$DP >= threshold
      if (any(ok)) {
        sub <- chunk[ok, , drop = FALSE]
        pos <- idx$col_idx[as.character(sub$.row_id)]
        gs <- sub$a1 + sub$a2   # 0/1/2 for hom-ref/het/hom-alt
        dp_pass <- dp_pass + tabulate(pos, nbins = idx$n_var)
        mn <- tapply(gs, pos, min)
        mx <- tapply(gs, pos, max)
        p <- as.integer(names(mn))
        min_gs[p] <- pmin(min_gs[p], as.integer(mn))
        max_gs[p] <- pmax(max_gs[p], as.integer(mx))
      }
    }
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # Pass: at least 2 DP-passing called genotypes AND not monomorphic
  pass <- dp_pass >= 2L & min_gs != max_gs
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variants \\
     (polymorphic with DP >= {threshold})"
  )
  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
