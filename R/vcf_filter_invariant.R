#' @title vcf_filter_invariant
#'
#' @description
#' Remove invariant loci from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes invariant loci from a VCFArrow object,
#' returning a new VCFArrow object.
#' This might be desirable after subsetting a VCFArrow object by individuals.
#'
#' @examples
#' vcf_filter_invariant(vcf_arrow = my_vcf)
#' vcf_filter_invariant(my_vcf)
#'
#' @export
#'

vcf_filter_invariant <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  min_gs <- rep( 3L, idx$n_var)  # above max possible allele_sum (2)
  max_gs <- rep(-1L, idx$n_var)  # below min possible allele_sum (0)
  n_pass <- integer(idx$n_var)  # called genotypes per variant
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_progress_bar("Scanning chunk (invariant filter)", total = length(ffiles))
  for (fpath in ffiles) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    chunk <- chunk[idx$lv[chunk$.row_id] & chunk$sample %in% idx$samples, , drop = FALSE]
    if (nrow(chunk) > 0L) {
      called <- !is.na(chunk$a1) & !is.na(chunk$a2)
      if (any(called)) {
        sub <- chunk[called, , drop = FALSE]
        pos <- idx$col_idx[as.character(sub$.row_id)]
        gs <- sub$a1 + sub$a2  # 0 / 1 / 2
        n_pass <- n_pass + tabulate(pos, nbins = idx$n_var)
        mn<- tapply(gs, pos, min)
        mx<- tapply(gs, pos, max)
        p <- as.integer(names(mn))
        min_gs[p] <- pmin(min_gs[p], as.integer(mn))
        max_gs[p] <- pmax(max_gs[p], as.integer(mx))
      }
    }
    chunk <- NULL; gc(verbose = FALSE, full = FALSE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  pass <- n_pass > 0L & min_gs != max_gs
  keep <- vcf_arrow@variants$.row_id[pass]

  n_removed <- idx$n_var - length(keep)
  if (n_removed > 0L)
    cli::cli_alert_info(
      "Removed {n_removed} invariant variant{?s}; {length(keep)} retained."
    )

  # apply filter using unified API
  vcf_arrow <- vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
