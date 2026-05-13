#' @title vcf_filter_rows
#'
#' @description
#' Unified API for filtering VCFArrow objects by row IDs.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep -> rows to keep (numeric or logical)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#' @param verbose -> report removed samples and final % missing data, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes all rows in the 'keep' parameter from a VCFArrow object,
#' returning a new VCFArrow object.
#' Optionally will remove any loci that may have become invariant as the
#' result of the removal of samples.
#' Optionally will report removed samples, final % missing data, and number of
#' retained samples after sample filtering.
#'
#' @examples
#' vcf_filter_rank(vcf_arrow = my_vcf, keep = rows_to_keep)
#' vcf_filter_rank(my_vcf, rows_to_keep)
#'

vcf_filter_columns <- function(vcf_arrow, keep, f_invar, verbose) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  s <- vcf_arrow@samples

  # normalize input → row_ids
  if (is.logical(keep)) {
    if (length(keep) != lengths(s)) {
      cli::cli_abort("Logical filter must match number of samples")
    }
    keep_ids <- s[keep]
  } else if (is.numeric(keep)) {
    keep_ids <- s[keep]
  } else if (is.character(keep)) {
    keep_ids <- intersect(s, keep)
  } else {
    cli::cli_abort("The 'keep' parameter must be character, logical or numeric")
  }

  # enforce character type
  keep_ids <- as.character(keep_ids)

  # get removed samples
  removed_samples <- setdiff(s, keep)

  # update metadata
  vcf_arrow@samples <- keep_ids
  idx <- match(keep_ids, s)
  vcf_arrow@groups <- vcf_arrow@groups[idx]

  # Arrow filter
  keep_tbl <- arrow::arrow_table(sample = keep_ids)

  vcf_arrow@gt <- vcf_arrow@gt |>
    dplyr::semi_join(keep_tbl, by = "sample")

  # remove invariant loci
  if (f_invar) {
    vcf_arrow <- vcf_filter_invariant(vcf_arrow)
  }

  # reporting
  if (verbose) {

    if (length(removed_samples) > 0) {
      cli::cli_alert_info("Removed samples: {removed_samples}")
    }

    total_stats <- vcf_arrow@gt |>
      dplyr::summarise(
        total_na = sum(is.na(a1)),
        total_n = dplyr::n()
      ) %>%
      dplyr::collect()

    p_missing <- total_stats$total_na / total_stats$total_n

    cli::cli_alert_info("Final % missing data in VCF is {sprintf('%.2f', 100 * p_missing)}%")
    cli::cli_alert_info("Number of samples in VCF after filtering is {length(keep_ids)}")
  }

  return(vcf_arrow)
}
