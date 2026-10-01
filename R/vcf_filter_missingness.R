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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_missingness(vcf, threshold = 0.2)
#'
#' @export
#'

vcf_filter_missingness <- function(vcf_arrow, threshold = 0.1, verbose = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying locus missingness filter")

  parts <- .map_chunks(ffiles, .missingness_chunk, shared = idx[c("pos", "samples")],
                       label = "Scanning chunk")
  total_n <- .merge_ranges(parts, "n", integer(idx$n_var))
  miss_n <- .merge_ranges(parts, "miss", integer(idx$n_var))

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

# Per-chunk counts of genotypes and of missing (a1 NA) genotypes per variant
.missingness_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, "a1", shared)
  list(n = .range_count(ch$pos), miss = .range_count(ch$pos[is.na(ch$a1)]))
}
