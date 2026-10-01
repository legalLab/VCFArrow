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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_coverage(vcf, threshold = 10)
#'
#' @export
#'

vcf_filter_coverage <- function(vcf_arrow, threshold = 10) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying read coverage filter")

  parts <- .map_chunks(ffiles, .coverage_chunk,
                       shared = c(idx[c("pos", "samples")], threshold = threshold),
                       label = "Scanning chunk")
  # Per variant: count of DP-passing called genotypes, and sum and sum of
  # squares of their allele_sum: all values of a variant are equal
  # (min == max) exactly when n * sum(x^2) == sum(x)^2.
  dp_pass <- .merge_ranges(parts, "n", integer(idx$n_var))
  s1 <- .merge_ranges(parts, "s1", numeric(idx$n_var))
  s2 <- .merge_ranges(parts, "s2", numeric(idx$n_var))

  # Pass: at least 2 DP-passing called genotypes AND not monomorphic
  pass <- dp_pass >= 2L & dp_pass * s2 != s1 * s1
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} \\
     (polymorphic with DP >= {threshold})"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}

# Per-chunk count, sum and sum of squares of allele_sum per variant, over
# called genotypes with DP >= threshold
.coverage_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2", "DP"), shared)
  ok <- !is.na(ch$a1) & !is.na(ch$a2) & !is.na(ch$DP) & ch$DP >= shared$threshold
  pos <- ch$pos[ok]
  gs <- ch$a1[ok] + ch$a2[ok]   # 0/1/2 for hom-ref/het/hom-alt
  list(n = .range_count(pos), s1 = .range_sum(pos, gs), s2 = .range_sum(pos, gs * gs))
}
