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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf |>
#'   vcf_filter_missingness(threshold = 0.2) |>
#'   vcf_filter_missing(threshold = 0.3)
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

  parts <- .map_chunks(ffiles, .missing_chunk, shared = idx[c("pos", "samples")],
                       label = "Scanning chunk")
  for (p in parts) {
    total_n <- total_n + p$n
    miss_n <- miss_n + p$miss
  }

  p_miss <- ifelse(total_n > 0L, miss_n / total_n, 1)
  keep <- samples[p_miss < threshold]

  if (length(keep) == 0L)
    cli::cli_abort("All samples removed by vcf_filter_missing(threshold = {threshold})")

  vcf_arrow <- .vcf_filter_columns(vcf_arrow, keep, f_invar, verbose)

  return(vcf_arrow)
}

# Per-chunk counts of genotypes and of missing (a1 NA) genotypes per sample
.missing_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, "a1", shared)
  n_s <- length(shared$samples)
  list(n = tabulate(ch$s, n_s), miss = tabulate(ch$s[is.na(ch$a1)], n_s))
}
