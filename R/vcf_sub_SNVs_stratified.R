#' @title vcf_sub_SNVs_stratified
#'
#' @description
#' Randomly subsets SNVs from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param n_SNVs -> number of SNVs to subset, default 1000 (integer)
#' @param seed -> random number generator seed, default NULL (integer)
#'
#' @return VCFArrow object
#'
#' @details
#' This function subsets a VCFArrow object to specific number of SNVs,
#' returning new VCFArrow object.
#' The subsampling is stratified, i.e. the same proportion of SNVs per CHROM.
#' Each CHROM contributes its proportional share of n_SNVs, rounded down;
#' the SNVs still needed to reach n_SNVs come one each from the CHROMs with
#' the largest remainders (ties broken at random), so exactly n_SNVs SNVs
#' are returned. When most CHROMs hold a single variant (e.g. DiscoSnp-RAD
#' output), this amounts to simple random sampling.
#' The seed for random number generator is automatically generated unless
#' specified.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_sub_SNVs_stratified(vcf, n_SNVs = 1000, seed = 42)
#'
#' @export
#'

vcf_sub_SNVs_stratified <- function(vcf_arrow, n_SNVs = 1000, seed = NULL) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  variants <- vcf_arrow@variants
  n_vars <- nrow(variants)

  if (n_SNVs >= n_vars) {
    cli::cli_alert_warning("Number of SNVs to subsample ({n_SNVs}) is greater than number of variants available ({n_vars});
                           returning original VCFArrow object")
    return(vcf_arrow)
  }

  # with a seed, sample reproducibly without changing the session's random
  # number state (.with_seed() restores it afterwards)
  draw <- function() {
    # number of SNVs per CHROM: the whole part of its proportional share
    # n_SNVs * size / n_vars, plus one for the CHROMs with the largest
    # remainders (ties at random) until the total is n_SNVs. Rounding each
    # share down alone (as slice_sample(prop =) does) gives CHROMs with fewer
    # than n_vars / n_SNVs variants none, i.e. nothing at all when nearly
    # every variant has its own CHROM (e.g. DiscoSnp). Exact in doubles.
    # vctrs::vec_count() rather than dplyr::count(): counting by grouping is
    # slow with millions of CHROMs (15 s vs 0.5 s for 2 million)
    quota <- vctrs::vec_count(variants$CHROM, sort = "location") |>
      dplyr::rename(CHROM = key, size = count) |>
      dplyr::mutate(
        share = n_SNVs * as.numeric(size),
        take = share %/% n_vars,
        u = stats::runif(dplyr::n())
      ) |>
      dplyr::arrange(dplyr::desc(share %% n_vars), u) |>
      dplyr::mutate(take = take + (dplyr::row_number() <= n_SNVs - sum(take))) |>
      dplyr::select(CHROM, take)

    # take variants of each CHROM, chosen at random: order the variants by
    # CHROM and a random key, and keep the first `take` of each CHROM.
    # slice_sample() cannot do this (its n is the same for every group), and
    # the rank within CHROM is computed ungrouped: per-group evaluation is
    # slow when CHROM has nearly one value per variant (see .snv_blocks()).
    variants |>
      dplyr::select(CHROM, .row_id) |>
      dplyr::mutate(u = stats::runif(dplyr::n())) |>
      dplyr::arrange(CHROM, u) |>
      dplyr::mutate(rank = dplyr::row_number() - match(CHROM, CHROM) + 1L) |>
      dplyr::inner_join(quota, by = "CHROM") |>
      dplyr::filter(rank <= take) |>
      dplyr::pull(.row_id) |>
      sort()
  }
  keep <- if (is.null(seed)) draw() else .with_seed(seed, draw())

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}
