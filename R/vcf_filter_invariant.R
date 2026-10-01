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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_invariant(vcf)
#'
#' @export
#'

vcf_filter_invariant <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying invariant filter")

  parts <- .map_chunks(ffiles, .invariant_chunk, shared = idx[c("pos", "samples")],
                       label = "Scanning chunk")
  # Per variant: count of called genotypes, and sum and sum of squares of
  # their allele_sum: all values of a variant are equal (min == max) exactly
  # when n * sum(x^2) == sum(x)^2.
  n_pass <- .merge_ranges(parts, "n", integer(idx$n_var))
  s1 <- .merge_ranges(parts, "s1", numeric(idx$n_var))
  s2 <- .merge_ranges(parts, "s2", numeric(idx$n_var))

  pass <- n_pass > 0L & n_pass * s2 != s1 * s1
  keep <- vcf_arrow@variants$.row_id[pass]

  n_removed <- idx$n_var - length(keep)
  if (n_removed > 0L)
    cli::cli_alert_info(
      "Removed {n_removed} invariant variant{?s}; {length(keep)} retained."
    )

  # Record the removed variants, with their metadata and INFO, so they stay
  # recoverable (e.g. by vcf_bind(), which treats them as the object's own).
  if (n_removed > 0L) {
    removed <- as.data.frame(vcf_arrow@variants)[!pass, , drop = FALSE]
    if (length(vcf_arrow@info) == nrow(vcf_arrow@variants))
      removed$.info_str <- vcf_arrow@info[!pass]
    else
      removed$.info_str <- rep(NA_character_, nrow(removed))
    prev <- vcf_arrow@invariant_removed
    rec <- if (nrow(prev) > 0L) rbind(prev, removed[names(prev)]) else removed
    vcf_arrow@invariant_removed <- rec[!duplicated(rec$.row_id), , drop = FALSE]
  }

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}

# Per-chunk count, sum and sum of squares of allele_sum per variant, over
# called genotypes
.invariant_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2"), shared)
  called <- !is.na(ch$a1) & !is.na(ch$a2)
  pos <- ch$pos[called]
  gs <- ch$a1[called] + ch$a2[called]  # 0 / 1 / 2
  list(n = .range_count(pos), s1 = .range_sum(pos, gs), s2 = .range_sum(pos, gs * gs))
}
