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
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_hets(vcf, threshold = 0.5)
#'
#' @export
#'

vcf_filter_hets <- function(vcf_arrow, threshold = 0.5) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying heterozygosity filter")

  parts <- .map_chunks(ffiles, .hets_chunk, shared = idx[c("pos", "samples")],
                       label = "Scanning chunk")
  n_called <- .merge_ranges(parts, "n", integer(idx$n_var))
  n_het <- .merge_ranges(parts, "het", integer(idx$n_var))

  het_rate <- ifelse(n_called > 0L, n_het / n_called, NA_real_)
  # n_called == 0 → excluded (matches the original's implicit-omission behaviour)
  pass <- n_called > 0L & het_rate < threshold
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} \\
     (heterozygosity rate < {threshold})"
  )
  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}

# Per-chunk counts of called and heterozygous genotypes per variant
.hets_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2"), shared)
  called <- !is.na(ch$a1) & !is.na(ch$a2)
  pos <- ch$pos[called]
  list(n = .range_count(pos), het = .range_count(pos[ch$a1[called] != ch$a2[called]]))
}
