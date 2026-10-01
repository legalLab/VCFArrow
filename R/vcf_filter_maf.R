#' @title vcf_filter_maf
#'
#' @description
#' Remove loci below MAF threshold from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> decimal MAF threshold, default 0.05 (numeric)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes loci below MAF threshold from a VCFArrow object,
#' returning a new VCFArrow object.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_maf(vcf, threshold = 0.05)
#'
#' @export
#'

vcf_filter_maf <- function(vcf_arrow, threshold = 0.05) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying MAF filter")

  parts <- .map_chunks(ffiles, .maf_chunk, shared = idx[c("pos", "samples")],
                       label = "Scanning chunk")
  n_called <- .merge_ranges(parts, "n", integer(idx$n_var))
  alt_sum <- .merge_ranges(parts, "alt", numeric(idx$n_var))  # sum of (a1+a2)

  af <- ifelse(n_called > 0L, alt_sum / (2L * n_called), NA_real_)
  maf <- pmin(af, 1 - af, na.rm = FALSE)
  pass <- !is.na(maf) & maf >= threshold
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} (MAF >= {threshold})"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}

# Per-chunk counts of called genotypes and sums of (a1 + a2) per variant
.maf_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2"), shared)
  called <- !is.na(ch$a1) & !is.na(ch$a2)
  pos <- ch$pos[called]
  list(n = .range_count(pos), alt = .range_sum(pos, ch$a1[called] + ch$a2[called]))
}
