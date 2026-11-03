#' @title vcf_filter_coverage
#'
#' @description
#' Remove genotypes and variants below a read depth (DP) threshold from a
#' VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param threshold -> DP threshold, default 10 (integer)
#' @param strict -> set called genotypes with DP below the threshold to
#'   missing, default TRUE (Boolean)
#' @param strict_poly -> polymorphism at the genotype level (TRUE) or the
#'   allele level (FALSE), default FALSE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function keeps the variants that have at least two called genotypes
#' with DP >= threshold and are polymorphic among them, returning a new
#' VCFArrow object.
#' With strict = TRUE (the default), called genotypes with DP below the
#' threshold, or without DP, are also set to missing (as vcftools --minDP
#' does), so every called genotype of the result has DP >= threshold; the
#' number of genotypes set to missing is reported.  With strict = FALSE the
#' genotypes are not changed, and the kept variants can still contain called
#' genotypes with DP below the threshold.
#' With strict_poly = FALSE (the default), polymorphic means that not all
#' of these genotypes' alleles are the same, e.g. both REF and ALT alleles
#' occur.  With strict_poly = TRUE, it means that not all of these genotypes
#' are the same, so a variant whose genotypes with DP >= threshold are all
#' heterozygous is removed although both alleles occur.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_coverage(vcf, threshold = 10)
#' # keep the genotypes, filter variants only
#' vcf_filter_coverage(vcf, threshold = 10, strict = FALSE)
#' # also remove variants whose well-covered genotypes are all heterozygous
#' vcf_filter_coverage(vcf, threshold = 10, strict_poly = TRUE)
#'
#' @export
#'

vcf_filter_coverage <- function(vcf_arrow, threshold = 10, strict = TRUE,
                                strict_poly = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying read coverage filter")

  shared <- c(idx[c("pos", "samples")], threshold = threshold,
              strict_poly = strict_poly)
  if (strict) {
    tmp_dir <- tempfile("arrow_vcf_coverage_")
    dir.create(tmp_dir)
    tasks <- lapply(seq_along(ffiles), function(i) list(
      fpath = ffiles[[i]], out = file.path(tmp_dir, paste0("chunk_", i, ".arrow"))
    ))
    parts <- .map_chunks(tasks, .coverage_strict_chunk, shared = shared,
                         label = "Rewriting chunk")
    masked <- sum(vapply(parts, `[[`, numeric(1), "masked"))
    called <- sum(vapply(parts, `[[`, numeric(1), "called"))
    cli::cli_alert_info(
      "Set {masked} / {called} called genotype{?s} with DP < {threshold} to missing"
    )
    vcf_arrow <- .new_vcfarrow(
      header = vcf_arrow@header,
      info = vcf_arrow@info,
      format = vcf_arrow@format,
      variants = vcf_arrow@variants,
      gt = suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather")),
      samples = vcf_arrow@samples,
      groups = vcf_arrow@groups,
      path = tmp_dir,
      invariant_removed = vcf_arrow@invariant_removed,
      loci = vcf_arrow@loci
    )
  } else {
    parts <- .map_chunks(ffiles, .coverage_chunk, shared = shared,
                         label = "Scanning chunk")
  }
  # Per variant: count of DP-passing called genotypes, sum of their
  # allele_sum and sum of squares (of genotypes or alleles, see
  # .is_polymorphic())
  dp_pass <- .merge_ranges(parts, "n", integer(idx$n_var))
  s1 <- .merge_ranges(parts, "s1", numeric(idx$n_var))
  s2 <- .merge_ranges(parts, "s2", numeric(idx$n_var))

  # Pass: at least 2 DP-passing called genotypes AND not monomorphic
  pass <- dp_pass >= 2L & .is_polymorphic(dp_pass, s1, s2, strict_poly)
  keep <- vcf_arrow@variants$.row_id[pass]

  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} \\
     (polymorphic at the {if (strict_poly) 'genotype' else 'allele'} level \\
     with DP >= {threshold})"
  )

  # apply filter using unified API
  vcf_arrow <- .vcf_filter_rows(vcf_arrow, keep)

  return(vcf_arrow)
}

# Per-chunk count, sum of allele_sum and sum of squares (of genotypes or
# alleles, see .poly_squares()) per variant, over called genotypes with
# DP >= threshold
.coverage_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2", "DP"), shared)
  ok <- !is.na(ch$a1) & !is.na(ch$a2) & !is.na(ch$DP) & ch$DP >= shared$threshold
  pos <- ch$pos[ok]
  a1 <- ch$a1[ok]
  a2 <- ch$a2[ok]
  list(n = .range_count(pos), s1 = .range_sum(pos, a1 + a2),
       s2 = .range_sum(pos, .poly_squares(a1, a2, shared$strict_poly)))
}

# As .coverage_chunk(), and also rewrite the chunk to task$out with the called
# genotypes with DP < threshold (or without DP) set to missing, so that the
# called genotypes left are exactly those counted.  All rows are rewritten;
# the counts (and the numbers of called and masked genotypes) cover only the
# rows of live variants and retained samples.
.coverage_strict_chunk <- function(task, shared) {
  chunk <- arrow::read_feather(task$fpath, as_data_frame = FALSE)
  a1 <- as.vector(chunk$a1)
  a2 <- as.vector(chunk$a2)
  dp <- as.vector(chunk$DP)
  called <- !is.na(a1) & !is.na(a2)
  mask <- called & (is.na(dp) | dp < shared$threshold)
  pos <- .match_row_id(as.vector(chunk$.row_id), shared$pos)
  live <- !is.na(pos) & !is.na(.sample_index(chunk, shared$samples))
  ok <- live & called & !mask
  out <- list(n = .range_count(pos[ok]), s1 = .range_sum(pos[ok], a1[ok] + a2[ok]),
              s2 = .range_sum(pos[ok], .poly_squares(a1[ok], a2[ok], shared$strict_poly)),
              masked = sum(live & mask), called = sum(live & called))

  a1[mask] <- NA_integer_
  a2[mask] <- NA_integer_
  chunk[["a1"]] <- a1
  chunk[["a2"]] <- a2
  chunk <- .set_fmt_gt(chunk, mask, a1)
  arrow::write_feather(chunk, task$out)
  chunk <- NULL
  .release_arrow_memory()
  out
}
