#' @title vcf_filter_invariant
#'
#' @description
#' Remove invariant loci from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param strict_poly -> polymorphism at the genotype level (TRUE) or the
#'   allele level (FALSE), default FALSE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function removes invariant loci from a VCFArrow object,
#' returning a new VCFArrow object.
#' This might be desirable after subsetting a VCFArrow object by individuals.
#' A variant is kept if it has at least one called genotype and is
#' polymorphic among its called genotypes.  With strict_poly = FALSE (the
#' default), polymorphic means that not all called alleles are the same,
#' e.g. both REF and ALT alleles occur; a variant with a single called
#' heterozygote is kept.  With strict_poly = TRUE, polymorphic means that not
#' all called genotypes are the same, so a variant where every called sample
#' is heterozygous is removed although both alleles occur (such sites are
#' often paralogs).
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_filter_invariant(vcf)
#' # also remove variants where all called genotypes are heterozygous
#' vcf_filter_invariant(vcf, strict_poly = TRUE)
#'
#' @export
#'

vcf_filter_invariant <- function(vcf_arrow, strict_poly = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  idx <- .vcf_filter_index(vcf_arrow)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)

  cli::cli_alert_info("Applying invariant filter")

  parts <- .map_chunks(ffiles, .invariant_chunk,
                       shared = c(idx[c("pos", "samples")], strict_poly = strict_poly),
                       label = "Scanning chunk")
  # Per variant: count of called genotypes, sum of their allele_sum and sum
  # of squares (of genotypes or alleles, see .is_polymorphic())
  n_pass <- .merge_ranges(parts, "n", integer(idx$n_var))
  s1 <- .merge_ranges(parts, "s1", numeric(idx$n_var))
  s2 <- .merge_ranges(parts, "s2", numeric(idx$n_var))

  pass <- n_pass > 0L & .is_polymorphic(n_pass, s1, s2, strict_poly)
  keep <- vcf_arrow@variants$.row_id[pass]
  
  cli::cli_alert_info(
    "Retained {length(keep)} / {idx$n_var} variant{?s} \\
     (polymorphic at the {if (strict_poly) 'genotype' else 'allele'} level)"
  )

  # Record the removed variants, with their metadata and INFO, so they stay
  # recoverable (e.g. by vcf_bind(), which treats them as the object's own).
  n_removed <- idx$n_var - length(keep)
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

# Per-chunk count, sum of allele_sum and sum of squares (of genotypes or
# alleles, see .poly_squares()) per variant, over called genotypes
.invariant_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2"), shared)
  called <- !is.na(ch$a1) & !is.na(ch$a2)
  pos <- ch$pos[called]
  a1 <- ch$a1[called]
  a2 <- ch$a2[called]
  list(n = .range_count(pos), s1 = .range_sum(pos, a1 + a2),
       s2 = .range_sum(pos, .poly_squares(a1, a2, shared$strict_poly)))
}
