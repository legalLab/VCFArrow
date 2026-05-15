#' @title vcf_filter_adr
#'
#' @description
#' Correct or remove genotypes with > % ALT/REF ratio from a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param mode -> switch between 'correct' or 'remove' mode (character)
#' @param threshold -> decimal missing threshold, default 0.1 (numeric)
#' @param f_invar -> filter invariant loci flag, default TRUE (Boolean)
#'
#' @return subsetted VCFArrow object
#'
#' @details
#' This function either changes genotypes or makes genotypes missing in
#' a VCFArrow object if they have above/below threshold normalized ADR ratio,
#' returning a new VCFArrow object.
#' The ADR is calculated as ADR = ALT / (REF + ALT).
#' If ADR > threshold, the genotype becomes ALT homozygous.
#' If ADR < threshold, the genotype becomes REF homozygous.
#' By default will remove any loci that may have become invariant as the
#' result of the change/removal of genotypes
#'
#' @examples
#' vcf_filter_adr(vcf_arrow = my_vcf, mode = "correct", threshold = my_threshold, f_invar = TRUE)
#' vcf_filter_adr(my_vcf, "correct", my_threshold, TRUE)
#' vcf_filter_adr(my_vcf, "correct")
#'

vcf_filter_adr <- function(vcf_arrow, mode = c("correct", "remove"), threshold = 0.1, f_invar = TRUE) {

  mode <- match.arg(mode)

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  gt <- vcf_arrow@gt

  # set allele ratio imbalance flag and precompute direction
  gt <- gt |>
    dplyr::mutate(
      adr_flag = ADR < threshold | ADR > (1 - threshold),
      adr_dir = dplyr::case_when(
        ADR <= threshold ~ 0L,
        ADR >= (1 - threshold) ~ 1L,
        TRUE ~ NA_integer_
      )
    )

  # apply genotype correction based on direction of bias
  if (mode == "correct") {
    vcf_arrow@gt <- gt |>
      dplyr::mutate(
        a1 = dplyr::if_else(adr_flag, adr_dir, a1),
        a2 = dplyr::if_else(adr_flag, adr_dir, a2)
      )
  }

  # apply 'missing' mask
  if (mode == "remove") {
    vcf_arrow@gt <- gt |>
      dplyr::mutate(
        a1 = dplyr::if_else(adr_flag, NA_integer_, a1),
        a2 = dplyr::if_else(adr_flag, NA_integer_, a2)
      )
  }

  # remove invariant loci
  if (f_invar) {
    vcf_arrow <- vcf_filter_invariant(vcf_arrow)
  }

  return(vcf_arrow)
}
