#' @title vcf_bind_sparse
#'
#' @description
#' Bind two or more VCFArrow objects into a new VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param ... -> a collection of VCFArrow objects
#'
#' @return VCFArrow object
#'
#' @details
#' This function binds two or more VCFArrow objects, returning new VCFArrow object.
#' The VCFArrow objects must have the same SNPs, and must have unique individuals.
#'
#' @examples
#' vcf_bind_sparse(my_vcf1, my_vcf2, other_vcf, ...)
#'

vcf_bind_sparse <- function(...) {

  vcfs <- list(...)

  # make unified variants
  unified_variants <- dplyr::bind_rows(
    lapply(vcfs, function(x) x@variants[, c("CHROM","POS","REF","ALT")])
  ) |>
    dplyr::distinct() |>
    dplyr::arrange(CHROM, POS)

  unified_variants$.row_id <- seq_len(nrow(unified_variants))

  # remap each VCF
  remap_gt <- function(vcf) {

    mapping <- vcf@variants |>
      dplyr::select(.row_id, CHROM, POS, REF, ALT) |>
      dplyr::inner_join(unified_variants,
                        by = c("CHROM","POS","REF","ALT"),
                        suffix = c("_old", "_new")) |>
      dplyr::select(old_id = .row_id_old, new_id = .row_id_new)

    vcf@gt |>
      dplyr::inner_join(mapping, by = c(".row_id" = "old_id")) |>
      dplyr::mutate(.row_id = new_id) |>
      dplyr::select(-new_id)
  }

  merged_gt <- dplyr::bind_rows(lapply(vcfs, remap_gt))

  merged_samples <- unlist(lapply(vcfs, function(x) x@samples), use.names = FALSE)

  # check duplicates
  if (anyDuplicated(vcf_new@samples)) {
    cli::cli_abort("Duplicate sample names detected after binding")
  }

  vcf_new <- vcfs[[1]]
  vcf_new@variants <- unified_variants
  vcf_new@gt <- merged_gt
  vcf_new@samples <- merged_samples

  return(vcf_new)
}
