#' @title vcf_bind
#'
#' @description
#' Bind two or more VCFArrow objects into a new VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param ... -> a collection of VCFArrow objects
#' @param .check -> check that VCFArrow objects are compatible, default TRUE (Boolean)
#'
#' @return VCFArrow object
#'
#' @details
#' This function binds two or more VCFArrow objects, returning new VCFArrow object.
#' The VCFArrow objects must have the same SNPs, and must have unique individuals.
#' Default is to check compatibility between VCFArrow objects;
#' this default can but is not meant to be changed.
#' If VCFArrow objects are not compatible but need to be merged,
#' use vcf_bind_sparse().
#'
#' @examples
#' vcf_bind(my_vcf1, my_vcf2, other_vcf, ..., check = TRUE)
#' vcf_bind(my_vcf1, my_vcf2, other_vcf)
#'

vcf_bind <- function(..., .check = TRUE) {

  vcfs <- list(...)

  # validate inputs
  if (!all(vapply(vcfs, function(x) methods::is(x, "VCFArrow"), logical(1)))) {
    cli::cli_abort("At least one input is not a VCFArrow object")
  }

  if (length(vcfs) < 2) {
    cli::cli_abort("Provide at least two VCFArrow objects to merge")
  }

  # check variant compatibility
  if (.check) {
    ref_var <- vcfs[[1]]@variants

    for (i in seq_along(vcfs)[-1]) {
      variants <- vcfs[[i]]@variants

      # sanity check
      if (!all(variants$.row_id == ref_var$.row_id)) {
        stop("Variants (.row_id) do not match between VCFArrow objects")
      }
      if (!all(variants$CHROM == ref_var$CHROM & variants$POS == ref_var$POS)) {
        cli::cli_abort("Variant coordinates do not match across inputs")
      }
    }
  }

  # template
  vcf_new <- vcfs[[1]]

  # bind genotype tables
  vcf_new@gt <- do.call(
    dplyr::union_all,
    lapply(vcfs, function(x) x@gt)
  )

  # merge samples
  vcf_new@samples <- unlist(
    lapply(vcfs, function(x) x@samples),
    use.names = FALSE
  )

  # merge groups
  vcf_new@groups <- unlist(
    lapply(vcfs, function(x) x@groups),
    use.names = FALSE
  )

  # check duplicates
  if (anyDuplicated(vcf_new@samples)) {
    cli::cli_abort("Duplicate sample names detected after binding")
  }

  return(vcf_new)
}
