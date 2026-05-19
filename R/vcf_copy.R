#' @title vcf_copy
#'
#' @description
#' Copy a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#'
#' @return VCFArrow object
#'
#' @details
#' This function copies a VCFArrow object.
#' It makes an actual physical copy of the object.
#'
#' @examples
#' vcf_copy(vcf_arrow = my_vcf)
#' vcf_copy(my_vcf)
#'

vcf_copy <- function(vcf_arrow) {

  .register_vcfarrow(vcf_arrow@path)

  vcf_arrow2 <- vcf_arrow
  vcf_arrow2 <- .attach_finalizer(vcf_arrow2)

  return(vcf_arrow2)
}
