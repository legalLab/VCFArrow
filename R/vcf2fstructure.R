#' @title vcf2fstructure
#'
#' @description
#' Converts a VCFArrow object to a FastStructure infile
#'
#' @author Tomas Hrbek November 2026
#'
#' @inheritParams vcf2structure
#' @param method -> flag for Structure/FastStructure formats, default 'F'
#'   (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function is an alias of vcf2structure() that writes the FastStructure
#' format by default (method = 'F'); see vcf2structure() for the file layout.
#' Only biallelic SNVs are exported (multiallelic sites and indels are left
#' out), and every sample must be assigned to a group (see set_vcf_groups()),
#' also when keep_groups is NULL.
#' FastStructure adds the extension '.str' to the name given with --input,
#' so give out_file the extension .str (e.g. "data.str", read with
#' --input=data --format=str).
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2fstructure(vcf, out_file = tempfile(fileext = ".str"))
#'
#' @export
#'

vcf2fstructure <- function(vcf_arrow, out_file, keep_groups = NULL,
                           method = "F") {
  vcf2structure(vcf_arrow, out_file, keep_groups = keep_groups, method = method)
}
