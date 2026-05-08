#' @description helper VCFArrow class functions
#'
#' @author Tomas Hrbek April 2026
#'
#' @details
#' This function defines the VCFArrow S4 class
#'

setClass(
  "VCFArrow",
  slots = list(
    header = "character", # complete VCF header
    info = "character", # INFO field
    format = "data.frame", # FORMAT field
    variants = "data.frame", # CHROM, POS, ID, REF, ALT, QUAL, FILTER
    gt = "ANY", # Arrow dataset
    samples = "character", # sample names from gt columns
    groups = "character", # sample group names
    path = "character" # dataset location - lazy loading
  )
)
