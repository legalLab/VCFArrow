#' @title VCFArrow S4 class
#'
#' @description
#' S4 class holding VCF data. Variant metadata are kept in memory while
#' genotypes are stored on disk as an Apache Arrow dataset and loaded lazily.
#' Objects are created with [read_vcf()].
#'
#' @author Tomas Hrbek April 2026
#'
#' @slot header Complete VCF header (character).
#' @slot info INFO field of each variant (character).
#' @slot format FORMAT field of each variant (data.frame).
#' @slot variants CHROM, POS, ID, REF, ALT, QUAL, FILTER and precalculated
#'   per-variant metrics (data.frame).
#' @slot gt Arrow dataset of genotypes in long format.
#' @slot samples Sample names (character).
#' @slot groups Group assignment of each sample (character).
#' @slot path Location of the on-disk Arrow dataset (character).
#' @slot finalizer_env Environment used to clean up `path` on garbage collection.
#' @slot invariant_removed Variants removed as invariant, with their metadata
#'   (data.frame).
#'
#' @exportClass VCFArrow
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
    path = "character", # dataset location - lazy loading
    finalizer_env = "environment", # slot for maintaining info for GC
    invariant_removed = "data.frame" # for maintaining removed variants and their metadata
  )
)
