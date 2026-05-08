#' @description helper VCFArrow class functions to run on load
#'
#' @author Tomas Hrbek April 2026
#'
#' @details
#' These functions perform GC and register new VCFArrow class instance
#'

.cleanup_old_vcfarrow <- function() {
  tmp <- tempdir()
  dirs <- list.dirs(tmp, recursive = FALSE)

  vcf_dirs <- dirs[grepl("arrow_vcf_", dirs)]

  for (d in vcf_dirs) {
    unlink(d, recursive = TRUE, force = TRUE)
  }
}

.onLoad <- function(libname, pkgname) {
  .cleanup_old_vcfarrow()
}

.onLoad <- function(libname, pkgname) {
  assign(".vcfarrow_registry", new.env(parent = emptyenv()), envir = parent.env(environment()))
}
