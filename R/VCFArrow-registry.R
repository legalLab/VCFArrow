#' @description helper VCFArrow class functions
#'
#' @author Tomas Hrbek April 2026
#'
#' @details
#' These functions register and deregister new VCFArrow class instance
#'

.vcfarrow_registry <- new.env(parent = emptyenv())

.register_vcfarrow <- function(path) {
  if (!exists(path, envir = .vcfarrow_registry)) {
    .vcfarrow_registry[[path]] <- 1L
  } else {
    .vcfarrow_registry[[path]] <- .vcfarrow_registry[[path]] + 1L
  }
}

.deregister_vcfarrow <- function(path) {
  if (!exists(path, envir = .vcfarrow_registry)) return()
  .vcfarrow_registry[[path]] <- .vcfarrow_registry[[path]] - 1L
  if (.vcfarrow_registry[[path]] <= 0) {
    unlink(path, recursive = TRUE, force = TRUE)
    rm(list = path, envir = .vcfarrow_registry)
  }
}

.attach_finalizer <- function(vcf_obj) {
  path <- vcf_obj@path
  ptr <- new.env()  # lightweight holder
  reg.finalizer(ptr, function(e) {
    .deregister_vcfarrow(path)
  }, onexit = TRUE)
  attr(vcf_obj, "finalizer_env") <- ptr
  vcf_obj
}
