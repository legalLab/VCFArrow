#' @description helper VCFArrow class functions
#'
#' @author Tomas Hrbek April 2026
#'
#' @details
#' This function makes a new VCFArrow class instance
#'

.new_vcfarrow <- function(header, info, format, variants, gt, samples, groups, path) {

  .register_vcfarrow(path)

  obj <- new("VCFArrow",
             header = header,
             info = info,
             format = format,
             variants = variants,
             gt = gt,
             samples = samples,
             groups = groups,
             path = path
  )

  obj <- .attach_finalizer(obj)

  return(obj)
}
