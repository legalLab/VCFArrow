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
#' @export
#'

vcf_copy <- function(vcf_arrow) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  old_files <- list.files(vcf_arrow@path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(old_files) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")

  new_path <- tempfile("arrow_vcf_copy_")
  dir.create(new_path)

  new_files <- file.path(new_path, basename(old_files))
  ok <- file.copy(old_files, new_files)

  if (!all(ok)) {
    # Clean up the partially-copied directory before erroring, so a failed
    # vcf_copy() doesn't leave an orphaned, unregistered temp directory.
    unlink(new_path, recursive = TRUE, force = TRUE)
    cli::cli_abort(
      "Failed to copy {sum(!ok)} of {length(old_files)} feather file{?s} \\
       to {.path {new_path}} — check disk space and permissions."
    )
  }

  # suppressWarnings: Arrow's "Invalid metadata$r" warning when re-parsing
  # R-specific IPC schema annotations on freshly copied files — benign, see
  # the note on this same pattern in vcf_bind_sparse().
  gt_arrow <- suppressWarnings(arrow::open_dataset(new_path, format = "feather"))

  new_vcfarrow <- .new_vcfarrow(
    header = vcf_arrow@header,
    info = vcf_arrow@info,
    format = vcf_arrow@format,
    variants = vcf_arrow@variants,
    gt = gt_arrow,
    samples = vcf_arrow@samples,
    groups = vcf_arrow@groups,
    path = new_path
  )

  return(new_vcfarrow)
}
