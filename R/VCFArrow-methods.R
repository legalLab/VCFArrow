#' @title Show method for VCFArrow
#'
#' @description Print a summary of a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param object A VCFArrow object
#'
#' @return Invisibly returns `object`; called for its side effect of printing
#'   a summary.
#'
#' @details
#' This function is a method of the VCFArrow S4 class
#' Method to show object content summary
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f)
#' show(vcf)
#'
#' @export
#'

setMethod(
  "show",
  "VCFArrow",
  function(object) {

    cat("\nAn object of class \"VCFArrow\"\n")

    # --- basic dimensions ---
    n_var <- tryCatch(nrow(object@variants), error = function(e) NA_integer_)
    n_samp <- length(object@samples)

    cat("\nDimensions:\n")
    cat("  Variants:", n_var, "\n")
    cat("  Samples: ", n_samp, "\n")

    cat("\nQuick stats:\n")
    cat("  Non-missing variants:", sum(!is.na(object@variants$POS)), "\n")

    has_phased <- any(grepl("_phased$", names(object@gt$schema)))
    cat("\nPhased genotypes:", has_phased, "\n")

    # --- path info ---
    cat("\nStorage:\n")
    cat("  Path:", object@path, "\n")

    # --- Arrow dataset summary (no collect) ---
    cat("\nGenotype storage (Arrow):\n")
    print(object@gt)

    # --- variants preview ---
    if (!is.null(object@variants) && nrow(object@variants) > 0) {
      cat("\nVariants (first 5 rows):\n")
      print(utils::head(object@variants, 5))
    } else {
      cat("\nVariants: <empty>\n")
    }

    if (!is.null(object@variants) && nrow(object@variants) > 5) {
      cat("  ...", nrow(object@variants) - 5, "more\n")
    }

    # --- INFO preview ---
    if (length(object@info) > 0) {
      cat("\nINFO (first 5):\n")
      print(utils::head(object@info, 5))
    }

    # --- FORMAT preview ---
    if (length(object@info) > 0) {
      cat("\nFORMAT (first 5):\n")
      print(utils::head(object@format, 5))
    }

    # --- samples preview ---
    cat("\nSamples (first 5):\n")
    print(utils::head(object@samples, 5))

    if (length(object@samples) > 5) {
      cat("  ...", length(object@samples) - 5, "more\n")
    }

    invisible(object)
  }
)

#' @title Subset method for VCFArrow
#'
#' @description Subset a VCFArrow object by variants (rows) and samples (columns)
#'
#' @author Tomas Hrbek April 2026
#'
#' @param x A VCFArrow object
#' @param i Variant (row) positions, numeric or logical
#' @param j Column indices: numeric, logical, or sample name character vector
#' @param ... Ignored
#' @param drop Ignored; kept for S4 compatibility
#'
#' @return A new VCFArrow object containing the selected variants and samples
#'
#' @details
#' This function is a method of the VCFArrow S4 class
#' Method to subset by row and column of GT
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f)
#' vcf[1:100, 1:5]
#'
#' @export
#'

setMethod(
  "[",
  signature(x = "VCFArrow", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, ..., drop = FALSE) {

    if (!missing(i)) {
      if (is.character(i)) cli::cli_abort("Row subsetting by character is not supported")
      x <- .vcf_filter_rows(x, x@variants$.row_id[i])
    }

    if (!missing(j)) {
      if (is.character(j) && anyNA(match(j, x@samples))) {
        cli::cli_abort("Some sample names not found")
      }
      x <- .vcf_filter_columns(x, j, f_invar = FALSE, verbose = FALSE)
    }

    x
  }
)
