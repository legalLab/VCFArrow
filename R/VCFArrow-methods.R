#' @title Show method for VCFArrow
#'
#' @description Print a summary of a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param object A VCFArrow object
#'
#' @details
#' This function is a method of the VCFArrow S4 class
#' Method to show object content summary
#' Method to subset by row and column of GT
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
#' @param i Row indices (numeric or logical)
#' @param j Column indices: numeric, logical, or sample name character vector
#' @param drop Ignored; kept for S4 compatibility
#'
#' @return A new VCFArrow object containing the selected variants and samples
#'
#' @details
#' This function is a method of the VCFArrow S4 class
#' Method to subset by row and column of GT
#'
#' @export
#'

setMethod(
  "[",
  signature(x = "VCFArrow", i = "ANY", j = "ANY", drop = "ANY"),
  function(x, i, j, ..., drop = FALSE) {

    # --- Handle missing indices ---
    if (missing(i)) i <- seq_len(nrow(x@variants))
    if (missing(j)) j <- seq_along(x@samples)

    # --- Normalize row index ---
    if (is.logical(i)) i <- which(i)
    if (is.character(i)) stop("Row subsetting by character not supported")

    # --- Normalize column index ---
    if (is.logical(j)) j <- which(j)

    if (is.character(j)) {
      j <- match(j, x@samples)
      if (any(is.na(j))) {
        stop("Some sample names not found")
      }
    }

    # --- Subset samples ---
    new_samples <- x@samples[j]

    # --- Build Arrow column selection ---
    col_names <- unlist(lapply(new_samples, function(s) {
      c(
        paste0(s, "_a1"),
        paste0(s, "_a2"),
        paste0(s, "_phased")
      )
    }))

    # always keep variant_index
    col_names <- c(col_names, "variant_index")

    # --- Subset Arrow dataset ---
    new_gt <- x@gt %>%
      dplyr::select(dplyr::all_of(col_names))

    # row filtering (lazy)
    if (isTRUE(all(diff(i) == 1))) {
      new_gt <- new_gt %>%
        dplyr::filter(
          variant_index >= min(i),
          variant_index <= max(i)
        )
    }

    # --- Subset metadata ---
    new_variants <- x@variants[i, , drop = FALSE]
    new_info <- x@info[i]

    # --- IMPORTANT: reindex variant_index ---
    new_index <- seq_along(i)

    new_gt <- new_gt %>%
      dplyr::mutate(variant_index = new_index)

    # --- Create new object (SAFE) ---
    new_obj <- .new_vcfarrow(
      header = x@header,
      info = new_info,
      format = x@format,
      variants = new_variants,
      gt = new_gt,
      samples = new_samples,
      groups = x@groups,
      path = x@path  # shared dataset
    )

    return(new_obj)
  }
)
