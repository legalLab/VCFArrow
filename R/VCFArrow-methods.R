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
#' This function is a method of the VCFArrow S4 class.
#' It prints the numbers of variants and samples, the numbers of biallelic
#' SNVs, multiallelic variants and indels, the samples per group, the
#' percentage of missing genotypes and the number of phased genotypes, the
#' storage location, and the first five variants with their INFO and FORMAT
#' fields and the first five samples.
#' Missing and phased genotypes are counted in one genotype chunk file only,
#' the first that holds current variants and samples, so that printing stays
#' fast.
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

    v <- object@variants
    n_var <- nrow(v)
    ffiles <- list.files(object@path, pattern = "\\.arrow$", full.names = TRUE)

    cat("\nAn object of class \"VCFArrow\"\n")

    # --- basic dimensions ---
    cat("\nDimensions:\n")
    cat("  Variants:", n_var, "\n")
    cat("  Samples: ", length(object@samples), "\n")

    # --- variant types (from @variants, no genotype reading) ---
    if (n_var > 0) {
      cat("\nVariant types:\n")
      cat("  Biallelic SNVs:", sum(v$is_biallelic & !v$is_indel, na.rm = TRUE), "\n")
      cat("  Multiallelic:  ", sum(!v$is_biallelic, na.rm = TRUE), "\n")
      cat("  Indels:        ", sum(v$is_indel, na.rm = TRUE), "\n")
      if (any(is.na(v$POS))) cat("  Without a valid POS:", sum(is.na(v$POS)), "\n")
    }

    # --- groups ---
    if (all(is.na(object@groups))) {
      cat("\nGroups: not set (use set_vcf_groups())\n")
    } else {
      grp <- table(object@groups, useNA = "ifany")
      cat("\nGroups:", length(grp), "\n")
      print(grp)
    }

    # --- missing and phased genotypes: counted in one chunk file only, the
    # first that holds current variants and samples (filtering can leave the
    # first files without any), so that printing stays fast (counting all
    # genotypes would read every chunk file) ---
    if (n_var > 0 && length(ffiles) > 0) {
      live <- list(pos = .row_id_pos(v$.row_id), samples = object@samples)
      for (fpath in .get_sorted_feather_files(object@path)) {
        ch <- .read_live_chunk(fpath, c("a1", "a2", "phased"), live)
        if (length(ch$phased)) break
      }
      if (length(ch$phased)) {
        n_gt <- length(ch$phased)
        cat("\nGenotypes (first chunk with variants):", n_gt, "\n")
        cat("  Missing:", sprintf("%.2f%%", 100 * mean(is.na(ch$a1) | is.na(ch$a2))), "\n")
        cat("  Phased: ", sum(ch$phased, na.rm = TRUE), "\n")
      }
    }

    # --- storage ---
    cat("\nStorage:\n")
    cat("  Path:", object@path, "\n")
    cat("  Genotype chunk files:", length(ffiles), "\n")

    # --- variants preview ---
    if (n_var > 0) {
      cat("\nVariants (first 5 rows):\n")
      v5 <- utils::head(v, 5)
      v5$POS <- format(v5$POS, scientific = FALSE, trim = TRUE)  # no 2.5e+09
      print(v5)
      if (n_var > 5) cat("  ...", n_var - 5, "more\n")
    } else {
      cat("\nVariants: <empty>\n")
    }

    # --- INFO and FORMAT of the same 5 variants (@format is not subset by
    # filters; it is looked up by .row_id) ---
    if (length(object@info) > 0) {
      cat("\nINFO (first 5):\n")
      print(utils::head(object@info, 5))
    }
    if (nrow(object@format) > 0 && n_var > 0) {
      ids <- utils::head(v$.row_id, 5)
      fmt <- object@format$FORMAT[.match_row_id(ids, .row_id_pos(object@format$.row_id))]
      cat("\nFORMAT (first 5):\n")
      print(data.frame(.row_id = ids, FORMAT = fmt), row.names = FALSE)
    }

    # --- samples preview ---
    cat("\nSamples (first 5):\n")
    print(utils::head(object@samples, 5))
    if (length(object@samples) > 5)
      cat("  ...", length(object@samples) - 5, "more\n")

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
