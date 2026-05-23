#' @title vcf2smartsnp
#'
#' @description
#' Converts a VCFArrow object to smartsnp tabular format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'smartsnp_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2smartsnp(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "smartsnp_infile.txt")
#' vcf2smartsnp(vcf_arrow, my_groups, out_file = "smartsnp_infile.txt")
#' vcf2smartsnp(vcf_arrow)
#'

vcf2smartsnp <- function(vcf_arrow, keep_groups = NULL,
                         out_file = "smartsnp_infile.txt") {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  all_groups <- vcf_arrow@groups

  # keep all groups if no groups specified
  if (is.null(keep_groups)) {
    keep_groups <- unique(all_groups)
  }

  # define encoder
  # 0 = hom-ref, 1 = het, 2 = hom-alt, 9 = missing
  value_fn <- function(a1, a2, REF, ALT) {
    geno <- rep(9L, length(a1))
    ok <- !is.na(a1) & !is.na(a2)
    geno[ok & a1 == 0L & a2 == 0L] <- 0L
    geno[ok & a1 == 1L & a2 == 1L] <- 2L
    geno[ok & a1 != a2] <- 1L
    return(as.character(geno))
  }

  res <- vcf_build_wide(
    vcf_arrow = vcf_arrow,
    keep_groups = keep_groups,
    value_fn = value_fn,
    n_rows_per_cell = 1
  )

  con <- file(out_file, open = "wt")
  on.exit(close(con))

  writeLines(paste(res$samples, collapse = " "), con)

  # geno_list is (sample → n_var vector); transpose to (variant → all samples)
  mat <- do.call(rbind, res$data) # n_samples × n_variants
  mat[is.na(mat)] <- "9"

  # buffered writing - reduces syscall overhead
  buffer <- character(ncol(mat))
  for (i in seq_len(ncol(mat))) {
    buffer[i] <- paste(mat[, i], collapse = " ")
  }
  writeLines(buffer, con)

  invisible(vcf_arrow)
}
