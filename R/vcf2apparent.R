#' @title vcf2apparent
#'
#' @description
#' Converts a VCFArrow object to an Apparent format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param key -> relationship type (All, Pa, Mo, Fa, Off), default All (character)
#' @param out_file -> name of file to output, default 'apparent_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external SmartSNP formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Possible relationships defined by the parameter 'kee' are All, Pa, Mo, Fa, Off.
#'
#' @examples
#' vcf2apparent(vcf_arrow = my_vcf, keep_groups = my_groups, key = my_key, out_file = "apparent_infile.txt")
#' vcf2apparent(vcf_arrow, my_groups, my_key, out_file = "apparent_infile.txt")
#' vcf2apparent(vcf_arrow)
#'

vcf2apparent <- function(vcf_arrow, keep_groups = NULL, key = "All",
                         out_file = "apparent_infile.txt") {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  all_groups <- vcf_arrow@groups

  # keep all groups if no groups specified
  if (is.null(keep_groups)) {
    keep_groups <- unique(all_groups)
  }

  # define encoder
  value_fn <- function(a1, a2, REF, ALT) {
    allele1 <- ifelse(a1 == 0L, REF, ifelse(a1 == 1L, ALT, "?"))
    allele2 <- ifelse(a2 == 0L, REF, ifelse(a2 == 1L, ALT, "?"))
    # canonical ordering (in principle needed)
    swap <- !is.na(code1) & !is.na(code2) & code1 > code2
    tmp <- allele1
    allele1[swap] <- allele2[swap]
    allele2[swap] <- tmp[swap]
    return(paste0(allele1, "/", allele2))
  }

  res <- vcf_build_wide(
    vcf_arrow = vcf_arrow,
    keep_groups = keep_groups,
    value_fn = value_fn,
    n_rows_per_cell = 1
  )

  geno_list <- res$data
  samples <- res$samples

  con <- file(out_file, open = "wt")
  on.exit(close(con))

  # buffered writing - reduces syscall overhead
  buffer <- character(length(samples))
  for (i in seq_along(samples)) {
    buffer[i] <- paste(key, paste(geno_list[[i]], collapse = "\t"), sep = "\t")
  }
  writeLines(buffer, con)

  invisible(vcf_arrow)
}
