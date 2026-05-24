#' @title vcf2related
#'
#' @description
#' Converts a VCFArrow object to Related format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'related_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external Related formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2related(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "related_infile.txt")
#' vcf2related(vcf_arrow, my_groups, out_file = "related_infile.txt")
#' vcf2related(vcf_arrow)
#'

vcf2related <- function(vcf_arrow, keep_groups = NULL,
                        out_file = "related_infile.txt") {

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
    # nucleotide → numeric code
    nuc_code <- c(A = "01", C = "02", G = "03", T = "04")
    allele1 <- ifelse(a1 == 0L, REF, ifelse(a1 == 1L, ALT, NA_integer_))
    allele2 <- ifelse(a2 == 0L, REF, ifelse(a2 == 1L, ALT, NA_integer_))
    # numeric encoding
    code1 <- nuc_code[allele1]
    code2 <- nuc_code[allele2]
    # canonical ordering (in principle needed)
    swap <- !is.na(code1) & !is.na(code2) & code1 > code2
    tmp <- code1
    code1[swap] <- code2[swap]
    code2[swap] <- tmp[swap]
    code1[is.na(code1)] <- "NA"
    code2[is.na(code2)] <- "NA"
    # return list-of-length-2 vectors
    return(Map(c, code1, code2))
  }

  res <- vcf_build_wide(
    vcf_arrow = vcf_arrow,
    keep_groups = keep_groups,
    value_fn = value_fn,
    n_rows_per_cell = 2
  )

  geno_list <- res$data
  samples <- res$samples

  con <- file(out_file, open = "wt")
  on.exit(close(con))

  # buffered writing - reduces syscall overhead
  n_ind <- length(samples)
  buffer <- character(n_ind)
  for (i in seq_len(n_ind)) {
    res <- geno_list[[i]]
    # replace NA with "NA"
    res[is.na(res)] <- "NA"
    # interleave rows: a1 a2 a1 a2 ...
    vec <- as.vector(rbind(res[1, ], res[2, ]))
    buffer[i] <- paste(vec, collapse = "\t")
  }
  writeLines(buffer, con)

  invisible(vcf_arrow)
}
