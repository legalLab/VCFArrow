#' @title vcf2structure
#'
#' @description
#' Converts a VCFArrow object to a Structure or FastStructure infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'structure.str' (character)
#' @param method -> flag for Structure/FastStructure formats, default 'S' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external Structure
#' or FastStructure formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The flag parameter controls whether Structure (flag = 'S') or
#' FastStructure (flag = 'F') formatted output is written out.
#'
#' @examples
#' vcf2structure(vcf_arrow = my_vcf, keep_groups = NULL, out_file = "structure.str", method = "S")
#' vcf2structure(my_vcf, keep_groups, out_file = "structure.str")
#' vcf2structure(my_vcf)
#'

vcf2structure <- function(vcf_arrow, keep_groups = NULL,
                          out_file = "structure.str", method = "S") {

  method <- match.arg(method, c("S", "F"))

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  all_groups <- vcf_arrow@groups

  # keep all groups if no groups specified
  if (is.null(keep_groups)) {
    keep_groups <- unique(all_groups)
  }

  value_fn <- function(a1, a2, REF, ALT) {
    # here: 0 -> 1, 1 -> 2, missing -> -9
    allele1 <- ifelse(is.na(a1), -9L, a1 + 1L)
    allele2 <- ifelse(is.na(a2), -9L, a2 + 1L)
    return(Map(c, allele1, allele2))
  }

  res <- vcf_build_wide(
    vcf_arrow = vcf_arrow,
    keep_groups = keep_groups,
    value_fn = value_fn,
    n_rows_per_cell = 2
  )

  geno_list <- res$data
  samples <- res$samples

  # population structure
  groups_keep <- all_groups[all_groups %in% keep_groups]
  group_members <- split(samples, groups_keep)

  con <- file(out_file, open = "wt")
  on.exit(close(con))

  fill <- if (method == "F") paste(rep(0, 4), collapse = "\t") else NULL

  # iterate populations
  for (group_id in seq_along(group_members)) {
    inds <- group_members[[group_id]]
    n_ind <- length(inds)

    # buffered writing - reduces syscall overhead
    # 2 rows per individual
    buffer <- character(2 * n_ind)
    idx <- 1

    for (ind in inds) {
      mat <- geno_list[[ind]]
      # replace NA just in case
      mat[is.na(mat)] <- -9L
      hap1 <- paste(mat[1, ], collapse = "\t")
      hap2 <- paste(mat[2, ], collapse = "\t")

      if (method == "S") {
        buffer[idx] <- paste(ind, group_id, hap1, sep = "\t")
        buffer[idx + 1] <- paste(ind, group_id, hap2, sep = "\t")
      } else {  # method == "F"
        buffer[idx] <- paste(ind, group_id, fill, hap1, sep = "\t")
        buffer[idx + 1] <- paste(ind, group_id, fill, hap2, sep = "\t")
      }
      idx <- idx + 2
    }
    writeLines(buffer, con)
  }

  invisible(vcf)
}
