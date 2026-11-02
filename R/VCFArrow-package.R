#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom methods new show
#' @importFrom stats setNames
#' @importFrom utils read.table write.table
#' @useDynLib VCFArrow, .registration = TRUE
## usethis namespace: end
NULL

# Column names used in dplyr / Arrow data-masking expressions
utils::globalVariables(c(
  ".row_id", "ALT", "CHROM", "DP", "POS", "REF", "block", "count",
  "first_pos", "group", "id", "interval", "key", "is_biallelic", "is_indel",
  "missing_n", "missing_p", "n_alt", "share", "size", "snvs_in_block",
  "take", "total_loci", "u", "variant_index"
))
