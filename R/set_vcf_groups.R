#' @title set_vcf_groups
#'
#' @description
#' Associate samples with group info based on "strt" file
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param data_path -> path to data (where strata are located)
#' @param strt -> file with 2+ columns, one id column and one pop column
#'   (tab- or space-separated, with a header line), default 'strata'
#'
#' @return vcf_arrow object with groups slot filled
#'
#' @details
#' This function assigns each sample of a VCFArrow object to a group, based
#' on a strata file.
#' The strt is a file with a header line and 2+ columns, two of them named
#' 'id' (sample name) and 'pop' (group); other columns are ignored.
#' The file is read as tab-separated when its header line contains a tab, so
#' fields may contain spaces and may be empty; otherwise its columns are
#' separated by spaces (any whitespace).
#' All columns are read as text, so ids such as '007' keep their leading
#' zeros.  Ids in the file that are not samples of the object are ignored.
#' An id listed more than once with different groups is an error.
#' Samples not in the file, or with an empty group, get the group NA and a
#' warning lists them; functions that need groups (e.g. the exporters) stop
#' until every sample has one.
#' By default strt is named 'strata' and found in the datapath where the
#' VCF file is found.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- set_vcf_groups(read_vcf(f), data_path = dirname(f),
#'                       strt = "strata")
#' table(vcf@groups)
#'
#' @export
#'

set_vcf_groups <- function(vcf_arrow, data_path, strt = "strata") {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  strata_file <- file.path(data_path, strt)
  if (!file.exists(strata_file))
    cli::cli_abort("Strata file {.file {strata_file}} not found.")

  # tab-separated if the header line has a tab, otherwise separated by
  # whitespace; every column read as text (so ids such as "007" are not
  # turned into numbers), empty fields as missing
  header <- readLines(strata_file, n = 1L, warn = FALSE)
  if (length(header) == 0L)
    cli::cli_abort("Strata file {.file {strata_file}} is empty.")
  strata_df <- utils::read.table(strata_file, header = TRUE,
                                 sep = if (grepl("\t", header, fixed = TRUE)) "\t" else "",
                                 colClasses = "character", strip.white = TRUE,
                                 na.strings = c("", "NA"), quote = "\"",
                                 comment.char = "")

  missing_cols <- setdiff(c("id", "pop"), names(strata_df))
  if (length(missing_cols))
    cli::cli_abort(c(
      "{.file {strata_file}} must have columns {.val id} and {.val pop}.",
      "x" = "Missing: {.val {missing_cols}}",
      "i" = "Columns found: {.val {names(strata_df)}}"
    ))

  # one group per id: repeated identical rows are harmless, conflicting
  # assignments are an error
  strata_df <- unique(strata_df[!is.na(strata_df$id), c("id", "pop")])
  dup <- unique(strata_df$id[duplicated(strata_df$id)])
  if (length(dup))
    cli::cli_abort(c(
      "Some ids are assigned to more than one group in {.file {strata_file}}.",
      "x" = "{.val {dup}}"
    ))

  # group of each sample, in sample order (NA if not in the file)
  groups <- strata_df$pop[match(vcf_arrow@samples, strata_df$id)]

  unassigned <- vcf_arrow@samples[is.na(groups)]
  if (length(unassigned))
    cli::cli_warn(c(
      "{length(unassigned)} sample{?s} not assigned to a group (group set to NA): \\
       {.val {unassigned}}",
      "i" = "Add {?it/them} to {.file {strata_file}}, or remove {?it/them} with \\
             {.code vcf_extract_samples(keep = FALSE)}."
    ))

  vcf_arrow@groups <- groups

  return(vcf_arrow)
}
