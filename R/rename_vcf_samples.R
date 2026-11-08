#' @title rename_vcf_samples
#'
#' @description
#' Rename samples of a VCFArrow object
#'
#' @author Tomas Hrbek November 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param new_names -> the new sample names: a named character vector
#'   (c(old = "new", ...)), a data frame, or the path of a tab- or
#'   space-separated file with a header line; the data frame or file has the
#'   columns id (current name) and new_id (new name), and optionally pop
#'   (group) (character, data.frame or file path)
#' @param update_groups -> also update the group assignments, from the pop
#'   column of new_names, default FALSE (keep the current group of each
#'   sample) (Boolean)
#'
#' @return VCFArrow object with the samples renamed
#'
#' @details
#' This function renames samples of a VCFArrow object, returning a new
#' VCFArrow object.
#' Samples not listed in new_names keep their names; names in new_names that
#' are not samples of the object are ignored, with a warning.  The new names
#' must be unique (also with respect to the samples not renamed) and must not
#' be empty or contain tabs or line breaks; a sample listed twice with
#' different new names is an error.  Names can also be exchanged (e.g.
#' c(A = "B", B = "A")).
#' By default (update_groups = FALSE) each sample keeps its group.  With
#' update_groups = TRUE, the samples listed get the group in the pop column
#' of new_names (a listed sample with an empty pop keeps its group); this
#' needs a data frame or file with a pop column.  A file is read as in
#' set_vcf_groups(): tab-separated if its header line has a tab (fields may
#' then contain spaces), otherwise space-separated, all columns as text.
#' Sample names are stored with every genotype, so the genotype files are
#' rewritten (in parallel with workers, see vcf_set_workers()); the genotypes
#' themselves are unchanged.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2 <- rename_vcf_samples(vcf, c(Pv120 = "WA_Pv120", Pv13 = "WA_Pv13"))
#' head(vcf2@samples)
#' # rename and move a sample to another group
#' new <- data.frame(id = "Pv120", new_id = "OG_Pv120", pop = "OG")
#' vcf3 <- rename_vcf_samples(vcf, new, update_groups = TRUE)
#'
#' @export
#'

rename_vcf_samples <- function(vcf_arrow, new_names, update_groups = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  if (!isTRUE(update_groups) && !isFALSE(update_groups))
    cli::cli_abort("{.arg update_groups} must be TRUE or FALSE.")

  map <- .rename_table(new_names)
  samples <- vcf_arrow@samples

  # one new name per sample: repeated identical rows are harmless
  map <- unique(map[!is.na(map$id), , drop = FALSE])
  dup <- unique(map$id[duplicated(map$id)])
  if (length(dup))
    cli::cli_abort(c(
      "Some samples are given more than one new name (or group).",
      "x" = "{.val {dup}}"
    ))
  bad <- is.na(map$new_id) | map$new_id == "" | grepl("[\t\r\n]", map$new_id)
  if (any(bad))
    cli::cli_abort(c(
      "New names must not be empty or contain tabs or line breaks.",
      "x" = "For sample{?s} {.val {map$id[bad]}}"
    ))

  unknown <- !map$id %in% samples
  if (any(unknown))
    cli::cli_warn(c(
      "{sum(unknown)} name{?s} in {.arg new_names} {?is/are} not a sample of the \\
       object and ignored: {.val {map$id[unknown]}}"
    ))
  map <- map[!unknown, , drop = FALSE]

  # new sample names, unique over all samples
  pos <- match(map$id, samples)
  new_samples <- samples
  new_samples[pos] <- map$new_id
  dup <- unique(new_samples[duplicated(new_samples)])
  if (length(dup))
    cli::cli_abort(c(
      "Sample names must be unique after renaming.",
      "x" = "Duplicated: {.val {dup}}"
    ))

  # groups: kept, or taken from the pop column for the samples listed
  groups <- vcf_arrow@groups
  if (update_groups) {
    if (!"pop" %in% names(map))
      cli::cli_abort(c(
        "{.code update_groups = TRUE} needs a {.val pop} column in {.arg new_names}.",
        "i" = "Give a data frame or file with columns id, new_id and pop."
      ))
    has_pop <- !is.na(map$pop)
    groups[pos[has_pop]] <- map$pop[has_pop]
  }

  renamed <- map$id != map$new_id
  if (!any(renamed)) {
    cli::cli_alert_info("No sample names changed.")
    vcf_arrow@groups <- groups
    return(vcf_arrow)
  }

  # Sample names are stored with every genotype: rewrite the sample column of
  # every chunk file to a new directory (renamed samples only; rows of other
  # samples, including ones filtered out, keep their names)
  ffiles <- .get_sorted_feather_files(vcf_arrow@path)
  tmp_dir <- tempfile("arrow_vcf_rename_")
  dir.create(tmp_dir)
  tasks <- lapply(seq_along(ffiles), function(i) list(
    fpath = ffiles[[i]], out = file.path(tmp_dir, paste0("chunk_", i, ".arrow"))
  ))
  .map_chunks(tasks, .rename_chunk,
              shared = list(old = map$id[renamed], new = map$new_id[renamed]),
              label = "Rewriting chunk")

  vcf_arrow <- .new_vcfarrow(
    header = vcf_arrow@header,
    info = vcf_arrow@info,
    format = vcf_arrow@format,
    variants = vcf_arrow@variants,
    gt = suppressWarnings(arrow::open_dataset(tmp_dir, format = "feather")),
    samples = new_samples,
    groups = groups,
    path = tmp_dir,
    invariant_removed = vcf_arrow@invariant_removed,
    loci = vcf_arrow@loci
  )

  n_group <- if (update_groups) sum(!is.na(map$pop)) else 0L
  cli::cli_alert_info(
    "Renamed {sum(renamed)} sample{?s}; {if (update_groups) paste(n_group, 'group assignment(s) updated') else 'group assignments kept'}."
  )

  return(vcf_arrow)
}

# Table of renames (columns id, new_id and, if given, pop; all character)
# from a named character vector, a data frame or a file path.
.rename_table <- function(new_names) {
  if (is.data.frame(new_names)) {
    map <- as.data.frame(lapply(new_names, as.character), stringsAsFactors = FALSE)
  } else if (is.character(new_names) && !is.null(names(new_names))) {
    map <- data.frame(id = names(new_names), new_id = unname(new_names),
                      stringsAsFactors = FALSE)
  } else if (is.character(new_names) && length(new_names) == 1L) {
    map <- .read_sample_table(new_names, "File")
  } else {
    cli::cli_abort(c(
      "{.arg new_names} must be a named character vector, a data frame or a file path.",
      "i" = "E.g. {.code c(old_name = \"new_name\")}, or columns id and new_id (and pop)."
    ))
  }
  missing_cols <- setdiff(c("id", "new_id"), names(map))
  if (length(missing_cols))
    cli::cli_abort(c(
      "{.arg new_names} must have columns {.val id} and {.val new_id}.",
      "x" = "Missing: {.val {missing_cols}}",
      "i" = "Columns found: {.val {names(map)}}"
    ))
  map[intersect(c("id", "new_id", "pop"), names(map))]
}

# Worker side of rename_vcf_samples(): replace the names shared$old by
# shared$new in the sample column of one chunk and write it to task$out.
# The chunk stays an Arrow Table: names are matched and replaced in Arrow,
# without converting the per-genotype strings to R.
.rename_chunk <- function(task, shared) {
  chunk <- arrow::read_feather(task$fpath, as_data_frame = FALSE)
  sample <- chunk$sample
  sample <- if (sample$num_chunks == 1L) sample$chunk(0L)
            else do.call(arrow::concat_arrays, sample$chunks)
  idx <- as.vector(arrow::call_function(
    "index_in", sample,
    options = list(value_set = arrow::Array$create(shared$old), skip_nulls = FALSE)
  ))
  if (any(!is.na(idx))) {
    new <- arrow::call_function("take", arrow::Array$create(shared$new),
                                arrow::Array$create(idx))
    chunk[["sample"]] <- arrow::call_function("coalesce", new, sample)
  }
  arrow::write_feather(chunk, task$out)
  chunk <- NULL
  .release_arrow_memory()
  invisible(TRUE)
}
