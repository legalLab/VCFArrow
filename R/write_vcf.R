#' @title write_vcf
#'
#' @description
#' Write a VCFArrow object to an external VCF file
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of the VCF file to be written to, no default (character)
#' @param gzip -> a flag to GZIP VCF when writing, default TRUE if out_file
#'   ends in ".gz" (Boolean)
#' @param update_info -> recompute the INFO fields AC, AN, AF and NS from the
#'   genotypes written, where present, default TRUE (Boolean)
#'
#' @return Invisibly returns the path of the written file.
#'
#' @details
#' This function writes a VCFArrow object to an external VCF file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' It writes both uncompressed and gz compressed files. Compressing increases
#' writing time be about 50%.
#' For large files, it is recommended to output an uncompressed VCF file,
#' and then compress with GZIP or PIGZ (or BGZIP, which tabix needs to
#' index the file).
#' The sample fields are written as stored; genotypes changed by VCFArrow
#' (e.g. by vcf_filter_adr() or vcf_filter_coverage()) have their GT updated.
#' With update_info = TRUE (the default), the INFO fields AC (count of each
#' ALT allele), AN (number of called alleles), AF (AC / AN) and NS (number of
#' samples with a called genotype) are recomputed from the samples and
#' genotypes written, wherever a variant's INFO has them, so they stay
#' correct after samples or genotypes were removed or changed; other INFO
#' fields are written unchanged.  Haploid genotypes, which read_vcf() stores
#' as homozygous, count as two alleles.
#' Variants without a valid position (POS) are not written, with a warning.
#' A warning is also given if the variants are not sorted by position (e.g.
#' after vcf_bind(mode = "union")), since indexing tools such as tabix need
#' sorted files; `bcftools sort` sorts them.
#' A '##source=VCFArrow' line with the package version is added to the header.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' write_vcf(vcf, out_file = tempfile(fileext = ".vcf"))
#' write_vcf(vcf, out_file = tempfile(fileext = ".vcf.gz"))
#'
#' @export
#'

write_vcf <- function(vcf_arrow, out_file, gzip = grepl("\\.gz$", out_file),
                      update_info = TRUE) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  # write header (with gzip, as its own gzip member: the chunks are appended
  # as further members, which together form a valid multi-member gzip file),
  # with a ##source line naming this package unless it has one already
  header <- vcf_arrow@header
  source_line <- paste0("##source=VCFArrow ", utils::packageVersion("VCFArrow"))
  meta_lines <- header[-length(header)]
  if (!source_line %in% meta_lines) meta_lines <- c(meta_lines, source_line)
  header <- c(meta_lines, paste(
    c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
      vcf_arrow@samples),
    collapse = "\t"
  ))
  con <- if (gzip) gzfile(out_file, "w") else file(out_file, "w")
  writeLines(header, con)
  close(con)

  variants <- vcf_arrow@variants

  # a VCF line needs a position: variants without one are not written
  no_pos <- is.na(variants$POS)
  if (any(no_pos)) {
    cli::cli_warn(
      "{sum(no_pos)} variant{?s} without a valid POS not written."
    )
    variants <- variants[!no_pos, , drop = FALSE]
  }

  # variants are written in .row_id order: warn if that is not sorted by
  # position (each CHROM in one block, POS not decreasing within it)
  o <- order(variants$.row_id)
  chrom <- variants$CHROM[o]
  pos <- variants$POS[o]
  if (length(o) > 1L) {
    same <- chrom[-1L] == chrom[-length(chrom)]
    chrom_runs <- chrom[c(TRUE, !same)]
    if (anyDuplicated(chrom_runs) > 0L || any(pos[-1L][same] < pos[-length(pos)][same]))
      cli::cli_warn(c(
        "The variants are not sorted by position (CHROM, POS), as after \\
         {.code vcf_bind(mode = \"union\")}.",
        "i" = "Sort the file (e.g. {.code bcftools sort}) before indexing it."
      ))
  }

  # locate and sort feather files by chunk number
  # natural numeric sort avoids chunk_10 < chunk_2 with plain sort()
  feather_files <- .get_sorted_feather_files(vcf_arrow@path)

  # Each chunk's task carries the metadata of the variants in its .row_id
  # range (CHROM ... INFO, and FORMAT from the per-variant lookup), so the
  # chunks can be written in parallel without sending all metadata to every
  # worker.
  fmt_lookup <- vcf_arrow@format
  info <- vcf_arrow@info
  if (any(no_pos)) info <- info[!no_pos]
  meta <- data.frame(
    .row_id = variants$.row_id,
    CHROM = variants$CHROM, POS = variants$POS, ID = variants$ID,
    REF = variants$REF, ALT = variants$ALT, QUAL = variants$QUAL,
    FILTER = variants$FILTER, INFO = info,
    FORMAT = fmt_lookup$FORMAT[.match_row_id(variants$.row_id,
                                             .row_id_pos(fmt_lookup$.row_id))],
    stringsAsFactors = FALSE
  )
  tasks <- lapply(feather_files, function(f) {
    id <- as.vector(arrow::read_feather(f, col_select = ".row_id",
                                        as_data_frame = FALSE)$.row_id)
    in_range <- if (length(id)) meta$.row_id >= min(id) & meta$.row_id <= max(id)
                else logical(nrow(meta))
    list(fpath = f, meta = meta[in_range, , drop = FALSE])
  })
  meta <- NULL

  # chunk size message
  cli::cli_alert_info("VCF is being written in {length(feather_files)} chunk{?s}")

  .write_chunks_ordered(tasks, .write_vcf_chunk,
                        shared = list(samples = vcf_arrow@samples, gzip = gzip,
                                      update_info = update_info),
                        out_file = out_file, label = "Writing VCF chunk")

  cli::cli_alert_success("VCFArrow object successfully written to {.file {out_file}}")

  invisible(out_file)
}

# Worker side of write_vcf() (see .write_chunks_ordered()): append the VCF
# lines of one chunk to task$part.  task$meta: metadata of the variants in the
# chunk's .row_id range.

.write_vcf_chunk <- function(task, shared) {
  samples <- shared$samples
  meta <- task$meta
  var_pos <- .row_id_pos(meta$.row_id)

  # only the per-sample FORMAT strings are written; the other genotype
  # columns are derived from them at read time (a1/a2 are read to update
  # INFO)
  cols <- c(".row_id", "sample", "fmt", if (shared$update_info) c("a1", "a2"))
  chunk <- arrow::read_feather(task$fpath, col_select = cols, as_data_frame = FALSE)

  # remove filtered variants and samples from chunks
  # filtering here mirrors the .reshape_chunk() logic used in all exporters
  sample_order <- .sample_index(chunk, samples)
  row_id <- as.vector(chunk$.row_id)
  keep <- !is.na(.match_row_id(row_id, var_pos)) & !is.na(sample_order)
  # nothing to write if no variant of this chunk is retained (e.g. all were
  # filtered out); this also avoids integer(0)[TRUE], which is NA, below
  if (!any(keep)) return(invisible(NULL))

  # sort so rows are: variant 1 sample 1, variant 1 sample 2, ..., variant 2 sample 1, ...
  # in practice read_vcf() writes them this way already, but sort defensively
  # (rows: positions in the chunk of the kept rows, in that order)
  rows <- which(keep)
  rows <- rows[order(row_id[rows], sample_order[rows])]
  row_id <- row_id[rows]

  # one entry per variant in chunk (row_id is sorted)
  first <- c(TRUE, row_id[-1L] != row_id[-length(row_id)])
  row_ids <- row_id[first]

  # variant metadata for this chunk (indexed by .row_id)
  vi <- .match_row_id(row_ids, var_pos)
  info <- meta$INFO[vi]
  if (shared$update_info)
    info <- .update_info(info, meta$ALT[vi], cumsum(first),
                         as.vector(chunk$a1)[rows], as.vector(chunk$a2)[rows])

  # per-sample FORMAT strings in that order (flat, row-major:
  # [var1_s1, var1_s2, ..., var2_s1, ...]), selected in Arrow and handed to
  # the C++ writer as an Arrow array, so they never become R strings
  fmt <- chunk$fmt$Take(arrow::Array$create(rows - 1L))
  fmt <- if (fmt$num_chunks == 1L) fmt$chunk(0L)
         else do.call(arrow::concat_arrays, fmt$chunks)
  fmt_c <- arrow_c_alloc_cpp()
  fmt$export_to_c(fmt_c$array, fmt_c$schema)
  fmt <- chunk <- NULL

  write_vcf_chunk_cpp(
    output_file = task$part,
    chrom = meta$CHROM[vi],
    pos = as.numeric(meta$POS[vi]),
    id = meta$ID[vi],
    ref = meta$REF[vi],
    alt = meta$ALT[vi],
    qual = meta$QUAL[vi],
    filter_col = meta$FILTER[vi],
    info = info,
    format_col = meta$FORMAT[vi],
    fmt_array = fmt_c$array,
    fmt_schema = fmt_c$schema,
    n_samples = length(samples),
    gzip = shared$gzip
  )
  invisible(NULL)
}

# INFO strings `info` (one per variant) with AC, AN, AF and NS replaced by
# values computed from the genotype rows a1/a2, where the keys are present
# (other keys are left unchanged; keys are not added).  var: the variant
# (1..length(info)) of each genotype row; alt: the ALT field of each variant.
#   NS  samples with a called genotype (both alleles)
#   AN  called alleles (2 * NS)
#   AC  copies of each ALT allele, comma-separated (allele codes 1..n_alt)
#   AF  AC / AN, "." when AN is 0

.update_info <- function(info, alt, var, a1, a2) {
  has_key <- function(key) grepl(paste0("(^|;)", key, "="), info)
  todo <- has_key("AC") | has_key("AN") | has_key("AF") | has_key("NS")
  if (!any(todo)) return(info)

  nv <- length(info)
  called <- !is.na(a1) & !is.na(a2)
  ns <- tabulate(var[called], nv)
  an <- 2 * ns
  n_alt <- stringr::str_count(alt, ",") + 1L
  k <- max(n_alt)
  # counts[v, j]: copies of ALT allele j at variant v
  alleles <- c(a1[called], a2[called])
  avar <- c(var[called], var[called])
  ok <- alleles >= 1L & alleles <= k
  counts <- matrix(tabulate((avar[ok] - 1L) * k + alleles[ok], nv * k),
                   nrow = nv, ncol = k, byrow = TRUE)

  fmt_af <- function(x) ifelse(is.finite(x), sprintf("%g", x), ".")
  if (k == 1L) {
    ac <- as.character(counts[, 1L])
    af <- ifelse(an > 0, fmt_af(counts[, 1L] / an), ".")
  } else {
    ac <- af <- character(nv)
    for (v in which(todo)) {
      cv <- counts[v, seq_len(n_alt[v])]
      ac[v] <- paste(cv, collapse = ",")
      af[v] <- if (an[v] > 0) paste(fmt_af(cv / an[v]), collapse = ",")
               else paste(rep(".", n_alt[v]), collapse = ",")
    }
  }

  replace_key <- function(info, key, value)
    stringr::str_replace(info, paste0("(^|;)", key, "=[^;]*"), paste0("\\1", key, "=", value))
  info[todo] <- replace_key(info[todo], "AC", ac[todo])
  info[todo] <- replace_key(info[todo], "AN", an[todo])
  info[todo] <- replace_key(info[todo], "AF", af[todo])
  info[todo] <- replace_key(info[todo], "NS", ns[todo])
  info
}
