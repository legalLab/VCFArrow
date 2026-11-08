#' @title vcf2eigenstrat
#'
#' @description
#' Converts a VCFArrow object to Eigenstrat format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> output file prefix; extensions are added automatically, no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param sex -> sex of the individual, default = U (undefined) (character)
#' @param chrom_code -> how CHROM is written in the .snp file: "auto" (keep
#'   CHROM if it is numeric, otherwise number the contigs 1, 2, 3, ...),
#'   "index" (number the contigs 1, 2, 3, ... in order of first appearance)
#'   or "keep" (verbatim), default "auto" (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external EIGENSTRAT formatted fileset (.geno, .ind and .snp).
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' In the .geno file, each digit is the number of copies of the reference
#' allele, as EIGENSTRAT defines it: the REF allele, given in column 5 of
#' the .snp file (2 = REF homozygote, 1 = heterozygote, 0 = ALT
#' homozygote, 9 = missing).
#' EIGENSOFT needs integer chromosome codes and removes SNPs on chromosome
#' 0, so contig names (e.g. from de novo SNP calling) are numbered 1, 2,
#' 3, ... in order of first appearance (chrom_code = "auto" or "index"),
#' and the correspondence is written to <out_file>.chrommap.  EIGENSOFT
#' reads chromosomes above its numchrom parameter (default 22) as X, Y, ...,
#' so with more than 22 codes set numchrom to at least the largest code in
#' the parameter file; a warning gives the value.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2eigenstrat(vcf, out_file = file.path(tempdir(), "eigenstrat_in"))
#'
#' @export
#'

vcf2eigenstrat <- function(vcf_arrow, out_file, keep_groups = NULL,
                           sex = NULL, chrom_code = c("auto", "index", "keep")) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  chrom_code <- match.arg(chrom_code)
  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  # ── sex vector ─────────────────────────────────────────────────────────────
  #
  # Accepted inputs
  #   NULL            → every retained sample is assigned "U" (unknown)
  #   length-1 string → recycled across all retained samples
  #   length-n_samples vector → used as-is (must match after group filtering)
  #
  # Allowed values follow the EIGENSTRAT convention: "M", "F", or "U".

  if (is.null(sex)) {
    sex_vec <- rep("U", setup$n_samples)
  } else {
    sex <- as.character(sex)
    if (length(sex) == 1L) {
      sex_vec <- rep(sex, setup$n_samples)
    } else if (length(sex) == setup$n_samples) {
      sex_vec <- sex
    } else {
      cli::cli_abort(c(
        "{.arg sex} must be {.code NULL}, a single string, or a vector whose \\
         length equals the number of retained samples.",
        "i" = "Retained samples after group filtering: {setup$n_samples}.",
        "i" = "{.arg sex} length supplied: {length(sex)}."
      ))
    }
  }

  # ── genetic position ───────────────────────────────────────────────────────
  #
  # EIGENSTRAT's .snp column 3 is the genetic position in Morgans.
  # Without a recombination map, dividing physical position (bp) by 1 000 000
  # yields position in Mb, which is a standard proxy for Morgans and is
  # accepted by all major EIGENSTRAT-compatible tools (ADMIXTOOLS, smartpca).
  #
  # setup$variants$POS is the physical position for every retained, filtered
  # variant, already arranged in .row_id order by .vcf_export_setup().

  rel_pos_vec <- setup$variants$POS / 1e6

  # ── .geno ──────────────────────────────────────────────────────────────────

  geno_file <- paste0(out_file, ".geno")
  write_eigenstrat_geno_header_cpp(geno_file)   # create / truncate

  cli::cli_alert_info(
    "Building EIGENSTRAT: {setup$n_var} variant{?s} x {setup$n_samples} sample{?s} \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)"
  )
  cli::cli_alert_info(
    "Writing EIGENSTRAT files..."
  )
  .write_chunks_ordered(.file_tasks(setup$feather_files), .eigenstrat_chunk,
                        shared = .reshape_shared(setup), out_file = geno_file,
                        label = "Writing chunk")

  # ── .ind ───────────────────────────────────────────────────────────────────
  #
  # Three whitespace-separated columns, no header:
  #   sample_name   sex   group_label

  ind_file <- paste0(out_file, ".ind")

  utils::write.table(
    data.frame(sample = setup$samples,
               sex = sex_vec,
               group = setup$samples_groups,
               stringsAsFactors = FALSE),
    file = ind_file,
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  # ── .snp ───────────────────────────────────────────────────────────────────
  #
  # Six whitespace-separated columns, no header:
  #   SNP_ID   CHROM   genetic_pos(Morgans)   physical_pos(bp)   REF   ALT

  snp_file <- paste0(out_file, ".snp")

  # chromosome codes: EIGENSOFT needs positive integers (chromosome 0 is
  # removed), so non-numeric CHROM (e.g. contig names) is numbered 1, 2, ...
  # in order of first appearance
  chrom <- as.character(setup$variants$CHROM)
  stripped <- sub("^(chr|Chr|CHR)[_.-]?", "", chrom)
  numeric <- all(grepl("^[0-9]+$", stripped)) &&
    all(suppressWarnings(as.numeric(stripped)) >= 1)
  chrom_out <- switch(chrom_code,
    keep = chrom,
    auto = if (numeric) as.character(as.numeric(stripped)) else as.character(match(chrom, unique(chrom))),
    index = as.character(match(chrom, unique(chrom))))
  if (chrom_code != "keep") {
    first <- !duplicated(chrom)
    .write_chrom_map(list(recoded = !identical(chrom_out, chrom),
                          chrom_map = data.frame(CHROM = chrom[first], code = chrom_out[first],
                                                 stringsAsFactors = FALSE)),
                     out_file)
    max_code <- suppressWarnings(max(as.numeric(chrom_out), 0))
    if (max_code > 22)
      cli::cli_warn(c(
        "Chromosome codes go up to {max_code}: EIGENSOFT reads chromosomes above \\
         its numchrom parameter (default 22) as X, Y, ...",
        "i" = "Set {.code numchrom: {max_code}} (or higher) in the parameter file."
      ))
  }

  utils::write.table(
    data.frame(ID = setup$loci,
               CHROM = chrom_out,
               # formatted here: write.table() would write e.g. 1e+05
               rel_pos = sprintf("%.8f", rel_pos_vec),
               POS = sprintf("%.0f", setup$variants$POS),
               REF = setup$variants$REF,
               ALT = setup$variants$ALT,
               stringsAsFactors = FALSE),
    file = snp_file,
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  cli::cli_alert_success(
    "EIGENSTRAT fileset written to \\
     {.file {geno_file}}, {.file {ind_file}}, {.file {snp_file}}"
  )

  invisible(vcf_arrow)
}

# Worker side of the chunk loop (see .write_chunks_ordered()): append one chunk to task$part
.eigenstrat_chunk <- function(task, shared) {
  rc <- .read_reshape_chunk(task$fpath, shared)
  if (!is.null(rc)) write_eigenstrat_chunk_cpp(rc$a1, rc$a2, task$part)
  invisible(NULL)
}
