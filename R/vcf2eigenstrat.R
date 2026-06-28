#' @title vcf2eigenstrat
#'
#' @description
#' Converts a VCFArrow object to Eigenstrat format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'eigenstrat_infile' (character)
#' @param sex -> sex of the individual, default = U (undefined) (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external BayesAss formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#'
#' @examples
#' vcf2eigenstrat(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "eigenstrat_infile")
#' vcf2eigenstrat(vcf_arrow, my_groups, out_file = "eigenstrat_infile")
#' vcf2eigenstrat(vcf_arrow)
#'
#' @export
#'

vcf2eigenstrat <- function(vcf_arrow, keep_groups = NULL,
                           out_file = "eigenstrat_infile",
                           sex = NULL) {

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
    "Building EIGENSTRAT: {setup$n_var} variants x {setup$n_samples}  \\
    ({.strong {format(round(2 * setup$n_var * setup$n_samples / 1024^2), big.mark=',')}} MiB raw storage)"
  )
  cli::cli_alert_info(
    "Writing EIGENSTRAT files..."
  )
  cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

  for (fpath in setup$feather_files) {
    chunk <- arrow::read_feather(fpath,
                                 col_select = c(".row_id", "sample", "a1", "a2"))
    rc <- .reshape_chunk(chunk, setup)
    if (!is.null(rc)) {
      write_eigenstrat_chunk_cpp(rc$a1, rc$a2, geno_file)
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

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

  utils::write.table(
    data.frame(ID = setup$loci,
               CHROM = setup$variants$CHROM,
               rel_pos = rel_pos_vec,
               POS = setup$variants$POS,
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
