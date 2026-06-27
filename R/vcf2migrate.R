#' @title vcf2migrate
#'
#' @description
#' Converts a VCFArrow object to a MIGRATE-N format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'migrateN_infile.txt' (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external MIGRATE-N formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' The function implements two different SNP formats,
#' new (un)linked SNPs (S) - fixed blocks, and new (un)linked SNPs (N) - chromosome intervals.
#' Alleles are sequence data when SNPs are mapped to a reference.
#' The S option generates blocks of a specific number of SNPs and treats them as linked;
#' this is appropriate if SNPs are extracted without a reference.
#' The N format extracts all SNPs within a chromosome, within the block size and treats them as linked;
#' this is appropriate if SNPs are mapped against a reference.
#' The size of the linked block is determined by the block_size parameter.
#' See https://peterbeerli.com/programs/migrate/distribution_4.x/migratedoc4.x.pdf for format detail.
#'
#' @examples
#' vcf2migrate(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "migrateN_infile.txt")
#' vcf2migrate(vcf_arrow, my_groups, out_file = "migrateN_infile.txt")
#' vcf2migrate(vcf_arrow)
#'
#' @export
#'

# method = "C"  (allele counts, population-level)
#   Format:  n_pops n_loci
#            max_obs pop_name
#            ref_count alt_count   (one line per variant)
#   inc_missing = FALSE drops variants where any sample has missing data.
#
# method = "S"  (nucleotide sequences, individual-level, fixed blocks)
#   Alleles are written as single nucleotide characters (REF or ALT).
#   Haplotypes are concatenated into a sequence string per line.
#   Block header uses (s{block_size}) for full blocks and (n{r}) for the
#   last partial block.  Missing allele: "?".
#
# method = "N"  (nucleotide sequences, individual-level, chromosome intervals)
#   Same output structure as "S" but block sizes are computed from CHROM/POS:
#   each block spans block_size bp within a chromosome.  All block labels
#   use (s{count}) regardless of whether the last interval is partial.

vcf2migrate <- function(vcf_arrow, keep_groups = NULL,
                        out_file = "migrateN_infile.txt",
                        block_size = 100L, method = "S") {

  method <- match.arg(method, c("C", "S", "N"))
  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  if (method == "C") {
    acc <- .accumulate_pops_lowmem(setup, "Migrate-N (C)")
    loci <- setup$loci

    cli::cli_alert_info("Writing Migrate-N ({method}) file...")

    write_migrate_cpp(acc$ref, acc$alt, acc$nobs,
                      setup$group_names, loci, out_file)
  } else {
    acc <- .accumulate_individuals(setup, paste0("Migrate-N (", method, ")"))

    # Compute block_labels: passed to C++ to write the block header line.
    # Method S: fixed window of block_size SNPs.
    #   Full blocks → "(s{block_size})", last partial block → "(n{r})".
    # Method N: window of block_size bp per chromosome.
    #   All block labels use "(s{count})" regardless of size.
    n_var <- setup$n_var
    block_size <- as.integer(block_size)

    if (method == "S") {
      n_full <- n_var %/% block_size
      remainder <- n_var %% block_size
      if (remainder == 0L) {
        block_labels <- rep(paste0("(s", block_size, ")"), n_full)
      } else {
        block_labels <- c(rep(paste0("(s", block_size, ")"), n_full),
                          paste0("(n", remainder, ")"))
      }
    } else {   # method == "N"
      # Compute SNP count per (CHROM, position-interval) pair.
      # Intervals are defined as floor((POS - min_POS_for_chrom) / block_size).
      # dense_rank() is applied within chromosome to ensure sequential labelling
      # when the position arithmetic produces gaps.
      interval_counts <- setup$variants |>
        dplyr::group_by(CHROM) |>
        dplyr::mutate(
          interval = dplyr::dense_rank(
            (POS - min(POS)) %/% block_size
          )
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(CHROM, interval) |>
        dplyr::summarize(count = dplyr::n(), .groups = "drop") |>
        dplyr::pull(count)

      block_labels <- paste0("(s", interval_counts, ")")
    }

    cli::cli_alert_info(
      "Writing Migrate-N ({method}) file..."
    )

    write_migrate_seq_cpp(
      acc$a1, acc$a2,
      setup$variants$REF, setup$variants$ALT,
      setup$samples, setup$group_names, setup$group_sizes,
      block_labels, out_file
    )
  }

  cli::cli_alert_success("Migrate-N file written to {.file {out_file}}")

  invisible(vcf_arrow)
}
