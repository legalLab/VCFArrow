#' @title vcf2eigenstrat
#' @description converts vcfR format data to Eigenstrat infiles
#' @description in part based on vcfR2migrate function (vcfR package)
#' @author Tomas Hrbek December 2022
#'
#' @param vcf -> vcfR object
#' @param ind_pop -> population assignment of individuals in vcf (factor)
#' @param keep_pop -> population(s) of interest to include in Eigenstrat infiles (factor)
#' @param inc_missing -> include missing data (logical)
#' @param out_file -> name of file to output (Eigenstrat infiles)
#'
#' @return nothing
#'
#' @details
#' This function converts the vcfR object to a Eigenstrat formatted input files
#' When list of sexes is not provided, lists all individuals as unknown
#' When relative position on chromosome (cM distance or similar) is not provides, list as 0
#' The function will remove indels, and multiallelic loci, and optionally loci with missing data
#'
#' @examples
#' vcf2eigenstrat(vcf = my_vcf, ind_pop = ind_pop, keep_pop = keepers, sex = list_of_sex, rel_pos = marker_cM_map, inc_missing = TRUE, out_file = "Eigenstrat")
#' vcf2eigenstrat(my_vcf, ind_pop, keepers, out_file = "Eigenstrat")
#' vcf2eigenstrat(my_vcf, ind_pop, keepers)
#'

vcf2eigenstrat <- function(vcf_arrow, keep_groups = NULL,
                           out_file = "eigenstrat_infile",
                           sex = "U", rel_pos = 0, inc_missing = TRUE) {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  # ── .geno ──────────────────────────────────────────────────────────────────
  geno_file <- paste0(out_file, ".geno")
  write_eigenstrat_geno_header_cpp(geno_file)   # create / truncate

  if (!inc_missing) {
    # Need complete-variant mask: accumulate nobs, then filter
    cli::cli_alert_info(
      "inc_missing = FALSE: scanning to identify complete variants..."
    )
    nobs_full <- matrix(0L, nrow = setup$n_pops, ncol = setup$n_var)

    cli::cli_progress_bar("Scanning chunk", total = length(setup$feather_files))
    for (fpath in setup$feather_files) {
      chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))
      rc <- .reshape_chunk(chunk, setup)
      if (!is.null(rc)) {
        pc <- .pop_counts_from_chunk(rc, setup)
        nobs_full[, rc$col_idx] <- nobs_full[, rc$col_idx] + pc$nobs
      }
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    complete_mask <- .complete_var_mask(nobs_full, setup)
    complete_ids  <- setup$valid_row_ids[complete_mask]

    # Second pass: write only complete variants
    cli::cli_alert_info(
      "Writing EIGENSTRAT .geno ({sum(complete_mask)} complete variants x \\
       {setup$n_samples} samples)..."
    )
    cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

    for (fpath in setup$feather_files) {
      chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))

      # Filter to complete row_ids only for this pass
      chunk <- chunk[chunk$.row_id %in% complete_ids &
                       chunk$sample  %in% setup$samples, ]
      if (nrow(chunk) == 0L) { cli::cli_progress_update(); next }

      sample_order  <- match(chunk$sample, setup$samples)
      chunk         <- chunk[order(chunk$.row_id, sample_order), ]
      chunk_row_ids <- unique(chunk$.row_id)
      n_cv          <- length(chunk_row_ids)

      a1 <- matrix(chunk$a1, nrow = setup$n_samples, ncol = n_cv)
      a2 <- matrix(chunk$a2, nrow = setup$n_samples, ncol = n_cv)
      write_eigenstrat_chunk_cpp(a1, a2, geno_file)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    # Restrict variant metadata to complete variants
    variants_out <- setup$variants[complete_mask, ]
    loci_out     <- setup$loci[complete_mask]
    rel_pos_out  <- if (length(rel_pos) == 1L) {
      rep_len(rel_pos, sum(complete_mask))
    } else {
      rel_pos[complete_mask]
    }

  } else {
    # Standard chunk-by-chunk path
    cli::cli_alert_info(
      "Writing EIGENSTRAT .geno: {setup$n_var} variants x {setup$n_samples} samples..."
    )
    cli::cli_progress_bar("Writing chunk", total = length(setup$feather_files))

    for (fpath in setup$feather_files) {
      chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"))
      rc <- .reshape_chunk(chunk, setup)
      if (!is.null(rc)) write_eigenstrat_chunk_cpp(rc$a1, rc$a2, geno_file)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    variants_out <- setup$variants
    loci_out     <- setup$loci
    rel_pos_out  <- rep_len(rel_pos, setup$n_var)
  }

  # ── .ind ───────────────────────────────────────────────────────────────────
  ind_file <- paste0(out_file, ".ind")
  sex_vec  <- rep_len(sex, setup$n_samples)

  utils::write.table(
    data.frame(sample = setup$samples,
               sex    = sex_vec,
               group  = setup$samples_groups,
               stringsAsFactors = FALSE),
    file      = ind_file,
    quote     = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE
  )

  # ── .snp ───────────────────────────────────────────────────────────────────
  snp_file <- paste0(out_file, ".snp")

  utils::write.table(
    data.frame(ID      = loci_out,
               CHROM   = variants_out$CHROM,
               rel_pos = rel_pos_out,
               POS     = variants_out$POS,
               REF     = variants_out$REF,
               ALT     = variants_out$ALT,
               stringsAsFactors = FALSE),
    file      = snp_file,
    quote     = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE
  )

  cli::cli_alert_success(
    "EIGENSTRAT files written: \\
     {.file {geno_file}}, {.file {ind_file}}, {.file {snp_file}}"
  )

  invisible(vcf_arrow)
}
