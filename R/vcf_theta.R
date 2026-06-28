#' @title vcf_theta
#'
#' @description
#' Calculates basic Watterson's theta and pi for all samples
#' and for sample groups from VCFArrow format data.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param gt -> VCFArrow gt in long format
#' @param sample -> samples to be included in analysis
#' @param grps -> assignment of individuals to groups
#'
#' @return list of statistics
#'
#' @details
#' This function calculates Watterson's theta and pi for the entire VCFArrow
#' object, and for groups of individuals whose grouping is indicated by
#' the groups slot in the VCFArrow object.
#'
#' @examples
#' vcf_theta(gt_matrix = my_gt_matrix, samples = vcf_arrow@samples, grps = vcf_arrow@groups)
#' vcf_theta(my_gt_matrix, vcf_arrow@samples, vcf_arrow@groups)
#'
#' @export
#'

vcf_theta <- function(vcf_arrow, keep_groups = NULL) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  # Setup + accumulation (shared with every other exporter) ────────────────────
  # .vcf_export_setup() applies the same is_biallelic & !is_indel filter used
  # consistently throughout this package, and restricts to the current
  # @variants$.row_id / @samples
  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  if (setup$n_var == 0L) {
    cli::cli_warn(c(
      "vcf_theta(): zero variants passed the is_biallelic & !is_indel filter \\
       — theta_w/pi will be NaN (written as \"NA\").",
      "i" = "vcf_arrow@variants has {nrow(vcf_arrow@variants)} total rows, but \\
             none satisfy is_biallelic == TRUE & is_indel == FALSE.",
      "i" = "Check: table(vcf_arrow@variants$is_biallelic, useNA = 'always')",
      "i" = "Check: table(vcf_arrow@variants$is_indel, useNA = 'always')"
    ))
  }

  counts <- .accumulate_pops_lowmem(setup, "theta/pi")

  if (setup$n_var > 0L && sum(counts$nobs) == 0L) {
    cli::cli_warn(c(
      "vcf_theta(): {setup$n_var} variants passed the filter, but zero \\
       genotype rows matched during accumulation — theta_w/pi will be NaN \\
       (written as \"NA\").",
      "i" = "This usually indicates a sample-name mismatch. Check:",
      "i" = "  f <- list.files(vcf_arrow@path, pattern = '[.]arrow$', full.names = TRUE)[1]",
      "i" = "  setdiff(vcf_arrow@samples, arrow::read_feather(f)$sample)"
    ))
  }

  n_pops <- setup$n_pops
  n_var <- setup$n_var
  alt <- counts$alt    # n_pops x n_var: ALT allele count per group per variant
  nobs <- counts$nobs   # n_pops x n_var: total observed alleles (= 2 * n_called)

  # Per-(group, variant) statistics — fully vectorised ─────────────────────────
  n_called <- nobs / 2
  p <- ifelse(nobs > 0, alt / nobs, NA_real_)
  valid <- n_called >= 2
  pi_mat <- ifelse(valid, 2 * p * (1 - p), NA_real_)
  seg <- (p > 0 & p < 1) & valid

  # Harmonic numbers for Watterson's a1, CORRECTLY indexed ─────────────────────
  # a1 for N = 2*n_called chromosomes = H_(N-1) = H_(2*n_called - 1).
  # harmonic_numbers[k] = H_k = sum_{i=1}^{k} 1/i  (no NA-offset trick, so no
  # off-by-one risk).
  a1_mat <- matrix(NA_real_, nrow = n_pops, ncol = n_var)
  if (any(valid)) {
    max_n <- max(n_called[valid], na.rm = TRUE)
    harmonic_numbers <- cumsum(1 / seq_len(max(2 * max_n - 1, 1)))
    a1_mat[valid] <- harmonic_numbers[2 * n_called[valid] - 1]
  }

  contrib <- matrix(0, nrow = n_pops, ncol = n_var)
  contrib[seg] <- 1 / a1_mat[seg]

  # Group-level theta_w / pi ───────────────────────────────────────────────────
  theta_w_g <- rowSums(contrib, na.rm = TRUE) / rowSums(valid, na.rm = TRUE)
  pi_g <- rowMeans(pi_mat, na.rm = TRUE)

  #  Global statistics — derived from the SAME accumulation, no re-scan ────────
  # Groups partition the current sample set, so summing the per-group
  # matrices column-wise is exactly equivalent to (and now correctly
  # consistent with) computing across all currently-retained samples.
  alt_total <- colSums(alt)
  nobs_total <- colSums(nobs)
  n_called_total <- nobs_total / 2
  p_total <- ifelse(nobs_total > 0, alt_total / nobs_total, NA_real_)
  valid_total <- n_called_total >= 2
  pi_total_vec <- ifelse(valid_total, 2 * p_total * (1 - p_total), NA_real_)
  seg_total <- (p_total > 0 & p_total < 1) & valid_total

  a1_total <- rep(NA_real_, n_var)
  if (any(valid_total)) {
    max_n_t <- max(n_called_total[valid_total], na.rm = TRUE)
    harmonic_numbers_t <- cumsum(1 / seq_len(max(2 * max_n_t - 1, 1)))
    a1_total[valid_total] <- harmonic_numbers_t[2 * n_called_total[valid_total] - 1]
  }

  contrib_total <- numeric(n_var)
  contrib_total[seg_total] <- 1 / a1_total[seg_total]

  theta_w_total <- sum(contrib_total) / sum(valid_total)
  pi_total <- mean(pi_total_vec, na.rm = TRUE)

  # Per-sample expansion of group-level stats ──────────────────────────────────
  # IMPORTANT: must align row-for-row with vcf_arrow@samples in its ORIGINAL
  # order, because vcf_stats() does cbind(stats, theta$theta_g) — a positional
  # bind, not a join.  .vcf_export_setup() reorders samples into
  # group-contiguous order internally for the matrix accumulation above, which
  # is why this expansion is done separately from vcf_arrow@samples / @groups
  # directly rather than from setup$samples.
  samples_out <- vcf_arrow@samples
  groups_out <- vcf_arrow@groups
  grp_match <- match(groups_out, setup$group_names)

  theta_g <- data.frame(
    sample = samples_out,
    group = groups_out,
    n_ind = setup$group_sizes[grp_match],
    theta_w = theta_w_g[grp_match],
    pi = pi_g[grp_match]
  )
  # Samples whose group was not in keep_groups (grp_match == NA) get NA stats,
  # consistent with the original's match()-based NA propagation.

  list(
    pi = pi_total,
    theta_w = theta_w_total,
    theta_g = theta_g
  )
}
