#' @title vcf_stats
#'
#' @description
#' Calculates basic stats of each samples from VCFArrow format data.
#' Includes average read depth per individual, missing data per individual,
#' the inbreeding coefficient F, Watterson's theta and pi.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param res_path -> directory where to write results
#' @param project -> base name of the project file
#' @param theta -> flag to perform theta and pi calculation, default FALSE (Boolean)
#' @param f_ref -> allele frequencies for the inbreeding coefficient F: "all"
#'   (all samples, default, as vcftools and PLINK) or "group" (the sample's
#'   own group) (character)
#'
#' @return table of statistics
#'
#' @details
#' This function calculates average read depth, heterozygosity
#' number of heterozygotes, number of reference and alternative homozygotes,
#' missing data and total number SNPs of each sample in an VCFArrow object.
#' Optionally calls vcf_theta() to get total and group Watterson's theta and pi.
#'
#' The inbreeding coefficient F of each sample is the method-of-moments
#' estimate of vcftools --het and plink --het, F = (O - E) / (N - E). It uses
#' the biallelic SNVs at which the sample is called, at least two samples of
#' the reference set are called and the reference set is polymorphic (N
#' sites; monomorphic sites would not change F); O is the number of homozygous
#' genotypes of the sample at these sites, and E the expected number, the sum
#' of 1 - n / (n - 1) * 2p(1 - p) over the sites, where p is the ALT allele
#' frequency and n the number of called alleles in the reference set. With
#' f_ref = "all" the reference set is all samples; with f_ref = "group" it is
#' the sample's group, which avoids the excess homozygosity that differences
#' in allele frequency between groups cause when they are pooled (Wahlund
#' effect). F is NA for samples without such sites (e.g. a group of one sample
#' with f_ref = "group", or a sample without a group). The table gives F with
#' N, O and E (F_sites, F_hom_obs, F_hom_exp), and for each sample the mean F
#' of its group (F_group_mean; over the group's samples with an F).
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_stats(vcf, res_path = tempdir(), project = "vaillantii",
#'           theta = TRUE)
#'
#' @export
#'

vcf_stats <- function(vcf_arrow, res_path, project, theta = FALSE,
                      f_ref = c("all", "group")) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")
  f_ref <- match.arg(f_ref)

  if (any(is.na(vcf_arrow@groups)))
    cli::cli_alert_warning(
      "Missing group assignments. Use set_vcf_groups() to fill the 'groups' slot."
    )

  samples <- vcf_arrow@samples
  groups <- vcf_arrow@groups
  group_df <- tibble::tibble(sample = samples, group = groups)

  valid_row_ids <- vcf_arrow@variants$.row_id
  n_samples <- length(samples)

  ffiles <- list.files(vcf_arrow@path, pattern = "\\.arrow$", full.names = TRUE)
  if (length(ffiles) == 0L)
    cli::cli_abort("No .arrow files found in {vcf_arrow@path}")
  chunk_nums <- as.integer(stringr::str_extract(basename(ffiles), "\\d+"))
  ffiles <- ffiles[order(chunk_nums)]

  # Per-sample accumulators (O(n_samples) memory — trivial) ────────────────────
  total_loci <- stats::setNames(integer(n_samples), samples)
  called_n <- stats::setNames(integer(n_samples), samples)
  het_n <- stats::setNames(integer(n_samples), samples)
  hom_ref_n <- stats::setNames(integer(n_samples), samples)
  hom_alt_n <- stats::setNames(integer(n_samples), samples)
  dp_n <- stats::setNames(integer(n_samples), samples)
  dp_sum <- stats::setNames(numeric(n_samples), samples)

  cli::cli_alert_info(
    "Computing per-sample stats: {length(valid_row_ids)} variant{?s} x \\
     {n_samples} sample{?s}, reading {length(ffiles)} chunk{?s} directly"
  )

  # inbreeding coefficient F: biallelic SNVs, allele frequencies over all
  # samples (one reference set) or per group (samples without a group skipped)
  f_ok <- vcf_arrow@variants$is_biallelic & !vcf_arrow@variants$is_indel
  f_ok[is.na(f_ok)] <- FALSE
  if (f_ref == "all") {
    f_grp <- rep(1L, n_samples); f_ng <- 1L
  } else {
    f_grp <- as.integer(factor(groups)); f_grp[is.na(f_grp)] <- 0L
    f_ng <- max(c(1L, f_grp))
  }
  idx <- list(pos = .row_id_pos(valid_row_ids), samples = samples,
              f_ok = f_ok, f_grp = f_grp, f_ng = f_ng)

  parts <- .map_chunks(ffiles, .stats_chunk, shared = idx, label = "Scanning chunk")
  f_n <- f_obs <- integer(n_samples); f_exp <- numeric(n_samples)
  seen <- raw(nrow(vcf_arrow@variants))
  split <- FALSE
  for (cnt in parts) {
    # F from per-chunk allele frequencies is valid only if every variant's
    # called genotypes are in a single chunk file
    if (!split && length(cnt$f_vars)) {
      if (any(seen[cnt$f_vars] != as.raw(0))) split <- TRUE
      seen[cnt$f_vars] <- as.raw(1)
    }
    f_n <- f_n + cnt$f$n_sites
    f_obs <- f_obs + cnt$f$obs_hom
    f_exp <- f_exp + cnt$f$exp_hom
    total_loci <- total_loci + cnt$total
    called_n <- called_n + cnt$called
    het_n <- het_n + cnt$het
    hom_ref_n <- hom_ref_n + cnt$hom_ref
    hom_alt_n <- hom_alt_n + cnt$hom_alt
    dp_n <- dp_n + cnt$dp_n
    dp_sum <- dp_sum + cnt$dp_sum
  }
  rm(seen)
  if (split) {
    cli::cli_alert_info("Variants are split across chunk files; computing F in two passes.")
    f <- .f_two_pass(ffiles, idx, nrow(vcf_arrow@variants))
    f_n <- f$n_sites; f_obs <- f$obs_hom; f_exp <- f$exp_hom
  }
  f_val <- (f_obs - f_exp) / (f_n - f_exp)
  f_val[!is.finite(f_val)] <- NA_real_

  # Assemble per-sample stats (equivalent to the original collect() output) ────
  sample_stats <- tibble::tibble(
    sample = samples,
    read_depth = dp_sum[samples] / dp_n[samples],  # NaN if dp_n == 0, matches
    # mean(DP, na.rm=TRUE) on
    # an all-NA input
    homo_ref = hom_ref_n[samples],
    homo_alt = hom_alt_n[samples],
    hetero = het_n[samples],
    missing = total_loci[samples] - called_n[samples],
    non_missing = called_n[samples],
    total_loci = total_loci[samples]
  )
  sample_stats$heterozygosity <-
    sample_stats$hetero / (sample_stats$homo_ref + sample_stats$homo_alt + sample_stats$hetero)
  sample_stats$missing_p <- sample_stats$missing / sample_stats$total_loci

  # since calculation of theta and pi is computationally intensive
  # making it conditional
  if (theta) {
    # theta / pi
    theta_ <- vcf_theta(vcf_arrow)

    # final table
    out <- data.frame(
      sample = samples,
      read_depth = sample_stats$read_depth,
      heterozygosity = sample_stats$heterozygosity,
      heterozygotes = sample_stats$hetero,
      homozygotes = sample_stats$homo_ref + sample_stats$homo_alt,
      homozygotes_ref = sample_stats$homo_ref,
      homozygotes_alt = sample_stats$homo_alt,
      missing_p = sample_stats$missing_p,
      missing_n = sample_stats$missing,
      non_missing = sample_stats$non_missing,
      total_loci = sample_stats$total_loci,
      F = f_val,
      F_sites = f_n,
      F_hom_obs = f_obs,
      F_hom_exp = f_exp,
      theta_total = theta_$theta_w,
      pi_total = theta_$pi
    )
    # explicit group join instead of relying on theta_$theta_g implicitly
    out <- dplyr::left_join(out, theta_$theta_g, by = "sample")
  } else {
    # final table
    out <- data.frame(
      sample = samples,
      read_depth = sample_stats$read_depth,
      heterozygosity = sample_stats$heterozygosity,
      heterozygotes = sample_stats$hetero,
      homozygotes = sample_stats$homo_ref + sample_stats$homo_alt,
      homozygotes_ref = sample_stats$homo_ref,
      homozygotes_alt = sample_stats$homo_alt,
      missing_p = sample_stats$missing_p,
      missing_n = sample_stats$missing,
      non_missing = sample_stats$non_missing,
      total_loci = sample_stats$total_loci,
      F = f_val,
      F_sites = f_n,
      F_hom_obs = f_obs,
      F_hom_exp = f_exp
    )
    # since calculation of theta is skipped, add group information
    out <- dplyr::left_join(out, group_df, by = "sample")
  }

  # mean F of each sample's group (over the group's samples with an F; NA for
  # samples without a group or groups without any F)
  f_mean <- tapply(out$F, out$group, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  out$F_group_mean <- unname(f_mean[match(out$group, names(f_mean))])

  out <- dplyr::arrange(out, group, sample)

  utils::write.table(
    out,
    file = file.path(res_path, paste0(project, "_stats.csv")),
    row.names = FALSE,
    quote = FALSE,
    sep = ","
  )

  invisible(vcf_arrow)
}

# Per-chunk, per-sample counts in one C++ pass (see src/sample_counts.cpp),
# plus the chunk's contributions to F from the chunk's own allele counts and
# the variants with called genotypes in it (to detect split variants)
.stats_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2", "DP"), shared)
  out <- sample_counts_cpp(ch$s, ch$a1, ch$a2, as.numeric(ch$DP), length(shared$samples))
  ac <- f_allele_counts_cpp(ch$s, ch$pos, ch$a1, ch$a2, shared$f_ok, shared$f_grp, shared$f_ng)
  out$f <- f_hom_cpp(ch$s, ch$pos, ch$a1, ch$a2, shared$f_ok, shared$f_grp, ac$alt, ac$nobs,
                     ac$lo, shared$f_ng, length(shared$samples))
  called <- if (ac$n_var) colSums(matrix(ac$nobs, nrow = shared$f_ng)) > 0 else logical(0)
  out$f_vars <- ac$lo - 1L + which(called)
  out
}

# F when some variant's genotypes are split across chunk files (e.g. after
# vcf_bind()): first the allele counts of every variant over all chunks, then
# each chunk's contributions from these counts (in the main session, as the
# counts are groups x variants).
.f_two_pass <- function(ffiles, idx, n_var) {
  ng <- idx$f_ng
  alt <- nobs <- integer(ng * n_var)
  counts <- .map_chunks(ffiles, .f_counts_chunk, shared = idx, label = "Counting alleles")
  for (ac in counts) {
    if (!ac$n_var) next
    cells <- (ac$lo - 1L) * ng + seq_len(ac$n_var * ng)
    alt[cells] <- alt[cells] + ac$alt
    nobs[cells] <- nobs[cells] + ac$nobs
  }
  rm(counts)
  n_samples <- length(idx$samples)
  out <- list(n_sites = integer(n_samples), obs_hom = integer(n_samples),
              exp_hom = numeric(n_samples))
  for (fpath in ffiles) {
    ch <- .read_live_chunk(fpath, c("a1", "a2"), idx)
    r <- f_hom_cpp(ch$s, ch$pos, ch$a1, ch$a2, idx$f_ok, idx$f_grp, alt, nobs, 1L,
                   ng, n_samples)
    for (k in names(out)) out[[k]] <- out[[k]] + r[[k]]
  }
  out
}

# Allele counts per (group, variant) of one chunk, for .f_two_pass()
.f_counts_chunk <- function(fpath, shared) {
  ch <- .read_live_chunk(fpath, c("a1", "a2"), shared)
  f_allele_counts_cpp(ch$s, ch$pos, ch$a1, ch$a2, shared$f_ok, shared$f_grp, shared$f_ng)
}
