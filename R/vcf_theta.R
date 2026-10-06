#' @title vcf_theta
#'
#' @description
#' Calculates basic Watterson's theta and pi for all samples
#' and for sample groups from VCFArrow format data.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#'
#' @return list of statistics
#'
#' @details
#' This function calculates Watterson's theta and pi for the entire VCFArrow
#' object, and for groups of individuals whose grouping is indicated by
#' the groups slot in the VCFArrow object.
#' The statistics are accumulated chunk by chunk as running sums per group,
#' so memory use does not grow with the number of variants, and chunks are
#' processed in parallel when workers are set with vcf_set_workers().
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf_theta(vcf)
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
       - theta_w/pi will be NaN (written as \"NA\").",
      "i" = "vcf_arrow@variants has {nrow(vcf_arrow@variants)} total rows, but \\
             none satisfy is_biallelic == TRUE & is_indel == FALSE.",
      "i" = "Check: table(vcf_arrow@variants$is_biallelic, useNA = 'always')",
      "i" = "Check: table(vcf_arrow@variants$is_indel, useNA = 'always')"
    ))
  }

  # Running sums, chunk by chunk ──────────────────────────────────────────────
  # theta_w and pi are means over the valid sites of per-site values (see
  # .theta_sums()), so per group only three sums are needed: the Watterson
  # contributions, the pi values and the number of valid sites. Each chunk's
  # sums are added as it is read (in parallel with workers), so memory does
  # not grow with the number of variants.
  sums <- .theta_stream(setup)
  if (is.null(sums)) {
    # a variant's genotypes were split across chunk files: its counts must be
    # added up before its statistics can be computed, so accumulate all counts
    # first and compute the sums from them in blocks of variants
    cli::cli_alert_info("Variants are split across chunk files; accumulating counts first.")
    sums <- .theta_from_counts(.accumulate_pops_lowmem(setup, "theta/pi"))
  }

  if (setup$n_var > 0L && sums$nobs == 0) {
    cli::cli_warn(c(
      "vcf_theta(): {setup$n_var} variants passed the filter, but zero \\
       genotype rows matched during accumulation - theta_w/pi will be NaN \\
       (written as \"NA\").",
      "i" = "This usually indicates a sample-name mismatch. Check:",
      "i" = "  f <- list.files(vcf_arrow@path, pattern = '[.]arrow$', full.names = TRUE)[1]",
      "i" = "  setdiff(vcf_arrow@samples, arrow::read_feather(f)$sample)"
    ))
  }

  # group-level and global theta_w / pi (NaN when a group has no valid site)
  theta_w_g <- sums$group[, "contrib"] / sums$group[, "n_valid"]
  pi_g <- sums$group[, "pi"] / sums$group[, "n_valid"]
  theta_w_total <- sums$total[["contrib"]] / sums$total[["n_valid"]]
  pi_total <- sums$total[["pi"]] / sums$total[["n_valid"]]

  # Per-sample expansion of group-level stats ──────────────────────────────────
  # IMPORTANT: must align row-for-row with vcf_arrow@samples in its ORIGINAL
  # order, because vcf_stats() does cbind(stats, theta$theta_g) — a positional
  # bind, not a join.  .vcf_export_setup() reorders samples into
  # group-contiguous order internally for the counting above, which
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


# ── Watterson's theta and pi as running sums ─────────────────────────────────
#
# For one group and one variant with nobs observed alleles (= 2 * n_called)
# and alt ALT alleles:
#   valid   n_called >= 2
#   p       alt / nobs
#   pi      2 * p * (1 - p)
#   contrib 1 / a1 if the site segregates in the group (0 < p < 1), else 0,
#           with a1 = H(2 * n_called - 1), the harmonic number for the
#           2 * n_called sampled chromosomes
# theta_w = sum(contrib) / n_valid and pi = sum(pi) / n_valid, both over the
# valid sites. The sums of contrib, pi and valid are therefore all that is
# kept, per group and for all samples together (counts summed over groups).

# Sums of contrib, pi and valid per row of the allele-count matrices `alt`
# and `nobs` (groups x variants); `harmonic` = H(1), H(2), ... (H(k) = sum of
# 1/i for i = 1..k), long enough for the largest sample size.
.theta_sums <- function(alt, nobs, harmonic) {
  n_called <- nobs / 2
  valid <- n_called >= 2
  p <- alt / nobs                       # NaN where nobs == 0 (never valid)
  pi <- 2 * p * (1 - p)
  pi[!valid] <- 0
  seg <- valid & p > 0 & p < 1
  contrib <- numeric(length(p))
  contrib[seg] <- 1 / harmonic[2 * n_called[seg] - 1]
  dim(contrib) <- dim(pi)
  cbind(contrib = rowSums(contrib), pi = rowSums(pi), n_valid = rowSums(valid))
}

.theta_harmonic <- function(n_samples) cumsum(1 / seq_len(max(2 * n_samples - 1, 1)))

# Per-chunk sums (run by workers when there are any), computed in one pass
# over the chunk's genotype rows by theta_chunk_cpp() (src/theta_counts.cpp),
# without reshaping the chunk into sample x variant matrices. col_idx (the
# chunk's variants) is returned so the main process can check that no variant
# is split across chunk files.
.theta_chunk <- function(fpath, shared) {
  chunk <- arrow::read_feather(fpath, col_select = c(".row_id", "sample", "a1", "a2"),
                               as_data_frame = FALSE)
  s <- .sample_index(chunk, shared$samples)
  vp <- .match_row_id(as.vector(chunk$.row_id), shared$var_pos)
  keep <- !is.na(vp) & !is.na(s)
  if (!any(keep)) return(NULL)
  vp <- vp[keep]
  col_idx <- unique(vp)
  r <- theta_chunk_cpp(match(vp, col_idx), shared$grp[s[keep]],
                       as.vector(chunk$a1)[keep], as.vector(chunk$a2)[keep],
                       length(col_idx), shared$n_pops, shared$harmonic)
  c(list(col_idx = col_idx), r)
}

# Streams the chunk files and adds up their sums. Returns NULL if a variant
# occurs in more than one chunk file, since its per-site values would then be
# computed from partial counts.
.theta_stream <- function(setup) {
  group <- matrix(0, nrow = setup$n_pops, ncol = 3,
                  dimnames = list(NULL, c("contrib", "pi", "n_valid")))
  total <- c(contrib = 0, pi = 0, n_valid = 0)
  nobs <- 0
  seen <- raw(setup$n_var)
  # group (row of P) of each sample in setup$samples, 0 if none
  grp <- integer(setup$n_samples)
  w <- which(setup$P != 0L, arr.ind = TRUE)
  grp[w[, 2]] <- w[, 1]
  shared <- list(samples = setup$samples, var_pos = setup$var_pos, grp = grp,
                 n_pops = setup$n_pops, harmonic = .theta_harmonic(setup$n_samples))
  it <- .chunk_iterator(setup$feather_files, .theta_chunk, shared = shared,
                        label = "Reading chunk")
  on.exit(it$done())
  while (!is.null(w <- it$next_wave())) {
    for (tc in w$res) {
      if (is.null(tc)) next
      if (any(seen[tc$col_idx] != as.raw(0))) return(NULL)
      seen[tc$col_idx] <- as.raw(1)
      group <- group + tc$group
      total <- total + tc$total
      nobs <- nobs + tc$nobs
    }
  }
  list(group = group, total = total, nobs = nobs)
}

# The same sums from complete count matrices (groups x variants), computed in
# blocks of variants so that no full-size double matrix is created.
.theta_from_counts <- function(counts, block = 100000L) {
  n_var <- ncol(counts$alt)
  harmonic <- .theta_harmonic(max(colSums(counts$nobs), 0) / 2)   # all-sample totals are largest
  group <- matrix(0, nrow = nrow(counts$alt), ncol = 3,
                  dimnames = list(NULL, c("contrib", "pi", "n_valid")))
  total <- c(contrib = 0, pi = 0, n_valid = 0)
  for (lo in seq(1L, max(n_var, 1L), by = block)) {
    if (n_var == 0L) break
    cols <- lo:min(lo + block - 1L, n_var)
    alt <- counts$alt[, cols, drop = FALSE]
    nobs <- counts$nobs[, cols, drop = FALSE]
    group <- group + .theta_sums(alt, nobs, harmonic)
    total <- total + .theta_sums(matrix(colSums(alt), nrow = 1),
                                 matrix(colSums(nobs), nrow = 1), harmonic)[1, ]
  }
  list(group = group, total = total, nobs = sum(counts$nobs))
}
