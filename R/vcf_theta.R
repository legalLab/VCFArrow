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

vcf_theta <- function(gt, samples, grps) {

  # map samples → groups
  sample_map <- data.frame(
    sample = samples,
    group = as.character(grps),
    stringsAsFactors = FALSE
  )

  map_tbl <- arrow::arrow_table(sample_map)

  # attach group info
  ds_g <- gt |>
    dplyr::inner_join(map_tbl, by = "sample")

  # per (variant × group) summaries
  by_var_grp <- ds_g |>
    dplyr::mutate(
      called = !(is.na(a1) | is.na(a2)),
      allele_sum = coalesce(a1, 0) + coalesce(a2, 0)
    ) |>
    dplyr::group_by(.row_id, group) |>
    dplyr::summarise(
      n_called = sum(called),
      alt_count = sum(allele_sum),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      p = dplyr::if_else(n_called > 0, alt_count / (2 * n_called), NA_real_),
      valid = n_called >= 2,
      pi = dplyr::if_else(valid, 2 * p * (1 - p), NA_real_),
      seg = (p > 0 & p < 1) & valid
    )

  by_var_grp <- by_var_grp |>
    dplyr::collect()

  # harmonic numbers
  max_n <- max(by_var_grp$n_called, na.rm = TRUE)
  harmonic <- c(NA, cumsum(1 / (1:(2 * max_n - 1))))

  # theta per variant
  by_var_grp$a1 <- NA_real_
  idx <- which(by_var_grp$valid)
  by_var_grp$a1[idx] <- harmonic[2 * by_var_grp$n_called[idx] - 1]

  by_var_grp$contrib <- 1 / by_var_grp$a1
  by_var_grp$contrib[!by_var_grp$seg] <- 0

  # group-level summaries
  theta_group <- by_var_grp |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      theta_w = sum(contrib, na.rm = TRUE) / sum(valid, na.rm = TRUE),
      pi = mean(pi, na.rm = TRUE),
      .groups = "drop"
    )

  # global (across all samples)
  global <- gt |>
    dplyr::mutate(
      called = !(is.na(a1) | is.na(a2)),
      allele_sum = coalesce(a1, 0) + coalesce(a2, 0)
    ) |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      n_called = sum(called),
      alt_count = sum(allele_sum),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      p = dplyr::if_else(n_called > 0, alt_count / (2 * n_called), NA_real_),
      valid = n_called >= 2,
      pi = dplyr::if_else(valid, 2 * p * (1 - p), NA_real_),
      seg = (p > 0 & p < 1) & valid
    ) |>
    dplyr::collect()

  max_n_g <- max(global$n_called, na.rm = TRUE)
  harmonic_g <- c(NA, cumsum(1 / (1:(2 * max_n_g - 1))))

  global$a1 <- NA_real_
  idx <- which(global$valid)
  global$a1[idx] <- harmonic_g[2 * global$n_called[idx] - 1]

  contrib <- 1 / global$a1
  contrib[!global$seg] <- 0

  theta_w <- sum(contrib, na.rm = TRUE) / sum(global$valid, na.rm = TRUE)
  theta_pi <- mean(global$pi, na.rm = TRUE)

  # expand per-sample output
  groups_char <- as.character(grps)

  theta_g <- data.frame(
    group = groups_char,
    n_ind = as.integer(table(groups_char)[groups_char]),
    theta_w = theta_group$theta_w[match(groups_char, theta_group$group)],
    pi = theta_group$pi[match(groups_char, theta_group$group)]
  )

  return(list(
    pi = theta_pi,
    theta_w = theta_w,
    theta_g = theta_g
  ))
}
