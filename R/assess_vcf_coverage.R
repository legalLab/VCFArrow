#' @title assess_vcf_coverage
#'
#' @description
#' Quantifying read depth of all samples in VCF.
#' Inspired by https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param res_path -> path to results (directory for output dataframe and plots)
#' @param species -> sample name for plot (character)
#' @param project -> project name / base output file name (character)
#' @param details -> flag for adding project name into figure title, default TRUE (Boolean)
#' @param max_points_per_sample -> maximum number of SNVs from which to plot coverage, default 5000 (integer)
#'
#' @return violin plot of read depth for each sample in VCF
#'
#' @details
#' This function generates a violin plot of read depths per sample.
#'
#' @examples
#' assess_vcf_coverage(vcf_arrow = my_vcf, res_path = my_res_path, species = species_name, project = project_name)
#' assess_vcf_coverage(my_vcf, my_res_path, species_name, project_name, details = TRUE, max_points_per_sample = 5000L)
#' assess_vcf_coverage(my_vcf, my_res_path, species_name, project_name)
#'

assess_vcf_coverage <- function(vcf_arrow, res_path, species, project,
                                details = TRUE,
                                max_points_per_sample = 5000L) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  group_df <- tibble::tibble(sample = vcf_arrow@samples, group = vcf_arrow@groups)

  # materialise a bounded (sample, DP) sample directly from feather chunks
  scan <- .scan_vcf_gt(vcf_arrow, collect_dp = TRUE,
                       max_points_per_sample = max_points_per_sample)

  dp_df <- scan$dp_df

  if (is.null(dp_df) || nrow(dp_df) == 0L)
    cli::cli_abort("No non-missing DP values found — cannot build coverage plot.")

  dp_df <- dp_df |>
    dplyr::left_join(group_df, by = "sample")

  # Ordering consistent with assess_vcf_missing_data()'s bar plot
  ordr <- group_df |>
    dplyr::arrange(group) |>
    dplyr::pull(sample)

  plot_title <- if (isTRUE(details)) {
    bquote(italic(.(species)) ~ "(" * .(project) * ")")
  } else {
    bquote(italic(.(species)))
  }

  # make coverage violin plot
  plt <- ggplot2::ggplot(
    dp_df,
    ggplot2::aes(x = factor(sample, levels = ordr), y = log10(DP), fill = group, color = group)
  ) +
    ggplot2::geom_violin(trim = FALSE, linewidth = 0.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1, size = ggplot2::rel(0.8)
      )
    ) +
    ggplot2::labs(
      x = "sample",
      y = "log10 coverage",
      title = plot_title
    )

  for (ext in c("pdf", "svg", "png")) {
    suppressWarnings(
      ggplot2::ggsave(
        plt,
        filename = file.path(res_path, paste0(project, "_coverage.", ext)),
        device = ext,
        width = 6,
        height = 4,
        bg = "transparent",
        limitsize = FALSE
      )
    )
  }

  invisible(vcf_arrow)
}
