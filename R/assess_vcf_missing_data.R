#' @title assess_vcf_missing_data
#'
#' @description
#' Quantifying missing data of all samples in VCF.
#' Inspired by https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param res_path -> path to results (directory for output dataframe and plots)
#' @param species -> sample name for plot (character)
#' @param project -> project name / base output file name (character)
#' @param details -> flag for adding project name into figure title (Boolean)
#'
#' @return dataframe and plot of missing data for each sample in VCF
#'
#' @details
#' This function generates a dataframe of absolute and relative missing data
#' per sample, and a plot of relative missing data % per sample.
#'
#' @examples
#' assess_vcf_missing_data(vcf_arrow = my_vcf, res_path = my_res_path, species = species_name, project = project_name)
#' assess_vcf_missing_data(my_vcf, my_res_path, species_name, project_name, details = TRUE)
#' assess_vcf_missing_data(my_vcf, my_res_path, species_name, project_name)
#'

assess_vcf_missing_data <- function(vcf_arrow, res_path, species, project, details = FALSE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  gt <- vcf_arrow@gt

  # ensure we have group info
  if (any(is.na(vcf_arrow@groups))) {
    cli::cli_alert_warning("Missing group assignments. Use set_vcf_groups() to fill the 'groups' slot.")
  }

  group_df <- tibble::tibble(
    sample = vcf_arrow@samples,
    group = vcf_arrow@groups
  )

  # compute missingness per sample
  miss_df <- gt |>
    dplyr::group_by(sample) |>
    dplyr::summarise(
      missing_n = sum(is.na(a1)),
      total_loci = dplyr::n(),
      missing_p = missing_n / total_loci,
      .groups = "drop"
    ) |>
    dplyr::collect() |>
    dplyr::left_join(group_df, by = "sample") |>
    dplyr::relocate(sample, group, missing_n, missing_p, .before = everything()) |>
    dplyr::arrange(group, sample)

  # save table
  write.table(
    miss_df,
    file = paste0(res_path, project, "_missingness.csv"),
    row.names = FALSE,
    quote = FALSE,
    sep = ","
  )

  # ordering for plot
  ordr <- miss_df |>
    dplyr::arrange(group) |>
    dplyr::pull(sample)

  # include project name into ggplot2 figure title
  plot_title <- if (details == TRUE) {
    bquote(italic(.(species)) ~ "(" * .(project) * ")")
  } else {
    bquote(italic(.(species)))
  }

  # make the plot
  plt <- ggplot2::ggplot(
    miss_df,
    ggplot2::aes(x = factor(sample, levels = ordr), y = missing_p, fill = group)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1, size = ggplot2::rel(0.8)
      )
    ) +
    ggplot2::labs(
      x = "sample",
      y = "missing",
      title = plot_title
    )

  # save plots
  ggplot2::ggsave(plt, filename = paste0(res_path, project, "_missingness.pdf"),
                  device = "pdf", width = 6, height = 4, bg = "transparent", limitsize = FALSE)

  ggplot2::ggsave(plt, filename = paste0(res_path, project, "_missingness.svg"),
                  device = "svg", width = 6, height = 4, bg = "transparent", limitsize = FALSE)

  ggplot2::ggsave(plt, filename = paste0(res_path, project, "_missingness.png"),
                  device = "png", width = 6, height = 4, bg = "transparent", limitsize = FALSE)

  # FUTURE - make violin plot of read depths

  invisible(vcf_arrow)
}
