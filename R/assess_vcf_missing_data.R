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
#' @param details -> flag for adding project name into figure title, default TRUE (Boolean)
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

assess_vcf_missing_data <- function(vcf_arrow, res_path, species, project,
                                    details = TRUE) {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # ensure we have group info
  if (any(is.na(vcf_arrow@groups))) {
    cli::cli_alert_warning("Missing group assignments. Use set_vcf_groups() to fill the 'groups' slot.")
  }

  group_df <- tibble::tibble(
    sample = vcf_arrow@samples,
    group = vcf_arrow@groups
  )

  # compute missingness per sample directly from feather chunks
  scan <- .scan_vcf_gt(vcf_arrow, collect_dp = FALSE)

  miss_df <- tibble::tibble(
    sample = vcf_arrow@samples,
    missing_n = scan$missing_n[vcf_arrow@samples],
    total_loci = scan$total_loci[vcf_arrow@samples]
    ) |>
    dplyr::mutate(missing_p = missing_n / total_loci) |>
    dplyr::left_join(group_df, by = "sample") |>
    dplyr::relocate(sample, group, missing_n, missing_p, .before = dplyr::everything()) |>
    dplyr::arrange(group, sample)

  # save table
  write.table(
    miss_df,
    file = file.path(res_path, paste0(project, "_missingness.csv")),
    row.names = FALSE,
    quote = FALSE,
    sep = ","
  )

  # ordering for plot
  ordr <- miss_df |>
    dplyr::arrange(group) |>
    dplyr::pull(sample)

  # include project name into ggplot2 figure title
  plot_title <- if (isTRUE(details)) {
    bquote(italic(.(species)) ~ "(" * .(project) * ")")
  } else {
    bquote(italic(.(species)))
  }

  # make missingness bar plot
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
  for (ext in c("pdf", "svg", "png")) {
    ggplot2::ggsave(
      plt,
      filename = file.path(res_path, paste0(project, "_missingness.", ext)),
      device = ext,
      width = 6,
      height = 4,
      bg = "transparent",
      limitsize = FALSE
    )
  }

  invisible(vcf_arrow)
}
