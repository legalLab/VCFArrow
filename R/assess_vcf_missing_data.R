#' @title assess_vcf_missing_data
#'
#' @description
#' Quantifying missing data of all samples in VCF
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param res_path -> path to results (directory for output dataframe and plots)
#' @param species -> species name, the plot title (character)
#' @param project -> project name / base output file name (character)
#' @param details -> flag for adding project name into figure title, default TRUE (Boolean)
#'
#' @return Invisibly returns the input VCFArrow object; called for its side
#'   effects of writing a table of the missing data of each sample to
#'   `<project>_missingness.csv` and a bar plot of it to
#'   `<project>_missingness.pdf`, `.svg` and `.png` in `res_path`.
#'
#' @details
#' This function generates a table of absolute and relative missing data
#' per sample (columns sample, group, missing_n, missing_p and total_loci,
#' missing_p being the proportion of missing genotypes), and a bar plot of
#' the proportion of missing data per sample.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' assess_vcf_missing_data(vcf, res_path = tempdir(),
#'                         species = "Phyllomedusa vaillantii",
#'                         project = "vaillantii")
#'
#' @export
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
  .save_plot(plt, file.path(res_path, paste0(project, "_missingness")))

  invisible(vcf_arrow)
}
