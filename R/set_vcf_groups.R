#' @title set_vcf_groups
#'
#' @description
#' Associate samples with group info based on "strt" file.
#' The "strt" is a tsv file with 2+ columns, with 2 columns name id and pop.
#' By default strt is named "strata" and found in the datapath where the
#' VCF file is found.
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param data_path -> path to data (where strata are located)
#' @param strt -> file with 2+ columns, one id column and one pop column (tsv file, header = TRUE)
#'
#' @return vcf_arrow object with groups slot filled
#'
#' @details
#' This function creates list of factors of sample-to-group assignments and
#' groups based on a strata and group file.
#'
#' @examples
#' get_vcf_group_info(vcf = my_vcf, data_path = my_data_path, strt = "strata")
#' get_vcf_group_info(my_vcf, my_data_path, strata)
#' get_vcf_group_info(my_vcf, my_data_path)
#'

set_vcf_groups <- function(vcf_arrow, data_path, strt = "strata") {

  if (!inherits(vcf_arrow, "VCFArrow")) {
    cli::cli_abort("Expecting a VCFArrow object")
  }

  # read sample to group assignment
  strata_df <- read.table(paste0(data_path, strt), header = TRUE) |>
    tibble::as_tibble() |>
    dplyr::mutate(id = as.character(id))

  # build sample table from VCFArrow
  sample_df <- tibble::tibble(id = vcf_arrow@samples)

  # join group info
  strata <- sample_df |>
    dplyr::left_join(strata_df, by = "id")

  # check if all samples in vcf are assigned to groups
  if (any(is.na(strata$pop) == TRUE)) {
    cli::cli_alert_warning("VCF has individuals not assigned to a group: {strata[is.na(strata$pop) == TRUE,]}")
  } else {
    vcf_arrow@groups <- as.character(strata$pop)
  }

  return(vcf_arrow)
}
