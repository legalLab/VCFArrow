#' @title vcf2genlight
#'
#' @description
#' Converts a VCFArrow object to Genlight format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'genlight.rds' (character)
#' @param ploidy -> ploidy level, default = 2 (integer)
#' @param save ->  save as R data object, default = FALSE (Boolean)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external Genlight formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Genlight objects can encode polyploid genomes, by default diploid genomes are assumed.
#' Genlight objects are in memory S4 objects, thus are returned as such, but
#' but optionally may be saved as R data objects.
#'
#' @examples
#' vcf2genlight(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "genlight.rds")
#' vcf2genlight(vcf_arrow, my_groups)
#' vcf2genlight(vcf_arrow)
#'
#' @export
#'

vcf2genlight <- function(vcf_arrow, keep_groups = NULL,
                         out_file = "genlight.rds",
                         ploidy = 2L, save = FALSE) {

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Genlight")

  cli::cli_alert_info("Building genlight object...")

  # 0+0=0 hom-ref, 0+1 or 1+0=1 het, 1+1=2 hom-alt; NA propagates naturally
  geno_mat <- acc$a1 + acc$a2
  storage.mode(geno_mat) <- "integer"
  rownames(geno_mat) <- setup$samples

  x <- suppressWarnings(
    new("genlight",
        gen = lapply(seq_len(nrow(geno_mat)), function(i) geno_mat[i, ]))
  )
  adegenet::indNames(x) <- setup$samples
  adegenet::chromosome(x) <- setup$variants$CHROM
  adegenet::position(x) <- setup$variants$POS
  adegenet::locNames(x) <- setup$loci
  adegenet::pop(x) <- factor(setup$samples_groups, levels = setup$group_names)
  adegenet::ploidy(x) <- ploidy
  adegenet::strata(x) <- data.frame(pop = adegenet::pop(x))

  if (save) {
    saveRDS(x, file = out_file)
    cli::cli_alert_success("Genlight object saved to {.file {out_file}}")
  }

  return(x)
}
