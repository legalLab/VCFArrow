#' @title vcf2genlight
#'
#' @description
#' Converts a VCFArrow object to Genlight format infile
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> name of .rds file to write, required when save = TRUE,
#'   default NULL (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param ploidy -> ploidy level, default = 2 (integer)
#' @param save ->  save as R data object, default = FALSE (Boolean)
#'
#' @importClassesFrom adegenet genlight
#'
#' @return An adegenet `genlight` object.
#'
#' @details
#' This function converts a VCFArrow object to an adegenet Genlight object.
#' Genotypes are read in chunks whose size is determined by the read_vcf() function.
#' If no groups are defined, the default behavior is to use all groups.
#' Genlight objects can encode polyploid genomes, by default diploid genomes are assumed.
#' Genlight objects are in-memory S4 objects and are always returned; if
#' save = TRUE, the object is also written to out_file with saveRDS().
#' Unlike the other exporters, out_file is optional because nothing is written
#' unless save = TRUE.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' gl <- vcf2genlight(vcf)
#' gl
#' gl <- vcf2genlight(vcf, out_file = tempfile(fileext = ".rds"),
#'                    save = TRUE)
#'
#' @export
#'

vcf2genlight <- function(vcf_arrow, out_file = NULL, keep_groups = NULL,
                         ploidy = 2L, save = FALSE) {

  if (save && is.null(out_file)) {
    cli::cli_abort("{.arg out_file} must be supplied when {.code save = TRUE}.")
  }

  setup <- .vcf_export_setup(vcf_arrow, keep_groups)
  acc <- .accumulate_individuals(setup, "Genlight")

  cli::cli_alert_info("Building Genlight object...")

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
    cli::cli_alert_info("Writing Genlight object...")
    saveRDS(x, file = out_file)
    cli::cli_alert_success("Genlight object written to {.file {out_file}}")
  }

  return(x)
}
