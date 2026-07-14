#' @title vcf2admixture
#'
#' @description
#' Converts a VCFArrow object to a PLINK .bed format infile plus ADMIXTURE-specific .pop file
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param out_file -> name of file to output, default 'plink_out' (character)
#' @param sex -> vector of sexes of samples, default NULL (character)
#' @param pheno -> vector of phenotypes of samples, default NULL (character)
#' @param supervised -> flag to generate .pop for use in ADMIXTURE's supervised mode, default FALSE (Boolean)
#' @param reference_groups -> vector of ancestry group labels, default NULL (character)
#'
#' @return NULL
#'
#' @details
#' This function converts a VCFArrow object to an external PLINK .bed formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Sex and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
#' For supervised analyses, the supervised flag needs to be true to generate
#' .pop file needed for ADMIXTURE's supervised mode.
#' The reference_groups vector specifies ancestral information of individuals.
#' When NULL, or for any missing individual, ancestry is inferred.
#'
#'
#' @examples
#' vcf2admixture(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "admix_out", sex = sex, pheno = pheno, supervised = TRUE, reference_groups = reference_groups)
#' vcf2admixture(vcf_arrow, my_groups, out_file = "admix_out", supervised = TRUE, reference_groups = c("XXX", "YYY"))
#' vcf2admixture(vcf_arrow, my_groups, out_file = "admix_out", supervised = TRUE)
#' vcf2admixture(vcf_arrow, my_groups, out_file = "admix_out")
#' vcf2admixture(vcf_arrow)
#'
#' @export
#'
# ADMIXTURE accepts binary PLINK (.bed/.bim/.fam), ordinary PLINK (.ped/.map),
# or EIGENSTRAT (.geno/.map) as input.  Binary PLINK is the recommended format
# (fastest I/O) and is what this function produces by delegating entirely to
# vcf2plink_bed() for the genotype/variant/sample files.
#
# The only ADMIXTURE-specific artifact not covered by vcf2plink_bed() is the
# .pop file required for supervised ancestry estimation.  This function adds
# that optionally via the supervised / reference_groups parameters.
#
# Files produced
# ──────────────────────────────────────────────────────────────────────────────
#   <out_file>.bed  binary genotype matrix      (via vcf2plink_bed)
#   <out_file>.bim  variant metadata            (via vcf2plink_bed)
#   <out_file>.fam  sample metadata             (via vcf2plink_bed)
#   <out_file>.pop  supervised ancestry labels  (only when supervised = TRUE)
#
# Parameters
# ──────────────────────────────────────────────────────────────────────────────
# supervised      Logical.  If TRUE, write a .pop file alongside the PLINK
#                 binary set.  ADMIXTURE's supervised mode (--supervised flag)
#                 requires this file.  Default FALSE.
#
# reference_groups
#                 Character vector of group labels to treat as known-ancestry
#                 reference populations.  Their group name appears in the .pop
#                 file.  All other groups get "-" (ancestry to be estimated).
#                 If NULL (default) and supervised = TRUE, ALL samples are
#                 treated as reference — i.e. every sample's group label is
#                 written, which is a fully supervised run.
#

vcf2admixture <- function(vcf_arrow, keep_groups = NULL,
                          out_file = "admixture_in",
                          sex = NULL, pheno = NULL,
                          supervised = FALSE,
                          reference_groups = NULL) {

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  # ── Validate supervised / reference_groups interaction ────────────────────
  if (!is.null(reference_groups) && !supervised)
    cli::cli_abort(c(
      "{.arg reference_groups} was supplied but {.arg supervised} is FALSE.",
      "i" = "Set {.code supervised = TRUE} to enable supervised mode."
    ))

  # ── Produce .bed / .bim / .fam via the existing PLINK binary exporter ─────
  # All setup (vcf_export_setup, variant filtering, sample ordering, chunk
  # reads, C++ writers) lives there — no duplication needed here.
  vcf2plink_bed(
    vcf_arrow = vcf_arrow,
    keep_groups = keep_groups,
    out_file = out_file,
    sex = sex,
    pheno = pheno
  )

  # ── .pop file for supervised mode ─────────────────────────────────────────
  if (supervised) {
    # Reconstruct the same sample order vcf2plink_bed() used (.fam row order).
    # .vcf_export_setup() reorders samples group-contiguously; the .pop file
    # must match that exact order row-for-row.
    setup <- .vcf_export_setup(vcf_arrow, keep_groups)

    pop_labels <- if (is.null(reference_groups)) {
      # Fully supervised: every sample's group is a known reference.
      setup$samples_groups
    } else {
      # Semi-supervised: only the nominated groups are reference; all others
      # get "-" (ancestry to be estimated by ADMIXTURE).
      unknown <- !(setup$samples_groups %in% reference_groups)
      labels <- setup$samples_groups
      labels[unknown] <- "-"
      labels
    }

    pop_file <- paste0(out_file, ".pop")
    writeLines(pop_labels, con = pop_file)
    cli::cli_alert_success(
      "ADMIXTURE .pop file written to {.file {pop_file}}"
    )

    n_ref <- sum(pop_labels != "-")
    n_unk <- sum(pop_labels == "-")
    cli::cli_alert_info(
      "Supervised mode: {n_ref} reference sample{?s}, \\
       {n_unk} sample{?s} with ancestry to be estimated."
    )
  }

  invisible(vcf_arrow)
}
