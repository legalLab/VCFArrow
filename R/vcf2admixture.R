#' @title vcf2admixture
#'
#' @description
#' Converts a VCFArrow object to a PLINK .bed format infile plus ADMIXTURE-specific .pop file
#'
#' @author Tomas Hrbek May 2026
#'
#' @param vcf_arrow -> VCFArrow object
#' @param out_file -> output file prefix; extensions are added automatically, 
#'   no default (character)
#' @param keep_groups -> groups to retain, default NULL (character)
#' @param sex -> vector of sexes of samples, default NULL (character)
#' @param pheno -> vector of phenotypes of samples, default NULL (character)
#' @param chrom_code -> how CHROM is written: "auto" (keep numeric CHROM,
#'   otherwise 0), "index" (1..n by first appearance), "zero" (all 0) or
#'   "keep" (verbatim, only if every CHROM is a number), default "zero"
#'   (character)
#' @param supervised -> flag to generate .pop for use in ADMIXTURE's 
#'   supervised mode, default FALSE (Boolean)
#' @param reference_groups -> groups whose samples are of known ancestry in a
#'   supervised run; their names must begin with a letter or digit, default
#'   NULL (all groups, a fully supervised run, with a warning) (character)
#'
#' @return Invisibly returns the input VCFArrow object; called for 
#'   its side effect of writing `out_file`.
#'
#' @details
#' This function converts a VCFArrow object to an external PLINK .bed formatted file.
#' Writing occurs in chunks whose size is determined by the read_vcf() function.
#' Larger chunks result in faster writing speeds.
#' If no groups are defined, the default behavior is to use all groups.
#' Only biallelic SNVs are exported (multiallelic sites and indels are left
#' out), and every sample must be assigned to a group (see set_vcf_groups()),
#' also when keep_groups is NULL.
#' Sex and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
#' For supervised analyses, the supervised flag needs to be true to generate
#' .pop file needed for ADMIXTURE's supervised mode.
#' The .pop file has one line per sample, in the order of the .fam file: the
#' sample's group, its known ancestry, if the group is in reference_groups,
#' and '-' (ancestry to be estimated) otherwise. With reference_groups = NULL
#' (the default) every sample's group is written, which is a fully supervised
#' run: the ancestry of every sample is fixed and only the allele frequencies
#' are estimated, so a warning is given. reference_groups can only be given
#' with supervised = TRUE, must be groups of the exported samples, and must
#' begin with a letter or digit, as ADMIXTURE reads a .pop line beginning
#' with any other character as unknown ancestry. In a supervised run, K
#' (given to ADMIXTURE) must equal the number of reference groups.
#' ADMIXTURE (from version 1.3) requires numeric chromosome names, so
#' chrom_code = "keep" is only allowed if every CHROM is a number; with
#' "zero" (the default), "auto" or "index" the names are always numeric.
#' Haploid analyses (ADMIXTURE's --haploid option) of sex chromosomes need
#' their numeric codes (e.g. 23 and 24 in humans), which "auto" and "keep"
#' retain.
#' ADMIXTURE does not model linkage disequilibrium, and its manual
#' recommends thinning linked markers, e.g. keeping one SNV per locus with
#' vcf_filter_oneSNV() before export.
#'
#' @examples
#' f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
#'                  package = "VCFArrow")
#' vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#' vcf2admixture(vcf, out_file = file.path(tempdir(), "admixture_in"))
#' vcf2admixture(vcf, out_file = file.path(tempdir(), "admixture_sup"),
#'               supervised = TRUE,
#'               reference_groups = c("GS", "BS", "WA"))
#'
#' @export
#'
# ADMIXTURE accepts binary PLINK (.bed/.bim/.fam), ordinary PLINK (.ped/.map),
# or EIGENSTRAT (.geno/.map) as input. Binary PLINK is the recommended format
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
# supervised      Logical. If TRUE, write a .pop file alongside the PLINK
#                 binary set. ADMIXTURE's supervised mode (--supervised flag)
#                 requires this file. Default FALSE.
#
# reference_groups
#                 Character vector of group labels to treat as known-ancestry
#                 reference populations. Their group name appears in the .pop
#                 file. All other groups get "-" (ancestry to be estimated).
#                 If NULL (default) and supervised = TRUE, ALL samples are
#                 treated as reference — i.e. every sample's group label is
#                 written, which is a fully supervised run.
#

vcf2admixture <- function(vcf_arrow, out_file, keep_groups = NULL,
                          sex = NULL, pheno = NULL,
                          chrom_code = c("zero", "auto", "index", "keep"),
                          supervised = FALSE,
                          reference_groups = NULL) {

  if (missing(out_file)) cli::cli_abort("{.arg out_file} must be supplied.")

  chrom_code <- match.arg(chrom_code)

  if (!inherits(vcf_arrow, "VCFArrow"))
    cli::cli_abort("Expecting a VCFArrow object")

  # ── Validate supervised / reference_groups interaction ────────────────────
  if (!is.null(reference_groups) && !supervised)
    cli::cli_abort(c(
      "{.arg reference_groups} was supplied but {.arg supervised} is FALSE.",
      "i" = "Set {.code supervised = TRUE} to enable supervised mode."
    ))

  # Same samples, sample order (.fam row order) and variants as
  # vcf2plink_bed(); checked before any file is written.
  setup <- .vcf_export_setup(vcf_arrow, keep_groups)

  # ADMIXTURE (>= 1.3) needs numeric chromosome names
  if (chrom_code == "keep") {
    chroms <- unique(as.character(setup$variants$CHROM))
    bad <- chroms[!grepl("^[0-9]+$", chroms)]
    if (length(bad) > 0L)
      cli::cli_abort(c(
        "ADMIXTURE requires numeric chromosome names, but {length(bad)} \\
         CHROM value{?s} {?is/are} not {?a number/numbers}, e.g. \\
         {.val {utils::head(bad, 3)}}.",
        "i" = "Use {.code chrom_code = \"zero\"}, {.val auto} or {.val index}."
      ))
  }

  if (supervised) {
    groups <- unique(setup$samples_groups)
    if (is.null(reference_groups)) {
      cli::cli_warn(c(
        "{.arg reference_groups} is NULL: every sample is a reference, a \\
         fully supervised run in which only the allele frequencies are \\
         estimated.",
        "i" = "Give {.arg reference_groups} to estimate the ancestry of the \\
               samples of the other groups."
      ))
      ref <- groups
    } else {
      unknown_ref <- setdiff(reference_groups, groups)
      if (length(unknown_ref) > 0L)
        cli::cli_abort(
          "{.arg reference_groups} {?is/are} not {?a group/groups} of the \\
           exported samples: {.val {unknown_ref}}."
        )
      ref <- unique(reference_groups)
    }
    bad_names <- ref[!grepl("^[[:alnum:]]", ref)]
    if (length(bad_names) > 0L)
      cli::cli_abort(c(
        "Reference group name{?s} {.val {bad_names}} {?does/do} not begin \\
         with a letter or digit; ADMIXTURE would read {?it/them} as unknown \\
         ancestry.",
        "i" = "{cli::qty(length(bad_names))}Rename the group{?s} with \\
               {.fn set_vcf_groups}."
      ))
  }

  # ── Produce .bed / .bim / .fam via the existing PLINK binary exporter ─────
  vcf2plink_bed(
    vcf_arrow = vcf_arrow,
    keep_groups = keep_groups,
    out_file = out_file,
    sex = sex,
    pheno = pheno,
    chrom_code = chrom_code
  )

  # ── .pop file for supervised mode ─────────────────────────────────────────
  if (supervised) {
    # one line per sample in .fam order: the group of reference samples,
    # "-" (ancestry to be estimated) for the others
    pop_labels <- setup$samples_groups
    pop_labels[!(pop_labels %in% ref)] <- "-"

    pop_file <- paste0(out_file, ".pop")
    writeLines(pop_labels, con = pop_file)
    cli::cli_alert_success(
      "ADMIXTURE .pop file written to {.file {pop_file}}"
    )

    n_ref <- sum(pop_labels != "-")
    n_unk <- sum(pop_labels == "-")
    cli::cli_alert_info(
      "Supervised mode: {n_ref} reference sample{?s} in {length(ref)} \\
       group{?s} (use K = {length(ref)}), {n_unk} sample{?s} with ancestry \\
       to be estimated."
    )
  }

  invisible(vcf_arrow)
}
