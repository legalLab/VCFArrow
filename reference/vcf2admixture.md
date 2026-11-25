# vcf2admixture

Converts a VCFArrow object to a PLINK .bed format infile plus
ADMIXTURE-specific .pop file

## Usage

``` r
vcf2admixture(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  sex = NULL,
  pheno = NULL,
  chrom_code = c("zero", "auto", "index", "keep"),
  supervised = FALSE,
  reference_groups = NULL
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> output file prefix; extensions are added automatically, no default
  (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- sex:

  -\> vector of sexes of samples, default NULL (character)

- pheno:

  -\> vector of phenotypes of samples, default NULL (character)

- chrom_code:

  -\> how CHROM is written: "auto" (keep numeric CHROM, otherwise 0),
  "index" (1..n by first appearance), "zero" (all 0) or "keep"
  (verbatim, only if every CHROM is a number), default "zero"
  (character)

- supervised:

  -\> flag to generate .pop for use in ADMIXTURE's supervised mode,
  default FALSE (Boolean)

- reference_groups:

  -\> groups whose samples are of known ancestry in a supervised run;
  their names must begin with a letter or digit, default NULL (all
  groups, a fully supervised run, with a warning) (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external PLINK .bed
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Only
biallelic SNVs are exported (multiallelic sites and indels are left
out), and every sample must be assigned to a group (see
set_vcf_groups()), also when keep_groups is NULL. Sex and phenotype
vectors are optional. If not defined sex = 0, pheno = -9. For supervised
analyses, the supervised flag needs to be true to generate .pop file
needed for ADMIXTURE's supervised mode. The .pop file has one line per
sample, in the order of the .fam file: the sample's group, its known
ancestry, if the group is in reference_groups, and '-' (ancestry to be
estimated) otherwise. With reference_groups = NULL (the default) every
sample's group is written, which is a fully supervised run: the ancestry
of every sample is fixed and only the allele frequencies are estimated,
so a warning is given. reference_groups can only be given with
supervised = TRUE, must be groups of the exported samples, and must
begin with a letter or digit, as ADMIXTURE reads a .pop line beginning
with any other character as unknown ancestry. In a supervised run, K
(given to ADMIXTURE) must equal the number of reference groups.
ADMIXTURE (from version 1.3) requires numeric chromosome names, so
chrom_code = "keep" is only allowed if every CHROM is a number; with
"zero" (the default), "auto" or "index" the names are always numeric.
Haploid analyses (ADMIXTURE's –haploid option) of sex chromosomes need
their numeric codes (e.g. 23 and 24 in humans), which "auto" and "keep"
retain. ADMIXTURE does not model linkage disequilibrium, and its manual
recommends thinning linked markers, e.g. keeping one SNV per locus with
vcf_filter_oneSNV() before export.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2admixture(vcf, out_file = file.path(tempdir(), "admixture_in"))
#> ℹ Building PLINK: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing PLINK files...
#> ℹ Chromosome names were recoded as integers; mapping written to /tmp/Rtmpq4egqB/admixture_in.chrommap
#> ✔ PLINK binary fileset written to /tmp/Rtmpq4egqB/admixture_in.bed, /tmp/Rtmpq4egqB/admixture_in.bim, /tmp/Rtmpq4egqB/admixture_in.fam
vcf2admixture(vcf, out_file = file.path(tempdir(), "admixture_sup"),
              supervised = TRUE,
              reference_groups = c("GS", "BS", "WA"))
#> ℹ Building PLINK: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing PLINK files...
#> ℹ Chromosome names were recoded as integers; mapping written to /tmp/Rtmpq4egqB/admixture_sup.chrommap
#> ✔ PLINK binary fileset written to /tmp/Rtmpq4egqB/admixture_sup.bed, /tmp/Rtmpq4egqB/admixture_sup.bim, /tmp/Rtmpq4egqB/admixture_sup.fam
#> ✔ ADMIXTURE .pop file written to /tmp/Rtmpq4egqB/admixture_sup.pop
#> ℹ Supervised mode: 15 reference samples in 3 groups (use K = 3), 3 samples with ancestry to be estimated.
```
