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
  (verbatim), default "zero" (character)

- supervised:

  -\> flag to generate .pop for use in ADMIXTURE's supervised mode,
  default FALSE (Boolean)

- reference_groups:

  -\> vector of ancestry group labels, default NULL (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external PLINK .bed
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Sex
and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
For supervised analyses, the supervised flag needs to be true to
generate .pop file needed for ADMIXTURE's supervised mode. The
reference_groups vector specifies ancestral information of individuals.
When NULL, or for any missing individual, ancestry is inferred.

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
#> ℹ Chromosome names were recoded as integers for PLINK compatibility; mapping written to /tmp/RtmpFwXb6F/admixture_in.chrommap
#> ✔ PLINK binary fileset written to /tmp/RtmpFwXb6F/admixture_in.bed, /tmp/RtmpFwXb6F/admixture_in.bim, /tmp/RtmpFwXb6F/admixture_in.fam
vcf2admixture(vcf, out_file = file.path(tempdir(), "admixture_sup"),
              supervised = TRUE,
              reference_groups = c("GS", "BS", "WA"))
#> ℹ Building PLINK: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing PLINK files...
#> ℹ Chromosome names were recoded as integers for PLINK compatibility; mapping written to /tmp/RtmpFwXb6F/admixture_sup.chrommap
#> ✔ PLINK binary fileset written to /tmp/RtmpFwXb6F/admixture_sup.bed, /tmp/RtmpFwXb6F/admixture_sup.bim, /tmp/RtmpFwXb6F/admixture_sup.fam
#> ✔ ADMIXTURE .pop file written to /tmp/RtmpFwXb6F/admixture_sup.pop
#> ℹ Supervised mode: 15 reference samples, 3 samples with ancestry to be estimated.
```
