# vcf2plink_ped

Converts a VCFArrow object to a PLINK .ped format infile

## Usage

``` r
vcf2plink_ped(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  sex = NULL,
  pheno = NULL,
  chrom_code = c("auto", "index", "zero", "keep")
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
  (verbatim), default "auto" (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external PLINK .ped
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Sex
and phenotype vectors are optional. If not defined sex = 0, pheno = -9.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2plink_ped(vcf, out_file = file.path(tempdir(), "plink_out"))
#> ℹ Accumulating PLINK .ped: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing PLINK file...
#> ℹ Chromosome names were recoded as integers for PLINK compatibility; mapping written to /tmp/RtmpFwXb6F/plink_out.chrommap
#> ✔ PLINK text fileset written to /tmp/RtmpFwXb6F/plink_out.ped, /tmp/RtmpFwXb6F/plink_out.map
```
