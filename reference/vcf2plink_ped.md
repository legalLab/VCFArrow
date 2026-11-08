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
  chrom_code = c("auto", "index", "zero", "keep"),
  allele_code = c("nucleotide", "12")
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

- allele_code:

  -\> how alleles are written in the .ped file: "nucleotide" (the REF
  and ALT bases) or "12" (1 = ALT, 2 = REF, as ADMIXTURE requires for
  .ped input), default "nucleotide" (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external PLINK .ped
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Sex
and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
Missing alleles are written as 0. In the "12" coding allele 1 is ALT and
allele 2 is REF, as in vcf2plink_bed(). With either coding, PLINK 1.9
makes the minor allele A1 when it loads a .ped file; PLINK 2 (and PLINK
1.9 with –keep-allele-order) take the "12" coding's 1 as A1. PLINK reads
both allele codings; ADMIXTURE reads .ped files only in the "12" coding
(allele_code = "12"), and needs integer chromosome codes (chrom_code),
as for .bed input.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
vcf2plink_ped(vcf, out_file = file.path(tempdir(), "plink_out"))
#> ℹ Accumulating PLINK .ped: 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing PLINK file...
#> ℹ Chromosome names were recoded as integers; mapping written to /tmp/RtmpVdrUAs/plink_out.chrommap
#> ✔ PLINK text fileset written to /tmp/RtmpVdrUAs/plink_out.ped, /tmp/RtmpVdrUAs/plink_out.map
```
