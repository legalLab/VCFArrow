# vcf2related

Converts a VCFArrow object to Related format infile

## Usage

``` r
vcf2related(vcf_arrow, out_file, keep_groups = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external Related
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Only
biallelic SNVs are exported (multiallelic sites and indels are left
out), and every sample must be assigned to a group (see
set_vcf_groups()), also when keep_groups is NULL. Each line holds a
sample name followed by two columns per locus, the alleles coded A = 1,
C = 2, G = 3, T = 4 (written 01-04); a missing genotype is '0 0', as
related's readgenotypedata() expects.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2related(vcf, out_file = tempfile(fileext = ".txt"))
#> ℹ Accumulating Related: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Related file...
#> ✔ Related file written to /tmp/Rtmpq4egqB/file1db21a02541c.txt
```
