# vcf2eigenstrat

Converts a VCFArrow object to Eigenstrat format infile

## Usage

``` r
vcf2eigenstrat(vcf_arrow, out_file, keep_groups = NULL, sex = NULL)
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

  -\> sex of the individual, default = U (undefined) (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external BayesAss
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2eigenstrat(vcf, out_file = file.path(tempdir(), "eigenstrat_in"))
#> ℹ Building EIGENSTRAT: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing EIGENSTRAT files...
#> ✔ EIGENSTRAT fileset written to /tmp/RtmpFwXb6F/eigenstrat_in.geno, /tmp/RtmpFwXb6F/eigenstrat_in.ind, /tmp/RtmpFwXb6F/eigenstrat_in.snp
```
