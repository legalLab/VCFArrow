# vcf2smartsnp

Converts a VCFArrow object to smartsnp tabular format infile

## Usage

``` r
vcf2smartsnp(vcf_arrow, out_file, keep_groups = NULL)
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

This function converts a VCFArrow object to an external SmartSNP
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
vcf2smartsnp(vcf, out_file = tempfile(fileext = ".txt"))
#> ℹ Building SmartSNP: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing SmartSNP file...
#> ✔ SmartSNP file written to /tmp/RtmpodaNXU/file1ffe6c79bb72.txt
```
