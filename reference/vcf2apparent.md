# vcf2apparent

Converts a VCFArrow object to an Apparent format infile

## Usage

``` r
vcf2apparent(vcf_arrow, out_file, keep_groups = NULL, key = "All")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- key:

  -\> relationship type (All, Pa, Mo, Fa, Off), default All (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external SmartSNP
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.
Possible relationships defined by the parameter 'kee' are All, Pa, Mo,
Fa, Off.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2apparent(vcf, out_file = tempfile(fileext = ".txt"))
#> ℹ Accumulating Apparent: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Apparent file...
#> ✔ Apparent file written to /tmp/RtmpSpcr9J/file1ff87fca44c2.txt
```
