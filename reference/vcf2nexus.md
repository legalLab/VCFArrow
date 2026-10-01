# vcf2nexus

Converts a VCFArrow object to a NEXUS format infile

## Usage

``` r
vcf2nexus(vcf_arrow, out_file, keep_groups = NULL)
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

This function converts a VCFArrow object to an external NEXUS formatted
file. Writing occurs in chunks whose size is determined by the
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
vcf2nexus(vcf, out_file = tempfile(fileext = ".nex"))
#> ℹ Accumulating Nexus: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Nexus file...
#> ✔ Nexus file written to /tmp/RtmpodaNXU/file1ffe47e5c99c.nex
```
