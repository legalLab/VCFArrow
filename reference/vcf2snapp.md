# vcf2snapp

Converts a VCFArrow object to a SNAPP encoded NEXUS format infile

## Usage

``` r
vcf2snapp(vcf_arrow, out_file, keep_groups = NULL)
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

This function converts a VCFArrow object to an external SNAPP/SNAPPER
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.
Genotypes are coded as SNAPP expects for diploid data: 0 (REF
homozygote), 1 (heterozygote), 2 (ALT homozygote), ? (missing), in a
NEXUS file declared as integer data with symbols 0, 1 and 2 (FORMAT
DATATYPE=INTEGERDATA SYMBOLS="012"), which BEAUti imports with 3 states.

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
vcf2snapp(vcf, out_file = tempfile(fileext = ".nex"))
#> ℹ Accumulating SNAPP: 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing SNAPP file...
#> ✔ SNAPP file written to /tmp/RtmpVdrUAs/file3153de55924.nex
```
