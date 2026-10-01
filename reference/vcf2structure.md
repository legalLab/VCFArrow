# vcf2structure

Converts a VCFArrow object to a Structure or FastStructure infile

## Usage

``` r
vcf2structure(vcf_arrow, out_file, keep_groups = NULL, method = "S")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- method:

  -\> flag for Structure/FastStructure formats, default 'S' (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external Structure or
FastStructure formatted file. Writing occurs in chunks whose size is
determined by the read_vcf() function. Larger chunks result in faster
writing speeds. If no groups are defined, the default behavior is to use
all groups. The flag parameter controls whether Structure (flag = 'S')
or FastStructure (flag = 'F') formatted output is written out.

## Author

Tomas Hrbek May 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2structure(vcf, out_file = tempfile(fileext = ".str"))
#> ℹ Accumulating Structure: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Structure file...
#> ✔ Structure file written to /tmp/RtmpodaNXU/file1ffe480ae0ad.str
vcf2structure(vcf, out_file = tempfile(fileext = ".fstr"),
              method = "F")
#> ℹ Accumulating Structure: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Structure file...
#> ✔ Structure file written to /tmp/RtmpodaNXU/file1ffe3418e73f.fstr
```
