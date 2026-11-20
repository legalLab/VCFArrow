# vcf2fstructure

Converts a VCFArrow object to a FastStructure infile

## Usage

``` r
vcf2fstructure(vcf_arrow, out_file, keep_groups = NULL, method = "F")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- method:

  -\> flag for Structure/FastStructure formats, default 'F' (character)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function is an alias of vcf2structure() that writes the
FastStructure format by default (method = 'F'); see vcf2structure() for
the file layout. FastStructure adds the extension '.str' to the name
given with –input, so give out_file the extension .str (e.g. "data.str",
read with –input=data –format=str).

## Author

Tomas Hrbek November 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf2fstructure(vcf, out_file = tempfile(fileext = ".str"))
#> ℹ Accumulating Structure: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing Structure file...
#> ✔ Structure file written to /tmp/RtmpKz15ZV/file1d3a61f5490.str
```
