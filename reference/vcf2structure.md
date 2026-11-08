# vcf2structure

Converts a VCFArrow object to a Structure or FastStructure infile

## Usage

``` r
vcf2structure(vcf_arrow, out_file, keep_groups = NULL, method = c("S", "F"))
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
or FastStructure (flag = 'F') formatted output is written out. Each
sample is written on two lines (one per allele copy): its name, its
group number and, for FastStructure, four columns of 0 (FastStructure
ignores the first six columns), then one allele per SNP: 0 (REF), 1
(ALT) or -9 (missing). FastStructure adds the extension '.str' to the
name given with –input, so give out_file the extension .str (e.g.
"data.str", read with –input=data –format=str).

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
vcf2structure(vcf, out_file = tempfile(fileext = ".str"))
#> ℹ Accumulating Structure: 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing Structure file...
#> ✔ Structure file written to /tmp/RtmpVdrUAs/file31535c7a1ff1.str
vcf2structure(vcf, out_file = tempfile(fileext = ".str"),
              method = "F")
#> ℹ Accumulating Structure: 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing Structure file...
#> ✔ Structure file written to /tmp/RtmpVdrUAs/file31535d1ee9cf.str
```
