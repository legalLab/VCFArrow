# write_vcf

Write a VCFArrow object to an external VCF file

## Usage

``` r
write_vcf(vcf_arrow, out_file, gzip = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of the VCF file to be written to, no default (character)

- gzip:

  -\> a flag to GZIP VCF when writing, default FALSE (Boolean)

## Value

Invisibly returns the path of the written file.

## Details

This function writes a VCFArrow object to an external VCF file. Writing
occurs in chunks whose size is determined by the read_vcf() function.
Larger chunks result in faster writing speeds. It writes both
uncompressed and gz compressed files. Compressing increases writing time
be about 50%. For large files, it is recommended to output an
uncompressed VCF file, and then compress with GZIP or PIGZ.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> Warning: write_feather(version = 2) has been superseded by write_ipc_file().
#> Warning: `format = "feather"` is deprecated; use `format = "ipc"` instead.
#> ✔ VCF successfully read into a VCFArrow object
write_vcf(vcf, out_file = tempfile(fileext = ".vcf"))
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ VCF is being written in 1 chunk
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ✔ VCFArrow object successfully written to /tmp/Rtmpmnning/file367058807b01.vcf
write_vcf(vcf, out_file = tempfile(fileext = ".vcf.gz"),
          gzip = TRUE)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ VCF is being written in 1 chunk
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ✔ VCFArrow object successfully written to /tmp/Rtmpmnning/file367043314b90.vcf.gz
```
