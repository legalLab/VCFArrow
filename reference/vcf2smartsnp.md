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
of writing `out_file` and its sample file.

## Details

This function converts a VCFArrow object to an external SmartSNP
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
genotype file has one line per SNP and one column per sample, separated
by spaces: the number of ALT alleles (0, 1, 2), 9 if missing. It has no
header line, as smartsnp requires ("without row names or column
headings"). The samples and their groups, in the order of the columns,
are written to a separate tab-separated file with a header line (sample,
group), named after out_file without its extension plus
'\_samples.smartsnp' (e.g. 'data.txt' gives 'data_samples.smartsnp').
Its group column is what smartsnp's analyses take as sample_group, e.g.
`s <- read.delim("data_samples.smartsnp")` then
`smartsnp::smart_pca("data.txt", sample_group = s$group)` or
`smartsnp::smart_permanova("data.txt", sample_group = s$group)`.

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
vcf2smartsnp(vcf, out_file = tempfile(fileext = ".txt"))
#> ℹ Building SmartSNP: 9313 variants x 18 samples (0 MiB raw storage)
#> ℹ Writing SmartSNP file...
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ✔ SmartSNP files written to /tmp/RtmpVdrUAs/file3153bf6b517.txt and /tmp/RtmpVdrUAs/file3153bf6b517_samples.smartsnp
```
