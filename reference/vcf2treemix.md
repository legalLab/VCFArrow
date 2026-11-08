# vcf2treemix

Converts a VCFArrow object to Treemix format infile

## Usage

``` r
vcf2treemix(
  vcf_arrow,
  out_file,
  keep_groups = NULL,
  gzip = grepl("\\.gz$", out_file)
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- gzip:

  -\> compress the file with gzip, as TreeMix expects, default TRUE if
  out_file ends in ".gz" (Boolean)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external Treemix
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
first line holds the group names; each further line one SNP, with the
REF and ALT allele counts of each group as "ref,alt", separated by
spaces. TreeMix reads gzip-compressed input (treemix -i file.gz): name
out_file with the extension .gz to compress it.

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
vcf2treemix(vcf, out_file = tempfile(fileext = ".gz"))
#> ℹ Building Treemix: 9313 variants x 4 pops (0 MiB raw storage)
#> ℹ Writing Treemix file...
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ✔ Treemix file written to /tmp/RtmpVdrUAs/file3153352e6a6d.gz
```
