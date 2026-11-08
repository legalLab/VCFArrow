# vcf2genepop

Converts a VCFArrow object to a Genepop format infile

## Usage

``` r
vcf2genepop(vcf_arrow, out_file, keep_groups = NULL, pop_names = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- out_file:

  -\> name of file to output, no default (character)

- keep_groups:

  -\> groups to retain, default NULL (character)

- pop_names:

  -\> write the group name after Pop on the line that starts each
  population ('pop WA') instead of Pop alone, default FALSE (Boolean)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing `out_file`.

## Details

This function converts a VCFArrow object to an external Genepop
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Each
population starts with a line holding the word Pop alone, as the GENEPOP
documentation shows it; GENEPOP itself names each population after its
last individual. Some readers, e.g. adegenet::read.genepop(), accept
only this form. With pop_names = TRUE, the line also holds the group
name ('pop WA'), which other programs may read as the population name.

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
vcf2genepop(vcf, out_file = tempfile(fileext = ".gen"))
#> ℹ Accumulating Genepop: 9313 variants x 18 samples (0 MiB raw storage)
#> Warning: 'arrow::read_feather' is deprecated.
#> Use 'read_ipc_file' instead.
#> See help("Deprecated")
#> ℹ Writing Genepop file...
#> ✔ Genepop file written to /tmp/RtmpVdrUAs/file315368e316f7.gen
```
