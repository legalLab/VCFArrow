# vcf2genlight

Converts a VCFArrow object to Genlight format infile

## Usage

``` r
vcf2genlight(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "genlight.rds",
  ploidy = 2L,
  save = FALSE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'genlight.rds' (character)

- ploidy:

  -\> ploidy level, default = 2 (integer)

- save:

  -\> save as R data object, default = FALSE (Boolean)

## Details

This function converts a VCFArrow object to an external Genlight
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.
Genlight objects can encode polyploid genomes, by default diploid
genomes are assumed. Genlight objects are in memory S4 objects, thus are
returned as such, but but optionally may be saved as R data objects.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2genlight(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "genlight.rds")
#> Error: object 'my_vcf' not found
vcf2genlight(vcf_arrow, my_groups)
#> Error: object 'vcf_arrow' not found
vcf2genlight(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
