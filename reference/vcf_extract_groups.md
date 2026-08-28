# vcf_extract_groups

Extract groups of samples from a VCFArrow object

## Usage

``` r
vcf_extract_groups(
  vcf_arrow,
  groups,
  keep = TRUE,
  f_invar = TRUE,
  verbose = TRUE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- groups:

  -\> groups to retain/drop (character)

- keep:

  -\> retain (TRUE) or drop (FALSE) individuals flag, default TRUE
  (Boolean)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report filtering stats, default TRUE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes groups of samples from a VCFArrow object,
returning a new VCFArrow object. It uses the 'keep' flag to either keep
or drop the groups in the list. By default will remove any loci that may
have become invariant as the result of the removal of samples. By
default will report removed samples, final % missing data, and number of
retained samples after sample filtering.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_extract_groups(vcf_arrow = my_vcf, groups = my_groups, keep = TRUE, f_invar = TRUE, verbose = TRUE)
#> Error: object 'my_vcf' not found
vcf_extract_groups(my_vcf, my_groups, TRUE, TRUE, TRUE)
#> Error: object 'my_vcf' not found
vcf_extract_groups(my_vcf, my_groups)
#> Error: object 'my_vcf' not found
```
