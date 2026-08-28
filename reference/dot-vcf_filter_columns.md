# vcf_filter_columns

Unified API for filtering VCFArrow objects by row IDs.

## Usage

``` r
.vcf_filter_columns(vcf_arrow, keep, f_invar = TRUE, verbose = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep:

  -\> rows to keep (numeric or logical)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report removed samples and final % missing data, default TRUE
  (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes all rows in the 'keep' parameter from a VCFArrow
object, returning a new VCFArrow object. Optionally will remove any loci
that may have become invariant as the result of the removal of samples.
Optionally will report removed samples, final % missing data, and number
of retained samples after sample filtering.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_rank(vcf_arrow = my_vcf, keep = rows_to_keep)
#> Error: object 'my_vcf' not found
vcf_filter_rank(my_vcf, rows_to_keep)
#> Error: object 'my_vcf' not found
```
