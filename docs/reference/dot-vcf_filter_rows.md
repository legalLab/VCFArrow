# vcf_filter_rows

Unified API for filtering VCFArrow objects by row IDs.

## Usage

``` r
.vcf_filter_rows(vcf_arrow, keep)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep:

  -\> rows to keep (numeric or logical)

## Value

subsetted VCFArrow object

## Details

This function removes all rows in the 'keep' parameter from a VCFArrow
object, returning a new VCFArrow object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_rank(vcf_arrow = my_vcf, keep = rows_to_keep)
#> Error: object 'my_vcf' not found
vcf_filter_rank(my_vcf, rows_to_keep)
#> Error: object 'my_vcf' not found
```
