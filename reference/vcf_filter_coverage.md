# vcf_filter_coverage

Remove genotypes below read DP threshold from a VCFArrow object

## Usage

``` r
vcf_filter_coverage(vcf_arrow, threshold = 10)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> DP threshold, default 10 (integer)

## Value

subsetted VCFArrow object

## Details

This function removes genotypes below a DP threshold from a VCFArrow
object, returning a new VCFArrow object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_coverage(vcf_arrow = my_vcf, threshold = 10)
#> Error: object 'my_vcf' not found
vcf_filter_coverage(my_vcf, 10)
#> Error: object 'my_vcf' not found
vcf_filter_coverage(my_vcf)
#> Error: object 'my_vcf' not found
```
