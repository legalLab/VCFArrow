# vcf_filter_invariant

Remove invariant loci from a VCFArrow object

## Usage

``` r
vcf_filter_invariant(vcf_arrow)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

## Value

subsetted VCFArrow object

## Details

This function removes invariant loci from a VCFArrow object, returning a
new VCFArrow object. This might be desirable after subsetting a VCFArrow
object by individuals.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_invariant(vcf_arrow = my_vcf)
#> Error: object 'my_vcf' not found
vcf_filter_invariant(my_vcf)
#> Error: object 'my_vcf' not found
```
