# vcf_copy

Copy a VCFArrow object

## Usage

``` r
vcf_copy(vcf_arrow)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

## Value

VCFArrow object

## Details

This function copies a VCFArrow object. It makes an actual physical copy
of the object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_copy(vcf_arrow = my_vcf)
#> Error: object 'my_vcf' not found
vcf_copy(my_vcf)
#> Error: object 'my_vcf' not found
```
