# vcf_bind

Bind two or more VCFArrow objects into a new VCFArrow object

## Usage

``` r
vcf_bind(..., .check = TRUE)
```

## Arguments

- ...:

  -\> a collection of VCFArrow objects

- .check:

  -\> check that VCFArrow objects are compatible, default TRUE (Boolean)

## Value

VCFArrow object

## Details

This function binds two or more VCFArrow objects, returning new VCFArrow
object. The VCFArrow objects must have the same SNPs, and must have
unique individuals. Default is to check compatibility between VCFArrow
objects; this default can but is not meant to be changed. If VCFArrow
objects are not compatible but need to be merged, use vcf_bind_sparse().

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_bind(my_vcf1, my_vcf2, other_vcf, ..., check = TRUE)
#> Error: '...' used in an incorrect context
vcf_bind(my_vcf1, my_vcf2, other_vcf)
#> Error: object 'my_vcf1' not found
```
