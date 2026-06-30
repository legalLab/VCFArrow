# vcf_filter_indels

Remove non-indels from a VCFArrow object

## Usage

``` r
vcf_filter_indels(vcf_arrow)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

## Value

subsetted VCFArrow object

## Details

This function removes non-indel loci from a VCFArrow object, returning a
new VCFArrow object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_indels(vcf_arrow = my_vcf)
#> Error: object 'my_vcf' not found
vcf_filter_indels(my_vcf)
#> Error: object 'my_vcf' not found
```
