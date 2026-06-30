# vcf_filter_hets

Remove loci above a heterozigosity threshold from a VCFArrow object

## Usage

``` r
vcf_filter_hets(vcf_arrow, threshold = 0.5)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> heterozigosity threshold, default 0.5 (numeric)

## Value

subsetted VCFArrow object

## Details

This function removes loci above a heterozigosity threshold from a
VCFArrow object, returning a new VCFArrow object. High heterozigosities
are indicative of potential paralogs.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_hets(vcf_arrow = my_vcf, threshold = 0.5)
#> Error: object 'my_vcf' not found
vcf_filter_hets(my_vcf, 0.5)
#> Error: object 'my_vcf' not found
vcf_filter_hets(my_vcf)
#> Error: object 'my_vcf' not found
```
