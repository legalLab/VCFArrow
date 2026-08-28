# vcf_filter_maf

Remove loci below MAF threshold from a VCFArrow object

## Usage

``` r
vcf_filter_maf(vcf_arrow, threshold = 0.05)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal MAF threshold, default 0.05 (numeric)

## Value

subsetted VCFArrow object

## Details

This function removes loci below MAF threshold from a VCFArrow object,
returning a new VCFArrow object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_maf(vcf_arrow = my_vcf, threshold = 0.05)
#> Error: object 'my_vcf' not found
vcf_filter_maf(my_vcf, 0.05)
#> Error: object 'my_vcf' not found
vcf_filter_maf(my_vcf)
#> Error: object 'my_vcf' not found
```
