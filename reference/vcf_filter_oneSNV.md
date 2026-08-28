# vcf_filter_oneSNV

Subset a VCFArrow object keeping only 1 SNV per locus

## Usage

``` r
vcf_filter_oneSNV(vcf_arrow, block_size = 10000)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- block_size:

  -\> size of linked SNV blocks, default 10000 bp (integer)

## Value

subsetted VCFArrow object

## Details

This function subsets a VCFArrow object keeping only 1 SNV per locus,
returning a new VCFArrow object. Locus is defined as a different
chromosome or a block of the 'block_size' parameter value within a
chromosome. The first SNV independent of quality is taken (may modify
this in the future).

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_oneSNV(vcf_arrow = my_vcf, block_size = 10000)
#> Error: object 'my_vcf' not found
vcf_filter_oneSNV(my_vcf, 10000)
#> Error: object 'my_vcf' not found
vcf_filter_oneSNV(my_vcf)
#> Error: object 'my_vcf' not found
```
