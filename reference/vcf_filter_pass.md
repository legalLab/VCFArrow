# vcf_filter_pass

Remove loci that did not PASS FILTER in a VCFArrow object

## Usage

``` r
vcf_filter_pass(vcf_arrow)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

## Value

subsetted VCFArrow object

## Details

This function removes loci that did not PASS the FILTER in a VCFArrow
object, returning a new VCFArrow object. PASS indicates that a variant
has successfully passed all applied quality control filters. When no
filters were applied and the FILTER field is '.' (Dot),
vcf_filter_pass() defaults to passing the locus.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_pass(vcf_arrow = my_vcf)
#> Error: object 'my_vcf' not found
vcf_filter_pass(my_vcf)
#> Error: object 'my_vcf' not found
```
