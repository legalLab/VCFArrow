# vcf_filter_quality

Remove loci below QUAL threshold from a VCFArrow object

## Usage

``` r
vcf_filter_quality(vcf_arrow, threshold = 30)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> QUAL threshold, default 30 (integer)

## Value

subsetted VCFArrow object

## Details

This function removes loci below QUAL threshold from a VCFArrow object,
returning a new VCFArrow object. QUAL is a Phred-scaled quality metric
assessing whether any variant exists at that position across all
samples. QUAL is different from PL which represents the confidence in a
genotype call of a specific sample. When QUAL not reported '.' (Dot),
QUAL is assumed to be above threshold.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_quality(vcf_arrow = my_vcf, threshold = 30)
#> Error: object 'my_vcf' not found
vcf_filter_quality(my_vcf, 30)
#> Error: object 'my_vcf' not found
vcf_filter_quality(my_vcf)
#> Error: object 'my_vcf' not found
```
