# vcf_filter_missingness

Remove loci above missingness threshold from a VCFArrow object

## Usage

``` r
vcf_filter_missingness(vcf_arrow, threshold = 0.1, verbose = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal missingness threshold, default 0.1 (numeric)

- verbose:

  -\> flag to report total % missing data after filtering, default TRUE
  (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes loci above missingness threshold from a VCFArrow
object, returning a new VCFArrow object. Missingness is a locus focused
metric, i.e. missing data per locus.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_missingness(vcf_arrow = my_vcf, threshold = p_miss)
#> Error: object 'my_vcf' not found
vcf_filter_missingness(my_vcf, p_miss, verbose = FALSE)
#> Error: object 'my_vcf' not found
vcf_filter_missingness(my_vcf, p_miss)
#> Error: object 'my_vcf' not found
```
