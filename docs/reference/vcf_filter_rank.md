# vcf_filter_rank

Remove loci below RANK threshold from a VCFArrow object

## Usage

``` r
vcf_filter_rank(vcf_arrow, threshold = 0.4, keep_na = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal rank threshold, default 0.4 (numeric)

## Value

subsetted VCFArrow object

## Details

This function removes loci below RANK threshold from a VCFArrow object,
returning a new VCFArrow object. RANK is calculated in DiscoSNP-RAD
(Gauthier et. al. 2020) and registered as Pk in INFO. RANK is calculated
as sqrt(chi-sqr/n) of allele read counts, and used for paralog detection
-\> very low rank values (\<0.4) are indicative of paralogs.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_rank(vcf_arrow = my_vcf, rank = 0.4)
#> Error in vcf_filter_rank(vcf_arrow = my_vcf, rank = 0.4): unused argument (rank = 0.4)
vcf_filter_rank(my_vcf, 0.4)
#> Error: object 'my_vcf' not found
vcf_filter_rank(my_vcf)
#> Error: object 'my_vcf' not found
```
