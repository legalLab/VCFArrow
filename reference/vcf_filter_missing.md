# vcf_filter_missing

Remove samples with \> % missing data from a VCFArrow object

## Usage

``` r
vcf_filter_missing(vcf_arrow, threshold = 0.5, f_invar = TRUE, verbose = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal missing threshold, default 0.5 (numeric)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

- verbose:

  -\> report filtering stats, default TRUE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes samples from a VCFArrow object if they have above
threshold missing loci, returning a new VCFArrow object. By default will
remove any loci that may have become invariant as the result of the
removal of samples. By default will report removed samples, final %
missing data, and number of retained samples after sample filtering.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_missing(vcf_arrow = my_vcf, threshold = my_threshold)
#> Error: object 'my_vcf' not found
vcf_filter_missing(my_vcf, my_threshold)
#> Error: object 'my_vcf' not found
vcf_filter_missing(my_vcf)
#> Error: object 'my_vcf' not found
```
