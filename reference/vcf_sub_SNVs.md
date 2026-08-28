# vcf_sub_loci

Randomly subsets SNVs from a VCFArrow object

## Usage

``` r
vcf_sub_SNVs(vcf_arrow, n_SNVs = 10000, seed = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- seed:

  -\> random number generator seed, default NULL (Boolean)

- n_loci:

  -\> number of SNVs to subset, default 10000 (integer)

## Value

VCFArrow object

## Details

This function subsets a VCFArrow object to specific number of SNVs,
returning new VCFArrow object. The seed for random number generator is
automatically generated unless otherwise specified.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_sub_SNVs(vcf_arrow = my_vcf, n_SNVs = n_SNVs, seed = my_seed)
#> Error: object 'my_vcf' not found
vcf_sub_SNVs(my_vcf, n_SNVs, 42)
#> Error: object 'my_vcf' not found
vcf_sub_SNVs(my_vcf)
#> Error: object 'my_vcf' not found
```
