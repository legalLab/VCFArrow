# vcf_sub_loci_stratified

Randomly subsets SNVs from a VCFArrow object

## Usage

``` r
vcf_sub_SNVs_stratified(vcf_arrow, n_SNVs = 1000, seed = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- n_SNVs:

  -\> number of SNVs to subset, default 1000 (integer)

- seed:

  -\> random number generator seed, default NULL (integer)

## Value

VCFArrow object

## Details

This function subsets a VCFArrow object to specific number of SNVs,
returning new VCFArrow object. The subsampling is stratified, i.e. the
same proportion of SNVs per CHROM. The seed for random number generator
is automatically generated unless specified.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_sub_SNVs_stratified(vcf, n_SNVs = 1000, seed = 42)
#> 
#> An object of class "VCFArrow"
#> 
#> Dimensions:
#>   Variants: 0 
#>   Samples:  18 
#> 
#> Quick stats:
#>   Non-missing variants: 0 
#> 
#> Phased genotypes: FALSE 
#> 
#> Storage:
#>   Path: /tmp/RtmpodaNXU/arrow_vcf_1ffe55cd2091 
#> 
#> Genotype storage (Arrow):
#> FileSystemDataset with 1 Feather file
#> 9 columns
#> .row_id: int32
#> sample: string
#> a1: int32
#> a2: int32
#> phased: bool
#> fmt: string
#> DP: double
#> GQ: double
#> ADR: double
#> 
#> See $metadata for additional Schema metadata
#> 
#> Variants: <empty>
#> 
#> Samples (first 5):
#> [1] "Pv120" "Pv126" "Pv13"  "Pv14"  "Pv27" 
#>   ... 13 more
```
