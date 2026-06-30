# vcf_filter_multiSNV

Subset a VCFArrow object keeping only loci with 2+ SNVs per locus

## Usage

``` r
vcf_filter_multiSNV(vcf_arrow, block_size = 10000, minSNV = 2, maxSNV = 5)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- block_size:

  -\> size of linked SNV blocks, default 10000 bp (integer)

- minSNV:

  -\> minimum linked block size, default 2 (integer)

- maxSNV:

  -\> maximum number of selected linked SNVs per block, default 5
  (integer)

## Value

subsetted VCFArrow object

## Details

This function subsets a VCFArrow object keeping only loci with between
min and max \# of SNVs per locus, returning a new VCFArrow object.
Default min = 2 and max = 5 SNVs per locus (recommended as input for
fineRADstructure analyses). Locus is defined as a different chromosome
or a block of the 'block_size' parameter value within a chromosome.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_multiSNV(vcf_arrow = my_vcf, block_size = 10000, minSNV = 2, maxSNV = 5)
#> Error: object 'my_vcf' not found
vcf_filter_multiSNV(my_vcf, 10000, 2, 5)
#> Error: object 'my_vcf' not found
vcf_filter_multiSNV(my_vcf)
#> Error: object 'my_vcf' not found
```
