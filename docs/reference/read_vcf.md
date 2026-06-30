# read_vcf

Read VCF file and store its content within a VCFArrow object

## Usage

``` r
read_vcf(vcf_file, chunk_size = 50000)
```

## Arguments

- vcf_file:

  -\> VCF file

- chunk_size:

  -\> number of variants to read in at a time, default 50000 (integer)

## Value

VCFArrow object

## Details

This function read a VCF file into an S4 class object in chunks,
returning a VCFArrow object. It accepts both uncompressed and gz
compressed files. The GT field is stored as an Apache Arrow in Long
format. Various metrics are precalculated for fast and easy filtering.
GT slot content stored in a TEMP directory for lazy loading.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_rank(vcf_file = my_vcf, chunk_size = 50000)
#> Error in vcf_filter_rank(vcf_file = my_vcf, chunk_size = 50000): unused arguments (vcf_file = my_vcf, chunk_size = 50000)
vcf_filter_rank(my_vcf, 50000)
#> Error: object 'my_vcf' not found
vcf_filter_rank(my_vcf)
#> Error: object 'my_vcf' not found
```
