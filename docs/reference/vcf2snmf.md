# vcf2snmf

Converts a VCFArrow object to an sNMF format infile

## Usage

``` r
vcf2snmf(vcf_arrow, keep_groups = NULL, out_file = "snmf_infile.geno")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'snmf_infile.geno' (character)

## Details

This function converts a VCFArrow object to an external sNMF formatted
file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2snmf(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "snmf_infile.geno")
#> Error: object 'my_vcf' not found
vcf2snmf(vcf_arrow, my_groups, out_file = "snmf_infile.geno")
#> Error: object 'vcf_arrow' not found
vcf2snmf(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
