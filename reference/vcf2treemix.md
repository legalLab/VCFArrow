# vcf2treemix

Converts a VCFArrow object to Treemix format infile

## Usage

``` r
vcf2treemix(vcf_arrow, keep_groups = NULL, out_file = "treemix_infile.txt")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'treemix_infile.txt' (character)

## Details

This function converts a VCFArrow object to an external Treemix
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2treemix(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "treemix_infile.txt")
#> Error: object 'my_vcf' not found
vcf2treemix(vcf_arrow, my_groups, out_file = "treemix_infile.txt")
#> Error: object 'vcf_arrow' not found
vcf2treemix(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
