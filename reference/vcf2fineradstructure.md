# vcf2fineradstructure

Converts a VCFArrow object to FineRadStructure format infile

## Usage

``` r
vcf2fineradstructure(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "fineradstructure_infile.txt"
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'fineradstructure_infile.txt'
  (character)

## Details

This function converts a VCFArrow object to an external Treemix
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2fineradstructure(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "fineradstructure_infile.txt")
#> Error: object 'my_vcf' not found
vcf2fineradstructure(vcf_arrow, my_groups, out_file = "fineradstructure_infile.txt")
#> Error: object 'vcf_arrow' not found
vcf2fineradstructure(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
