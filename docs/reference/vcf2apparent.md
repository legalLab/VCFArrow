# vcf2apparent

Converts a VCFArrow object to an Apparent format infile

## Usage

``` r
vcf2apparent(
  vcf_arrow,
  keep_groups = NULL,
  key = "All",
  out_file = "apparent_infile.txt"
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- key:

  -\> relationship type (All, Pa, Mo, Fa, Off), default All (character)

- out_file:

  -\> name of file to output, default 'apparent_infile.txt' (character)

## Details

This function converts a VCFArrow object to an external SmartSNP
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.
Possible relationships defined by the parameter 'kee' are All, Pa, Mo,
Fa, Off.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2apparent(vcf_arrow = my_vcf, keep_groups = my_groups, key = my_key, out_file = "apparent_infile.txt")
#> Error: object 'my_vcf' not found
vcf2apparent(vcf_arrow, my_groups, my_key, out_file = "apparent_infile.txt")
#> Error: object 'vcf_arrow' not found
vcf2apparent(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
