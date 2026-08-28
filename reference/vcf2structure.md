# vcf2structure

Converts a VCFArrow object to a Structure or FastStructure infile

## Usage

``` r
vcf2structure(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "structure.str",
  method = "S"
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'structure.str' (character)

- method:

  -\> flag for Structure/FastStructure formats, default 'S' (character)

## Details

This function converts a VCFArrow object to an external Structure or
FastStructure formatted file. Writing occurs in chunks whose size is
determined by the read_vcf() function. Larger chunks result in faster
writing speeds. If no groups are defined, the default behavior is to use
all groups. The flag parameter controls whether Structure (flag = 'S')
or FastStructure (flag = 'F') formatted output is written out.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2structure(vcf_arrow = my_vcf, keep_groups = NULL, out_file = "structure.str", method = "S")
#> Error: object 'my_vcf' not found
vcf2structure(my_vcf, keep_groups, out_file = "structure.str")
#> Error: object 'my_vcf' not found
vcf2structure(my_vcf)
#> Error: object 'my_vcf' not found
```
