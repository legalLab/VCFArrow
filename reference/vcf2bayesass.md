# vcf2bayesass

Converts a VCFArrow object to BayesAss format infile

## Usage

``` r
vcf2bayesass(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "bayesass3_infile.immanc"
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'bayesass3_infile.immanc'
  (character)

## Details

This function converts a VCFArrow object to an external BayesAss
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2bayesass(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "bayesass3_infile.immanc")
#> Error: object 'my_vcf' not found
vcf2bayesass(vcf_arrow, my_groups, out_file = "bayesass3_infile.immanc")
#> Error: object 'vcf_arrow' not found
vcf2bayesass(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
