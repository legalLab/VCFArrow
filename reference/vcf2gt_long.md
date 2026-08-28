# vcf2gt_long

Converts a VCFArrow object to tidy long format infile

## Usage

``` r
vcf2gt_long(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "gt_long",
  format = c("feather", "parquet", "csv"),
  col_select = NULL
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'gt_long' (character)

- format:

  -\> one of three output formats (arrow, parquet, CSV) (character)

- col_select:

  -\> optional selection of columns to save, default ALL

## Details

This function converts a VCFArrow object to an external SmartSNP
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. The
tidy data can be saved in either Arrow, Parquet or CSV formats. File
extension is added automatically if missing. Optionally, specific
columns can be saved, by default all columns are saved. The gt long slot
contains pre-calculated metrics in addition to just genotypes.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2gt_long(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "gt_long", format = "csv")
#> Error: object 'my_vcf' not found
vcf2gt_long(vcf_arrow, my_groups, format = "csv")
#> Error: object 'vcf_arrow' not found
vcf2gt_long(vcf_arrow, format = "csv")
#> Error: object 'vcf_arrow' not found
```
