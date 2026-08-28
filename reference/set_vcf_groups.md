# set_vcf_groups

Associate samples with group info based on "strt" file

## Usage

``` r
set_vcf_groups(vcf_arrow, data_path, strt = "strata")
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- data_path:

  -\> path to data (where strata are located)

- strt:

  -\> file with 2+ columns, one id column and one pop column (tsv file,
  header = TRUE)

## Value

vcf_arrow object with groups slot filled

## Details

This function creates list of factors of sample-to-group assignments and
groups based on a strata and group file. The strt is a tsv file with 2+
columns, with 2 columns named 'id' and 'pop'. By default strt is named
'strata' and found in the datapath where the VCF file is found.

## Author

Tomas Hrbek April 2026

## Examples

``` r
get_vcf_group_info(vcf = my_vcf, data_path = my_data_path, strt = "strata")
#> Error in get_vcf_group_info(vcf = my_vcf, data_path = my_data_path, strt = "strata"): could not find function "get_vcf_group_info"
get_vcf_group_info(my_vcf, my_data_path, "strata")
#> Error in get_vcf_group_info(my_vcf, my_data_path, "strata"): could not find function "get_vcf_group_info"
get_vcf_group_info(my_vcf, my_data_path)
#> Error in get_vcf_group_info(my_vcf, my_data_path): could not find function "get_vcf_group_info"
```
