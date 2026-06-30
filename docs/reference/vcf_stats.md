# vcf_stats

Calculates basic stats of each samples from VCFArrow format data.
Includes average read depth per individual, missing data per individual,
Watterson's theta and pi.

## Usage

``` r
vcf_stats(vcf_arrow, res_path, project, theta = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- res_path:

  -\> directory where to write results

- project:

  -\> base name of the project file

- theta:

  -\> flag to perform theta and pi calculation, default FALSE (Boolean)

## Value

table of statistics

## Details

This function calculates average read depth, heterozygosity number of
heterozygotes, number of reference and alternative homozygotes, missing
data and total number SNPs of each sample in an VCFArrow object.
Optionally calls vcf_theta() to get total and group Watterson's theta
and pi.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_stats(vcf_arrow = my_vcf, res_path = my_res_path, project = my_project, theta = FALSE)
#> Error: object 'my_vcf' not found
vcf_stats(my_vcf, my_res_path, my_project, FALSE)
#> Error: object 'my_vcf' not found
vcf_stats(my_vcf, my_res_path, my_project)
#> Error: object 'my_vcf' not found
```
