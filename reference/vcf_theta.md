# vcf_theta

Calculates basic Watterson's theta and pi for all samples and for sample
groups from VCFArrow format data.

## Usage

``` r
vcf_theta(vcf_arrow, keep_groups = NULL)
```

## Arguments

- gt:

  -\> VCFArrow gt in long format

- sample:

  -\> samples to be included in analysis

- grps:

  -\> assignment of individuals to groups

## Value

list of statistics

## Details

This function calculates Watterson's theta and pi for the entire
VCFArrow object, and for groups of individuals whose grouping is
indicated by the groups slot in the VCFArrow object.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_theta(gt_matrix = my_gt_matrix, samples = vcf_arrow@samples, grps = vcf_arrow@groups)
#> Error in vcf_theta(gt_matrix = my_gt_matrix, samples = vcf_arrow@samples,     grps = vcf_arrow@groups): unused arguments (gt_matrix = my_gt_matrix, samples = vcf_arrow@samples, grps = vcf_arrow@groups)
vcf_theta(my_gt_matrix, vcf_arrow@samples, vcf_arrow@groups)
#> Error in vcf_theta(my_gt_matrix, vcf_arrow@samples, vcf_arrow@groups): unused argument (vcf_arrow@groups)
```
