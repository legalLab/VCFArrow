# vcf_theta

Calculates basic Watterson's theta and pi for all samples and for sample
groups from VCFArrow format data.

## Usage

``` r
vcf_theta(vcf_arrow, keep_groups = NULL)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

## Value

list of statistics

## Details

This function calculates Watterson's theta and pi for the entire
VCFArrow object, and for groups of individuals whose grouping is
indicated by the groups slot in the VCFArrow object. The statistics are
accumulated chunk by chunk as running sums per group, so memory use does
not grow with the number of variants, and chunks are processed in
parallel when workers are set with vcf_set_workers().

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_theta(vcf)
#> $pi
#> [1] 0.3065724
#> 
#> $theta_w
#> [1] 0.3225428
#> 
#> $theta_g
#>    sample group n_ind   theta_w        pi
#> 1   Pv120    WA     5 0.2113146 0.1897392
#> 2   Pv126    GS     5 0.1916500 0.1751468
#> 3    Pv13    BS     5 0.1842761 0.1692308
#> 4    Pv14    GS     5 0.1916500 0.1751468
#> 5    Pv27    GS     5 0.1916500 0.1751468
#> 6    Pv28    GS     5 0.1916500 0.1751468
#> 7     Pv2    WA     5 0.2113146 0.1897392
#> 8    Pv31    WA     5 0.2113146 0.1897392
#> 9    Pv56    WA     5 0.2113146 0.1897392
#> 10   Pv62    BS     5 0.1842761 0.1692308
#> 11   Pv68    BS     5 0.1842761 0.1692308
#> 12   Pv73    BS     5 0.1842761 0.1692308
#> 13   Pv78    GS     5 0.1916500 0.1751468
#> 14   Pv79    BS     5 0.1842761 0.1692308
#> 15   Pv93    WA     5 0.2113146 0.1897392
#> 16  Pb2Jp    OG     3 0.1712101 0.1498553
#> 17 Pb2Scx    OG     3 0.1712101 0.1498553
#> 18  Pb1Rd    OG     3 0.1712101 0.1498553
#> 
```
