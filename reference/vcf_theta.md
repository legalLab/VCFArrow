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
indicated by the groups slot in the VCFArrow object. Only biallelic SNVs
are used, and a site is used for a group (valid) when at least two of
its samples are called. With n called alleles and ALT allele frequency p
at a site, the site's nucleotide diversity is n / (n - 1) \* 2p(1 - p),
the proportion of differing pairs among all pairs of called alleles (the
unbiased estimator, as in vcftools –site-pi), and its Watterson
contribution is 1 / a1 if the site is segregating, with a1 = sum(1 / i)
for i = 1 .. n - 1. pi and theta are the means of these values over the
valid sites. The statistics are accumulated chunk by chunk as running
sums per group, so memory use does not grow with the number of variants,
and chunks are processed in parallel when workers are set with
vcf_set_workers().

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
#> [1] 0.3415246
#> 
#> $theta_w
#> [1] 0.3225428
#> 
#> $theta_g
#>    sample group n_ind   theta_w        pi
#> 1   Pv120    WA     5 0.2113146 0.2231487
#> 2   Pv126    GS     5 0.1916500 0.2045097
#> 3    Pv13    BS     5 0.1842761 0.1979007
#> 4    Pv14    GS     5 0.1916500 0.2045097
#> 5    Pv27    GS     5 0.1916500 0.2045097
#> 6    Pv28    GS     5 0.1916500 0.2045097
#> 7     Pv2    WA     5 0.2113146 0.2231487
#> 8    Pv31    WA     5 0.2113146 0.2231487
#> 9    Pv56    WA     5 0.2113146 0.2231487
#> 10   Pv62    BS     5 0.1842761 0.1979007
#> 11   Pv68    BS     5 0.1842761 0.1979007
#> 12   Pv73    BS     5 0.1842761 0.1979007
#> 13   Pv78    GS     5 0.1916500 0.2045097
#> 14   Pv79    BS     5 0.1842761 0.1979007
#> 15   Pv93    WA     5 0.2113146 0.2231487
#> 16  Pb2Jp    OG     3 0.1712101 0.1854309
#> 17 Pb2Scx    OG     3 0.1712101 0.1854309
#> 18  Pb1Rd    OG     3 0.1712101 0.1854309
#> 
```
