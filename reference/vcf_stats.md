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
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_stats(vcf, res_path = tempdir(), project = "vaillantii",
          theta = TRUE)
#> ℹ Computing per-sample stats: 10000 variants x 18 samples, reading 1 chunk directly
```
