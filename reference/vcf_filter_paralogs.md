# vcf_filter_paralogs

Remove loci above REPEAT and below RANK threshold from a VCFArrow object

## Usage

``` r
vcf_filter_paralogs(
  vcf_arrow,
  threshold_rk = 0.2,
  threshold_rpt = 0.5,
  keep_na = FALSE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold_rk:

  -\> decimal RANK threshold, loci below it are removed, default 0.2
  (numeric)

- threshold_rpt:

  -\> decimal REPEAT threshold, loci above it are removed, default 0.5
  (numeric)

- keep_na:

  -\> retain loci without a RANK or REPEAT value, default FALSE
  (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes loci below the RANK threshold and above the REPEAT
threshold from a VCFArrow object, returning a new VCFArrow object. It
applies vcf_filter_rank() and then vcf_filter_rpt(). RANK and REPEAT
metrics focus on different paralog signals, and used jointly maximize
true positive and minimize false positive detections. Simulations
indicate that removing loci with Rk \< 0.2 or RPT \> 0.5 gives the best
results, with Rk \< 0.4 or RPT \> 0.5 being slightly more conservative.

## Author

Tomas Hrbek September 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_paralogs(vcf, threshold_rk = 0.2, threshold_rpt = 0.5)
#> Error in vcf_filter_paralogs(vcf, threshold_rk = 0.2, threshold_rpt = 0.5): could not find function "vcf_filter_paralogs"
```
