# vcf_filter_rpt

Remove loci above REPEAT threshold from a VCFArrow object

## Usage

``` r
vcf_filter_rpt(vcf_arrow, threshold = 0.5, keep_na = FALSE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- threshold:

  -\> decimal REPEAT threshold, loci above it are removed, default 0.5
  (numeric)

- keep_na:

  -\> retain loci without a REPEAT value, default FALSE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function removes loci above REPEAT threshold from a VCFArrow
object, returning a new VCFArrow object. REPEAT is calculated by
disco_haplotypes in DiscoSNP-RAD (Gauthier et. al. 2020) and registered
as RPT in INFO. REPEAT is calculated as max of overlap, depth score and
is 1 if in a giant-component bubble, and is used for paralog detection
-\> high repeat values (\>0.5) are indicative of paralogs.

## Author

Tomas Hrbek September 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
vcf_filter_rpt(vcf, threshold = 0.5)
#> Error in vcf_filter_rpt(vcf, threshold = 0.5): could not find function "vcf_filter_rpt"
```
