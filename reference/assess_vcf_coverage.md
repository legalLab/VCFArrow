# assess_vcf_coverage

Quantifying read depth of all samples in VCF

## Usage

``` r
assess_vcf_coverage(
  vcf_arrow,
  res_path,
  species,
  project,
  details = TRUE,
  max_points_per_sample = 5000L
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- res_path:

  -\> path to results (directory for output dataframe and plots)

- species:

  -\> species name, the plot title (character)

- project:

  -\> project name / base output file name (character)

- details:

  -\> flag for adding project name into figure title, default TRUE
  (Boolean)

- max_points_per_sample:

  -\> maximum number of SNVs from which to plot coverage, default 5000
  (integer)

## Value

Invisibly returns the input VCFArrow object; called for its side effect
of writing a violin plot of the read depth of each sample to
`<project>_coverage.pdf`, `.svg` and `.png` in `res_path`.

## Details

This function generates a violin plot of read depths (DP) per sample,
from the sample's called genotypes; the DP of missing genotypes is not
used. Each violin ends at the sample's lowest and highest plotted depth.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
assess_vcf_coverage(vcf, res_path = tempdir(),
                    species = "Phyllomedusa vaillantii",
                    project = "vaillantii")
```
