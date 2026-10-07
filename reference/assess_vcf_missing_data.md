# assess_vcf_missing_data

Quantifying missing data of all samples in VCF

## Usage

``` r
assess_vcf_missing_data(vcf_arrow, res_path, species, project, details = TRUE)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- res_path:

  -\> path to results (directory for output dataframe and plots)

- species:

  -\> sample name for plot (character)

- project:

  -\> project name / base output file name (character)

- details:

  -\> flag for adding project name into figure title, default TRUE
  (Boolean)

## Value

dataframe and plot of missing data for each sample in VCF

## Details

This function generates a dataframe of absolute and relative missing
data per sample, and a plot of relative missing data % per sample.

## Author

Tomas Hrbek April 2026

## Examples

``` r
f <- system.file("extdata", "vaillantii_discosnp_sub.vcf.gz",
                 package = "VCFArrow")
vcf <- read_vcf(f) |> set_vcf_groups(dirname(f))
#> ℹ VCF is being read in chunks of 50000 variants
#> ✔ VCF successfully read into a VCFArrow object
assess_vcf_missing_data(vcf, res_path = tempdir(),
                        species = "Phyllomedusa vaillantii",
                        project = "vaillantii")
```
