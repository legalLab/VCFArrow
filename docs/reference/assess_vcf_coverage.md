# assess_vcf_coverage

Quantifying read depth of all samples in VCF Inspired by
https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html

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

  -\> sample name for plot (character)

- project:

  -\> project name / base output file name (character)

- details:

  -\> flag for adding project name into figure title, default TRUE
  (Boolean)

- max_points_per_sample:

  -\> maximum number of SNVs from which to plot coverage, default 5000
  (integer)

## Value

violin plot of read depth for each sample in VCF

## Details

This function generates a violin plot of read depths per sample.

## Author

Tomas Hrbek April 2026

## Examples

``` r
assess_vcf_coverage(vcf_arrow = my_vcf, res_path = my_res_path, species = species_name, project = project_name)
#> Error: object 'my_vcf' not found
assess_vcf_coverage(my_vcf, my_res_path, species_name, project_name, details = TRUE, max_points_per_sample = 5000L)
#> Error: object 'my_vcf' not found
assess_vcf_coverage(my_vcf, my_res_path, species_name, project_name)
#> Error: object 'my_vcf' not found
```
