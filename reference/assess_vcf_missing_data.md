# assess_vcf_missing_data

Quantifying missing data of all samples in VCF Inspired by
https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html

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
assess_vcf_missing_data(vcf_arrow = my_vcf, res_path = my_res_path, species = species_name, project = project_name)
#> Error: object 'my_vcf' not found
assess_vcf_missing_data(my_vcf, my_res_path, species_name, project_name, details = TRUE)
#> Error: object 'my_vcf' not found
assess_vcf_missing_data(my_vcf, my_res_path, species_name, project_name)
#> Error: object 'my_vcf' not found
```
