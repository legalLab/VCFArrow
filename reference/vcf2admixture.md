# vcf2admixture

Converts a VCFArrow object to a PLINK .bed format infile plus
ADMIXTURE-specific .pop file

## Usage

``` r
vcf2admixture(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "admixture_in",
  sex = NULL,
  pheno = NULL,
  supervised = FALSE,
  reference_groups = NULL
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- keep_groups:

  -\> groups to retain, default NULL (character)

- out_file:

  -\> name of file to output, default 'plink_out' (character)

- sex:

  -\> vector of sexes of samples, default NULL (character)

- pheno:

  -\> vector of phenotypes of samples, default NULL (character)

- supervised:

  -\> flag to generate .pop for use in ADMIXTURE's supervised mode,
  default FALSE (Boolean)

- reference_groups:

  -\> vector of ancestry group labels, default NULL (character)

## Details

This function converts a VCFArrow object to an external PLINK .bed
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Sex
and phenotype vectors are optional. If not defined sex = 0, pheno = -9.
For supervised analyses, the supervised flag needs to be true to
generate .pop file needed for ADMIXTURE's supervised mode. The
reference_groups vector specifies ancestral information of individuals.
When NULL, or for any missing individual, ancestry is inferred.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2admixture(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "admix_out", sex = sex, pheno = pheno, supervised = TRUE, reference_groups = reference_groups)
#> Error: object 'my_vcf' not found
vcf2admixture(vcf_arrow, my_groups, out_file = "admix_out", supervised = TRUE, reference_groups = c("XXX", "YYY"))
#> Error: object 'vcf_arrow' not found
vcf2admixture(vcf_arrow, my_groups, out_file = "admix_out", supervised = TRUE)
#> Error: object 'vcf_arrow' not found
vcf2admixture(vcf_arrow, my_groups, out_file = "admix_out")
#> Error: object 'vcf_arrow' not found
vcf2admixture(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
