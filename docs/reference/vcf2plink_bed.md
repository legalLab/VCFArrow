# vcf2plink_bed

Converts a VCFArrow object to a PLINK .bed format infile

## Usage

``` r
vcf2plink_bed(
  vcf_arrow,
  keep_groups = NULL,
  out_file = "plink_out",
  sex = NULL,
  pheno = NULL
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

## Details

This function converts a VCFArrow object to an external PLINK .bed
formatted file. Writing occurs in chunks whose size is determined by the
read_vcf() function. Larger chunks result in faster writing speeds. If
no groups are defined, the default behavior is to use all groups. Sex
and phenotype vectors are optional. If not defined sex = 0, pheno = -9.

## Author

Tomas Hrbek May 2026

## Examples

``` r
vcf2plink_bed(vcf_arrow = my_vcf, keep_groups = my_groups, out_file = "plink_out", sex = sex, pheno = pheno)
#> Error: object 'my_vcf' not found
vcf2plink_bed(vcf_arrow, my_groups, out_file = "plink_out")
#> Error: object 'vcf_arrow' not found
vcf2plink_bed(vcf_arrow)
#> Error: object 'vcf_arrow' not found
```
