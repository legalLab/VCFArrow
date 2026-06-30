# vcf_filter_adr

Correct or remove genotypes with \> % ALT/REF ratio from a VCFArrow
object

## Usage

``` r
vcf_filter_adr(
  vcf_arrow,
  mode = c("correct", "remove"),
  threshold = 0.1,
  f_invar = TRUE
)
```

## Arguments

- vcf_arrow:

  -\> VCFArrow object

- mode:

  -\> switch between 'correct' or 'remove' mode (character)

- threshold:

  -\> decimal missing threshold, default 0.1 (numeric)

- f_invar:

  -\> filter invariant loci flag, default TRUE (Boolean)

## Value

subsetted VCFArrow object

## Details

This function either changes genotypes or makes genotypes missing in a
VCFArrow object if they have above/below threshold normalized ADR ratio,
returning a new VCFArrow object. The ADR is calculated as ADR = ALT /
(REF + ALT). If ADR \> threshold, the genotype becomes ALT homozygous.
If ADR \< threshold, the genotype becomes REF homozygous. By default
will remove any loci that may have become invariant as the result of the
change/removal of genotypes

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_filter_adr(vcf_arrow = my_vcf, mode = "correct", threshold = my_threshold, f_invar = TRUE)
#> Error: object 'my_vcf' not found
vcf_filter_adr(my_vcf, "correct", my_threshold, TRUE)
#> Error: object 'my_vcf' not found
vcf_filter_adr(my_vcf, "correct")
#> Error: object 'my_vcf' not found
```
