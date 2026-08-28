# vcf_bind_sparse

Bind two or more VCFArrow objects into a new VCFArrow object

## Usage

``` r
vcf_bind_sparse(
  ...,
  mode = c("intersect", "union"),
  absent_as = c("missing", "hom_ref")
)
```

## Arguments

- ...:

  -\> a collection of VCFArrow objects

## Value

VCFArrow object

## Details

This function binds two or more VCFArrow objects, returning new VCFArrow
object. The VCFArrow objects need not have the same SNPs, and must have
unique individuals.

## Author

Tomas Hrbek April 2026

## Examples

``` r
vcf_bind_sparse(my_vcf1, my_vcf2, other_vcf, ...)
#> Error: '...' used in an incorrect context
```
