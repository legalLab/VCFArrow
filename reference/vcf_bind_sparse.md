# vcf_bind_sparse

Deprecated: use vcf_bind(), which now binds objects with different
variants (see its mode and absent_as arguments). vcf_bind_sparse(...) is
vcf_bind(..., recover_loci = TRUE).

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

- mode:

  -\> see vcf_bind()

- absent_as:

  -\> see vcf_bind()

## Value

VCFArrow object
