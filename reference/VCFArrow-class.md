# VCFArrow S4 class

S4 class holding VCF data. Variant metadata are kept in memory while
genotypes are stored on disk as an Apache Arrow dataset and loaded
lazily. Objects are created with
[`read_vcf()`](https://legallab.github.io/VCFArrow/reference/read_vcf.md).

## Slots

- `header`:

  Complete VCF header (character).

- `info`:

  INFO field of each variant (character).

- `format`:

  FORMAT field of each variant (data.frame).

- `variants`:

  CHROM, POS, ID, REF, ALT, QUAL, FILTER and precalculated per-variant
  metrics (data.frame).

- `gt`:

  Arrow dataset of genotypes in long format.

- `samples`:

  Sample names (character).

- `groups`:

  Group assignment of each sample (character).

- `path`:

  Location of the on-disk Arrow dataset (character).

- `finalizer_env`:

  Environment used to clean up `path` on garbage collection.

- `loci`:

  CHROM, POS, REF and ALT of every variant of the source VCF, row i
  describing .row_id i; never filtered, so variants can be identified
  (e.g. by vcf_bind()) even after they were removed from `variants`
  (data.frame).

- `invariant_removed`:

  Variants removed as invariant, with their metadata (data.frame).

## Author

Tomas Hrbek April 2026
