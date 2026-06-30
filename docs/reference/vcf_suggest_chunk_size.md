# Suggest a chunk_size for read_vcf() given available RAM

Suggest a chunk_size for read_vcf() given available RAM

## Usage

``` r
vcf_suggest_chunk_size(available_gb, n_samples, n_columns = 5L)
```

## Arguments

- available_gb:

  RAM available in gigabytes.

- n_samples:

  Number of samples.

- n_columns:

  Number of columns per gt row (default 5: row_id, sample, a1, a2,
  phased).
